import sys, os, subprocess
import fileinput                                                      
import itertools                                                      
import io

# inline

ORD_M = ord('M')
ORD_S = ord('S')
ORD_H = ord('H')

def parse_cigar(CIGAR, strand):
    mapped_len = 0
    clip5 = 0
    clip3 = 0
    cur_num = 0
    total_len = 0
    for charval in CIGAR:
        #charval = ord(symbol)
        if charval >= 48 and charval <= 57:
            cur_num = cur_num * 10 + charval
        else:
            if charval == ORD_M: 
                mapped_len += cur_num
            elif charval == ORD_S:
                if mapped_len == 0:
                    clip5 = cur_num
                else: 
                    clip3 = cur_num
            elif charval == ORD_H:
                raise Exception('Unexpected hard-clipping!')

            total_len += cur_num
            cur_num = 0

    if strand == b'-':
        clip5, clip3 = clip3, clip5
                
    return clip5, clip3, mapped_len, total_len

# inline
def get_supp_alignment(samcols):
    for col in samcols:                                              
        if not col.startswith(b'SA:Z:'):                                     
            continue

        SAcols = col[5:].split(b',')
        chrom = SAcols[0]
        strand = SAcols[2]
        mapq = int(SAcols[4])

        clip5, clip3, mapped_len, total_len = parse_cigar(SAcols[3], strand)

        if strand == b'+':
            pos = int(SAcols[1]) + clip5
        else:
            pos = int(SAcols[1]) + total_len - clip5

        return {
            'chrom':chrom, 
            'pos':pos, 
            'strand':strand, 
            'mapq':mapq, 
            'clip5':clip5, 
            'clip3':clip3, 
            'mapped_len':mapped_len, 
            'total_len':total_len, 
            'is_mapped':True
            }
    return None

# inline
def get_algn_loc_mapq(samcols):                                         
    chrom = samcols[2]                                             
    strand = b'+' if ((int(samcols[1]) & 0x10) == 0) else b'-'
    mapq = int(samcols[4])

    is_mapped = (int(samcols[1]) & 0x04) == 0                       
    if is_mapped:
        clip5, clip3, mapped_len, total_len = parse_cigar(samcols[5], strand)
    else:
        clip5, clip3, mapped_len, total_len = 0, 0, 0, 0

    if strand == b'+':
        pos = int(samcols[3]) + clip5
    else:
        pos = int(samcols[3]) + total_len - clip5

    return {
        'chrom':chrom, 
        'pos':pos, 
        'strand':strand, 
        'mapq':mapq, 
        'clip5':clip5, 
        'clip3':clip3, 
        'mapped_len':mapped_len, 
        'total_len':total_len, 
        'is_mapped':is_mapped
        }

# inline in cython!
def get_chimera_position(
    repr_algn1, 
    repr_algn2, 
    sup_algn1, 
    sup_algn2, 
    min_mapq,
    max_chimera_dist):

    if (sup_algn1 is None) and (sup_algn2 is None):
        return repr_algn1, repr_algn2, False, True

    if (sup_algn1 is not None) and (sup_algn2 is not None):
        return None, None, True, False
        
    sup_algn = sup_algn1 if (sup_algn1 is not None) else sup_algn2
    # in normal chimeras, the non-unique alignment can be mapped non-uniquely
    if (sup_algn['mapq'] <= min_mapq):
        return repr_algn1, repr_algn2, True, True

    first_read_is_chimeric = (sup_algn2 is None)
    linear_algn = repr_algn2 if first_read_is_chimeric else repr_algn2
    repr_algn = repr_algn1 if first_read_is_chimeric else repr_algn2
    chim5_algn = repr_algn if (repr_algn['clip5'] < sup_algn['clip5']) else sup_algn
    chim3_algn = sup_algn  if (repr_algn['clip5'] < sup_algn['clip5']) else repr_algn

    is_normal = True
    # in normal chimeras, the supplemental alignment must be on the same chromosome with the linear alignment
    is_normal &= (chim3_algn['chrom'] == linear_algn['chrom']) 
    # in normal chimeras, the supplemental alignment must point along the linear alignment
    is_normal &= (chim3_algn['strand'] != linear_algn['strand'])
    
    # in normal chimeras, the supplemental alignment must be within
    # MAX_CHIMERA_DIST downstream of the linear alignment
    is_normal &= not (
        (linear_algn['strand']==b'+') and (
            (chim3_algn['pos'] - linear_algn['pos'] > max_chimera_dist)
            or (chim3_algn['pos'] - linear_algn['pos'] < 0)
        )
    )

    is_normal &= not (
        (linear_algn['strand']==b'-') and (
            (linear_algn['pos'] - chim3_algn['pos'] > max_chimera_dist)
            or (linear_algn['pos'] - chim3_algn['pos'] < 0)
        )
    )

    if is_normal:
        if first_read_is_chimeric:
            return chim5_algn, linear_algn, True, True
        else:
            return linear_algn, chim5_algn, True, True
    else:
        return None, None, True, False

def get_pair_order(chrm1, pos1, chrm2, pos2):
    if (chrm1 < chrm2):
        return -1
    elif (chrm1 > chrm2):
        return 1
    else:
        return int(pos1 < pos2) * 2 - 1

def write_pair_sam(algn1, algn2, sams1, sams2):
    # SAM is already tab-separated and
    # any printable character between ! and ~ may appear in the PHRED field! 
    # (http://www.ascii-code.com/)
    # Thus, use the vertical tab character to separate fields!

    sys.stdout.buffer.write(algn1['chrom'])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(str(algn1['pos']).encode())
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(algn1['strand'])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(algn2['chrom'])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(str(algn2['pos']).encode())
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(algn2['strand'])
    for sam in sams1:
        sys.stdout.buffer.write(b'\v')
        sys.stdout.buffer.write(sam[:-1])
    for sam in sams2:
        sys.stdout.buffer.write(b'\v')
        sys.stdout.buffer.write(sam[:-1])
    sys.stdout.buffer.write(b'\n')

def save_sams(out_file, sams1, sams2):
    for sam in sams1:
        out_file.write(sam)
    for sam in sams2:
        out_file.write(sam)

def process_sams(
    sams1, 
    sams2,
    unmapped_file, 
    singlesided_file, 
    multimapped_file, 
    abnormal_chimera_file):

    sam1_repr_cols = sams1[0].rstrip().split(b'\t')                         
    sam2_repr_cols = sams2[0].rstrip().split(b'\t')                         

    algn1 = get_algn_loc_mapq(sam1_repr_cols)
    algn2 = get_algn_loc_mapq(sam2_repr_cols)

    if (not algn1['is_mapped']) and (not algn2['is_mapped']):
        save_sams(unmapped_file, sams1, sams2)
        return

    elif (not algn1['is_mapped']) or (not algn2['is_mapped']):
        save_sams(singlesided_file, sams1, sams2)
        return

    elif (algn1['mapq'] < MIN_MAPQ) or (algn1['mapq'] < MIN_MAPQ):
        save_sams(multimapped_file, sams1, sams2)
        return

    sup_algn1 = get_supp_alignment(sam1_repr_cols)
    sup_algn2 = get_supp_alignment(sam2_repr_cols)

    algn1_5, algn2_5, is_chimera, is_normal_chimera = get_chimera_position(
        algn1, algn2, sup_algn1, sup_algn2, 
        MIN_MAPQ, MAX_CHIMERA_DIST)

    if is_chimera and (not is_normal_chimera):
        save_sams(abnormal_chimera_file, sams1, sams2)
        return

    pair_order = get_pair_order(
        algn1_5['chrom'], algn1_5['pos'], 
        algn2_5['chrom'], algn2_5['pos'])

    if pair_order >= 0:
        write_pair_sam(algn1_5, algn2_5, sams1, sams2)
    else:
        write_pair_sam(algn2_5, algn1_5, sams2, sams1)
    return


def append_sam(line, sams1, sam2):
    _, flag, _ = last_line.split(b'\t',2)
    flag = int(flag)

    if ((flag & 0x40) != 0):
        if ((flag & 0x800) == 0):
            sams1.insert(0,last_line)
        else:
            sams1.append(last_line)
    else:
        if ((flag & 0x800) == 0):
            sams2.insert(0,last_line)
        else:
            sams2.append(last_line)
    return 

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('Splits .sam entries into different'
    'read pair categories')
    parser.add_argument('infile', nargs='?', 
        type=argparse.FileType('rb'), 
        default=sys.stdin.buffer)
    parser.add_argument("--header", type=str, required=True)
    parser.add_argument("--unmapped", type=str, required=True)
    parser.add_argument("--singlesided", type=str, required=True)
    parser.add_argument("--multimapped", type=str, required=True)
    parser.add_argument("--abnormal-chimera", type=str, required=True)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--max-chimera-dist", type=int, default=10)

    args = parser.parse_args()
    MIN_MAPQ = args.min_mapq
    MAX_CHIMERA_DIST = args.max_chimera_dist
    IN_STREAM = args.infile

    header_file = open(args.header, 'wb')
    unmapped_file = open(args.unmapped, 'wb')
    singlesided_file = open(args.singlesided, 'wb')
    multimapped_file = open(args.multimapped, 'wb')
    abnormal_chimera_file = open(args.abnormal_chimera, 'wb')

    # The first few lines of a SAM file are header lines. These lines
    # begin with @ and they must be sent to all output files.
    last_line = IN_STREAM.readline()
    while last_line.startswith(b'@'):
        header_file.write(last_line)
        unmapped_file.write(last_line)
        singlesided_file.write(last_line)
        multimapped_file.write(last_line)
        abnormal_chimera_file.write(last_line)
        last_line = IN_STREAM.readline()

    header_file.flush()
    header_file.close()
    del(header_file)

    prev_read_id = b''
    sams1 = []
    sams2 = []
    while last_line:
        read_id, _ = last_line.split(b'\t', 1)
        if (read_id != prev_read_id) and (prev_read_id):
            process_sams(
                sams1, 
                sams2,
                unmapped_file, 
                singlesided_file, 
                multimapped_file, 
                abnormal_chimera_file)
            sams1.clear()
            sams2.clear()

        append_sam(last_line, sams1, sams2)
        prev_read_id = read_id
        last_line = IN_STREAM.readline()

    unmapped_file.close()
    singlesided_file.close()
    multimapped_file.close()
    abnormal_chimera_file.close()
                                                                                   
