import sys, os, subprocess, pipes
import fileinput                                                      
import itertools                                                      
import io

# inline

ORD_M = ord('M')
ORD_I = ord('I')
ORD_D = ord('D')
ORD_S = ord('S')
ORD_H = ord('H')

def parse_cigar(CIGAR):
    matched_bp = 0
    algn_ref_span = 0
    algn_read_span = 0
    read_len = 0
    clip5 = 0
    clip3 = 0

    if CIGAR != '*':
        cur_num = 0
        for char in CIGAR:
            charval = ord(char)
            if charval >= 48 and charval <= 57:
                cur_num = cur_num * 10 + (charval - 48)
            else:
                if charval == ORD_M: 
                    matched_bp += cur_num
                    algn_ref_span += cur_num
                    algn_read_span += cur_num
                    read_len += cur_num
                elif charval == ORD_I: 
                    algn_read_span += cur_num
                    read_len += cur_num
                elif charval == ORD_D: 
                    algn_ref_span += cur_num
                elif charval == ORD_S:
                    read_len += cur_num
                    if matched_bp == 0:
                        clip5 = cur_num
                    else: 
                        clip3 = cur_num

                cur_num = 0
                    
    return {
        'clip5':clip5, 
        'clip3':clip3, 
        'algn_ref_span':algn_ref_span, 
        'algn_read_span':algn_read_span,
        'read_len':read_len,
        'matched_bp':matched_bp,
    }

# inline
def get_supp_alignment(samcols):
    for col in samcols:                                              
        if not col.startswith('SA:Z:'):                                     
            continue

        SAcols = col[5:].split(',')
        chrom = SAcols[0]
        strand = SAcols[2]
        mapq = int(SAcols[4])

        cigar = parse_cigar(SAcols[3])

        if strand == '+':
            pos = int(SAcols[1])
        else:
            pos = int(SAcols[1]) + cigar['algn_ref_span']

        return {
            'chrom':chrom, 
            'pos':pos, 
            'strand':strand, 
            'mapq':mapq, 
            'is_mapped':True,
            'cigar':cigar,
            'dist_to_5':cigar['clip5'] if strand == '+' else cigar['clip3'],
            }
    return None

# inline
def get_algn_loc_mapq(samcols):                                         
    chrom = samcols[2]                                             
    strand = '+' if ((int(samcols[1]) & 0x10) == 0) else '-'
    mapq = int(samcols[4])

    is_mapped = (int(samcols[1]) & 0x04) == 0                       
    cigar = parse_cigar(samcols[5])
    pos = 0
    if is_mapped:
        if strand == '+':
            pos = int(samcols[3])
        else:
            pos = int(samcols[3]) + cigar['algn_ref_span']

    return {
        'chrom':chrom, 
        'pos':pos, 
        'strand':strand, 
        'mapq':mapq, 
        'is_mapped':is_mapped,
        'cigar':cigar,
        'dist_to_5':cigar['clip5'] if strand == '+' else cigar['clip3'],
    }
    
# inline in cython!
def get_chimera_position(
    repr_algn1, 
    repr_algn2, 
    sup_algn1, 
    sup_algn2, 
    min_mapq,
    max_molecule_size):

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

    chim5_algn = repr_algn if (repr_algn['dist_to_5'] < sup_algn['dist_to_5']) else sup_algn
    chim3_algn = sup_algn  if (repr_algn['dist_to_5'] < sup_algn['dist_to_5']) else repr_algn

    is_normal = True
    # in normal chimeras, the supplemental alignment must be on the same chromosome with the linear alignment
    is_normal &= (chim3_algn['chrom'] == linear_algn['chrom']) 
    # in normal chimeras, the supplemental alignment must point along the linear alignment
    is_normal &= (chim3_algn['strand'] != linear_algn['strand'])
    
    # in normal chimeras, the supplemental alignment must be within
    # MAX_CHIMERA_DIST downstream of the linear alignment
    unmapped_gap_width = (
        (chim3_algn['pos'] - linear_algn['pos'] 
        - linear_algn['cigar']['read_len'] + linear_algn['dist_to_5'])
        if linear_algn['strand'] == '+' 
        else (
        chim3_algn['pos'] - linear_algn['pos'] 
        - chim3_algn['cigar']['read_len'] + chim3_algn['dist_to_5']
        )
    )

    is_normal &= (unmapped_gap_width > 0)
    is_normal &= (unmapped_gap_width
                  + linear_algn['cigar']['read_len'] 
                  + repr_algn['cigar']['read_len'] 
                  < max_molecule_size
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

def open_bamsam(path, mode):
    if mode not in ['r','w']:
        raise Exception("mode can be either 'r' or 'w'")
    if path.endswith('.bam'):
        if mode =='w': 
            t = pipes.Template()
            t.append('samtools view -bS', '--')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append('samtools view -h', '--')
            f = t.open(path, 'r')
        else:
            raise Exception("Unknown mode : {}".format(mode))
        return f
    else:
        return open(path, mode)

def write_pair_sam(algn1, algn2, sams1, sams2, out_file=sys.stdout):
    # SAM is already tab-separated and
    # any printable character between ! and ~ may appear in the PHRED field! 
    # (http://www.ascii-code.com/)
    # Thus, use the vertical tab character to separate fields!

    out_file.write(algn1['chrom'])
    out_file.write('\v')
    out_file.write(str(algn1['pos']))
    out_file.write('\v')
    out_file.write(algn1['strand'])
    out_file.write('\v')
    out_file.write(algn2['chrom'])
    out_file.write('\v')
    out_file.write(str(algn2['pos']))
    out_file.write('\v')
    out_file.write(algn2['strand'])
    for sam in sams1:
        out_file.write('\v')
        out_file.write(sam[:-1])
    for sam in sams2:
        out_file.write('\v')
        out_file.write(sam[:-1])
    out_file.write('\n')


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

    sam1_repr_cols = sams1[0].rstrip().split('\t')                         
    sam2_repr_cols = sams2[0].rstrip().split('\t')                         

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
        MIN_MAPQ, MAX_MOLECULE_SIZE)
    is_normal_chimera = is_normal_chimera or (len(sams1) > 2) or (len(sams2) >2)

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


def append_sams(line, sams1, sam2):
    _, flag, _ = last_line.split('\t',2)
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
        type=argparse.FileType('r'), 
        default=sys.stdin)
    parser.add_argument("--out-basename", type=str, default='')
    parser.add_argument("--out-extension", type=str, default='bam', choices=['bam','sam'])
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--max-molecule-size", type=int, default=2000)
    parser.add_argument("--header", type=str, default='')
    parser.add_argument("--unmapped", type=str, default='')
    parser.add_argument("--singlesided", type=str, default='')
    parser.add_argument("--multimapped", type=str, default='')
    parser.add_argument("--abnormal-chimera", type=str, default='')

    args = parser.parse_args()
    MIN_MAPQ = args.min_mapq
    MAX_MOLECULE_SIZE = args.max_molecule_size
    IN_STREAM = args.infile

    if not(args.out_basename) and not(
        args.header 
        and args.unmapped 
        and args.singlesided
        and args.multimapped
        and args.abnormal_chimera):
        raise Exception(
            'Please, provide either out_basename or individual '
            'output paths.')

    if (args.out_basename) and (
        args.header 
        or args.unmapped 
        or args.singlesided
        or args.multimapped
        or args.abnormal_chimera):
        raise Exception(
            'Please, provide only either out_basename or individual '
            'output paths.')

    if args.out_basename:
        header_file = open_bamsam(args.out_basename + '.header.' + args.out_extension, 'w')
        unmapped_file  = open_bamsam(args.out_basename + '.unmapped.' + args.out_extension, 'w')
        singlesided_file = open_bamsam(args.out_basename + '.singlesided.' + args.out_extension, 'w')
        multimapped_file = open_bamsam(args.out_basename + '.multimapped.' + args.out_extension, 'w')
        abnormal_chimera_file = open_bamsam(args.out_basename + '.abnormal_chimera.' + args.out_extension, 'w')
    else:
        header_file = open(args.header, 'w')
        unmapped_file = open(args.unmapped, 'w')
        singlesided_file = open(args.singlesided, 'w')
        multimapped_file = open(args.multimapped, 'w')
        abnormal_chimera_file = open(args.abnormal_chimera, 'w')

    # The first few lines of a SAM file are header lines. These lines
    # begin with @ and they must be sent to all output files.
    last_line = IN_STREAM.readline()
    while last_line.startswith('@'):
        header_file.write(last_line)
        unmapped_file.write(last_line)
        singlesided_file.write(last_line)
        multimapped_file.write(last_line)
        abnormal_chimera_file.write(last_line)
        last_line = IN_STREAM.readline()

    header_file.flush()
    header_file.close()
    del(header_file)

    prev_read_id = ''
    sams1 = []
    sams2 = []
    while last_line:
        read_id, _ = last_line.split('\t', 1)
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

        append_sams(last_line, sams1, sams2)
        prev_read_id = read_id
        last_line = IN_STREAM.readline()

    # taking care of the last pair
    process_sams(
        sams1, 
        sams2,
        unmapped_file, 
        singlesided_file, 
        multimapped_file, 
        abnormal_chimera_file)

    unmapped_file.close()
    singlesided_file.close()
    multimapped_file.close()
    abnormal_chimera_file.close()
                                                                                   
