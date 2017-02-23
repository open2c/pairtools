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

        return (chrom, pos, strand, mapq, 
                clip5, clip3, mapped_len, total_len, True)
    return None


# inline
def get_algn_loc_mapq(samcols):                                         
    chrom = samcols[2]                                             
    strand = b'+' if ((int(samcols[1]) & 0x10) == 0) else b'-'
    mapq = int(samcols[4])

    ismapped = (int(samcols[1]) & 0x04) == 0                       
    if ismapped:
        clip5, clip3, mapped_len, total_len = parse_cigar(samcols[5], strand)
    else:
        clip5, clip3, mapped_len, total_len = 0, 0, 0, 0

    if strand == b'+':
        pos = int(samcols[3]) + clip5
    else:
        pos = int(samcols[3]) + total_len - clip5

    return (chrom, pos, strand, mapq, 
            clip5, clip3, mapped_len, total_len, ismapped)

# inline in cython!
def is_normal_chimera(
    repr_algn1, 
    repr_algn2, 
    sup_algn1, 
    sup_algn2, 
    min_mapq,
    max_chimera_dist):

    if (sup_algn1 is None) and (sup_algn2 is None):
        return True

    if (sup_algn1 is not None) and (sup_algn2 is not None):
        return False
        
    sup_algn = sup_algn1 if (sup_algn1 is not None) else sup_algn2
    # in normal chimeras, the non-unique alignment can be mapped non-uniquely
    if (sup_algn[3] <= min_mapq):
        return True

    linear_algn = repr_algn1 if (sup_algn1 is None) else repr_algn2
    repr_algn = repr_algn1 if (sup_algn1 is None) else repr_algn2
    chim5_algn = repr_algn if (repr_algn[4] < sup_algn[4]) else sup_algn
    chim3_algn = sup_algn  if (repr_algn[4] < sup_algn[4]) else repr_algn

    is_normal = True
    # in normal chimeras, the supplemental alignment must be on the same chromosome with the linear alignment
    is_normal &= (chim3_algn[0] == linear_algn[0]) 
    # in normal chimeras, the supplemental alignment must point along the linear alignment
    is_normal &= (chim3_algn[2] != linear_algn[2])
    
    # in normal chimeras, the supplemental alignment must be within
    # MAX_CHIMERA_DIST downstream of the linear alignment
    is_normal &= not (
        (linear_algn[2]==b'+') and (
            (chim3_algn[1] - linear_algn[1] > max_chimera_dist)
            or (chim3_algn[1] - linear_algn[1] < 0)
        )
    )

    is_normal &= not (
        (linear_algn[2]==b'-') and (
            (linear_algn[1] - chim3_algn[1] > max_chimera_dist)
            or (linear_algn[1] - chim3_algn[1] < 0)
        )
    )

    return is_normal

def cmp_chrms(chrm1, chrm2):
    if (chrm1 < chrm2):
        return -1
    elif (chrm1 > chrm2):
        return 1
    else:
        return 0

#def sam_bam_writer(out_path):                                                          
#    """                                                                            
#    Creates a writing process using samtools view -bS.
#    """    
#    if out_path.endswith('.bam'):
#        open(args.unmapped, 'wb')                                                                        
#        writer = ["samtools", "view", "-bS", "-"]                                          
#        pwrite = subprocess.Popen(
#            writer, 
#            stdin=subprocess.PIPE, 
#            stdout=out_stream, 
#            shell=False, 
#            bufsize=-1)
#    return pwrite
#
#def sam_writer(out_stream):                                                          
#    """                                                                            
#    Creates a writing process using samtools view -bS.
#    """                                                                            
#    writer = ["samtools", "view", "-bS", "-"]                                          
#    pwrite = subprocess.Popen(
#        writer, 
#        stdin=subprocess.PIPE, 
#        stdout=out_stream, 
#        shell=False, 
#        bufsize=-1)
#    return pwrite

def write_pair_sam(algn1, algn2, sam1, sam2):
    sys.stdout.buffer.write(algn1[0])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(str(algn1[1]).encode())
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(algn1[2])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(algn2[0])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(str(algn2[1]).encode())
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(algn2[2])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(sam1[:-1])
    sys.stdout.buffer.write(b'\v')
    sys.stdout.buffer.write(sam2[:-1])
    sys.stdout.buffer.write(b'\n')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('Splits .sam entries into different'
    'read pair categories')
    parser.add_argument("--input", type=str, default=None)
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

    if args.input:
        IN_STREAM = open(args.input, 'rb')
    else:
        IN_STREAM = fileinput.input(mode='rb')

    header_file = open(args.header, 'wb')
    header_writer = bam_writer(header_file)                                                
    header_pipe = header_writer.stdin   

    unmapped_file = open(args.unmapped, 'wb')
    unmapped_writer = bam_writer(unmapped_file)                                                
    unmapped_pipe = unmapped_writer.stdin   

    singlesided_file = open(args.singlesided, 'wb')
    singlesided_writer = bam_writer(singlesided_file)                                                
    singlesided_pipe = singlesided_writer.stdin   

    multimapped_file = open(args.multimapped, 'wb')
    multimapped_writer = bam_writer(multimapped_file)                                                
    multimapped_pipe = multimapped_writer.stdin   

    abnormal_chimera_file = open(args.abnormal_chimera, 'wb')
    abnormal_chimera_writer = bam_writer(abnormal_chimera_file)                                                
    abnormal_chimera_pipe = abnormal_chimera_writer.stdin   
    
    header = []

    while True:
        sam1 = IN_STREAM.readline()
        if not sam1:
            break

        if sam1.startswith(b'@'):
            header.append(sam1)

            unmapped_pipe.write(sam1)
            singlesided_pipe.write(sam1)
            multimapped_pipe.write(sam1)
            abnormal_chimera_pipe.write(sam1)
            continue

        sam2 = IN_STREAM.readline()

        samcols1 = sam1.rstrip().split(b'\t')                         
        samcols2 = sam2.rstrip().split(b'\t')                       

        algn1 = get_algn_loc_mapq(samcols1)
        algn2 = get_algn_loc_mapq(samcols2)

        if (not algn1[8]) and (not algn2[8]):
            unmapped_pipe.write(sam1)
            unmapped_pipe.write(sam2)
            continue

        elif (not algn1[8]) or (not algn2[8]):
            singlesided_pipe.write(sam1)
            singlesided_pipe.write(sam2)
            continue

        elif (algn1[3] < MIN_MAPQ) or (algn1[3] < MIN_MAPQ):
            multimapped_pipe.write(sam1)
            multimapped_pipe.write(sam2)
            continue

        sup_algn1 = get_supp_alignment(samcols1)
        sup_algn2 = get_supp_alignment(samcols2)

        if not is_normal_chimera(
            algn1, algn2, sup_algn1, sup_algn2, 
            MIN_MAPQ, MAX_CHIMERA_DIST):

            abnormal_chimera_pipe.write(sam1)
            abnormal_chimera_pipe.write(sam2)
            continue

        pair_order = cmp_chrms(algn1[0], algn2[0])
        if pair_order == 0:
            pair_order = int(algn1[1] > algn2[1]) * 2 - 1

        # SAM is already tab-separated and
        # any printable character between ! and ~ may appear in the PHRED field! 
        # (http://www.ascii-code.com/)
        # Thus, use the vertical tab character to separate fields!

        if pair_order == 1:
            write_pair_sam(algn1, algn2, sam1, sam2)

        if pair_order == -1:
            write_pair_sam(algn2, algn1, sam2, sam1)

    unmapped_writer.communicate()
    singlesided_writer.communicate()
    multimapped_writer.communicate()
    abnormal_chimera_writer.communicate()

    unmapped_file.close()
    singlesided_file.close()
    multimapped_file.close()
    abnormal_chimera_file.close()
                                                                                   
