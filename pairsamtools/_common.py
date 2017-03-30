import pipes

COL_READID = 0
COL_C1 = 1
COL_C2 = 2
COL_P1 = 3
COL_P2 = 4
COL_S1 = 5
COL_S2 = 6
COL_PTYPE = 7
COL_SAM1 = 8
COL_SAM2 = 9

PAIRSAM_SEP = '\t'
PAIRSAM_SEP_ESCAPE = r'\t'
SAM_SEP = '\031'
SAM_SEP_ESCAPE = r'\031'
INTER_SAM_SEP = '\031NEXT_SAM\031'

def open_sam_or_bam(path, mode):
    '''Opens a file as a bam file is `path` ends with .bam, otherwise 
    opens it as a sam.
    '''
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


def open_bgzip(path, mode):
    '''Opens a file as a bgzip file is `path` ends with .bam, otherwise 
    opens it as a text.
    '''
    if mode not in ['r','w']:
        raise Exception("mode can be either 'r' or 'w'")
    if path.endswith('.gz'):
        if mode =='w': 
            t = pipes.Template()
            t.append('bgzip -c', '--')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append('zcat', '--')
            f = t.open(path, 'r')
        else:
            raise Exception("Unknown mode : {}".format(mode))
        return f
    else:
        return open(path, mode)


