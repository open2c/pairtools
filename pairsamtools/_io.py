import pipes

def open_sam_or_bam(path, mode, nproc=8):
    '''Opens a file as a bam file is `path` ends with .bam, otherwise 
    opens it as a sam.
    '''
    if mode not in ['r','w']:
        raise Exception("mode can be either 'r' or 'w'")
    if path.endswith('.bam'):
        if mode =='w': 
            t = pipes.Template()
            t.append('samtools view -bS {}'.format(
                         '-@ '+str(nproc-1) if nproc>1 else ''),
                     '--')
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


def open_bgzip(path, mode, nproc=8):
    '''Opens a file as a bgzip file is `path` ends with .bam, otherwise 
    opens it as a text.
    '''
    if mode not in ['r', 'w', 'a']:
        raise Exception("mode can be either 'r', 'w' or 'a'")
    if path.endswith('.gz'):
        if mode =='w': 
            t = pipes.Template()
            t.append('pbgzip -c -n {}'.format(nproc), '--')
            f = t.open(path, 'w')
        elif mode =='a': 
            t = pipes.Template()
            t.append('pbgzip -c -n {} $IN >> $OUT'.format(nproc), 'ff')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append('pbgzip -dc -n {}'.format(nproc), '--')
            f = t.open(path, 'r')
        else:
            raise Exception("Unknown mode : {}".format(mode))
        return f
    else:
        return open(path, mode)


