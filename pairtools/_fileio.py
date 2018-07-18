import shutil
import pipes


class ParseError(Exception):
    pass

def auto_open(path, mode, nproc=1, command=None):
    '''Guess the file format from the extension and use the corresponding binary 
    to open it for reading or writing. If the extension is not known, open the
    file as text.

    If the binary allows parallel execution, specify the number of threads 
    with `nproc`.

    If `command` is supplied, use it to open the file instead of auto-guessing.
    The command must accept the filename as the last argument, accept input 
    through stdin and print output into stdout.

    Supported extensions and binaries (with comments):
    .bam - samtools view (allows parallel writing)
    .gz - pbgzip 
    .lz4 - lz4c (does not support parallel execution)
    '''
    if command:
        if mode =='w': 
            t = pipes.Template()
            t.append(command, '--')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append(command, '--')
            f = t.open(path, 'r')
        else:
            raise ValueError("Unknown mode : {}".format(mode))
        return f
    elif path.endswith('.bam'):
        if shutil.which('samtools') is None:
            raise ValueError({
                'w':'samtools is not found, cannot compress output',
                'r':'samtools is not found, cannot decompress input'
                    }[mode])
        if mode =='w': 
            t = pipes.Template()
            t.append('samtools view -bS {} -'.format(
                         '-@ '+str(nproc-1) if nproc>1 else ''),
                     '--')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append('samtools view -h', '--')
            f = t.open(path, 'r')
        else:
            raise ValueError("Unknown mode for .bam : {}".format(mode))
        return f
    elif path.endswith('.gz'):
        if shutil.which('pbgzip') is None:
            raise ValueError({
                'w':'pbgzip is not found, cannot compress output',
                'a':'pbgzip is not found, cannot compress output',
                'r':'pbgzip is not found, cannot decompress input'
                    }[mode])
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
            raise ValueError("Unknown mode for .gz : {}".format(mode))
        return f
    elif path.endswith('.lz4'):
        if shutil.which('lz4c') is None:
            raise ValueError({
                'w':'lz4c is not found, cannot compress output',
                'a':'lz4c is not found, cannot compress output',
                'r':'lz4c is not found, cannot decompress input'
                    }[mode])
        if mode =='w': 
            t = pipes.Template()
            t.append('lz4c -cz', '--')
            f = t.open(path, 'w')
        elif mode =='a': 
            t = pipes.Template()
            t.append('lz4c -cz $IN >> $OUT', 'ff')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append('lz4c -cd', '--')
            f = t.open(path, 'r')
        else:
            raise ValueError("Unknown mode : {}".format(mode))
        return f
    else:
        return open(path, mode)


