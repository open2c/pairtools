import sys, subprocess

def open_bgzip(path, mode):
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


headers = []
for path in sys.argv[1:]:
    header = []
    f = open_bgzip(path, mode='r')
    for line in f.readlines():
        if line and not line.isspace():
            if line.strip().startswith('#'):
                header.append(line.strip())
            else:
                break
    f.close()
    headers.append(header)

SQ_headers = [set(line for line in header if line.startswith('#@SQ'))
              for header in headers]
common_sq_header = set.intersection(*SQ_headers)
sq_headers_same = all([len(header) == len(common_sq_header) 
                       for header in SQ_headers])

if not sq_headers_same:
    raise Exception('The SQ (sequence) lines of the sam headers are not identical')

other_headers = []
for header in headers:
    for line in header:
        line = line.strip()
        if line.startswith('#@') and not line.startswith('#@SQ'):
            if line not in other_headers:
                other_headers.append(line)

for header in headers:
    for line in header:
        line = line.strip()
        if line.startswith('#') and not line.startswith('#@'):
            if line not in other_headers:
                other_headers.append(line)
 
for line in headers[0]:
    if line.startswith('#@SQ'):
        print(line, flush=True)
for line in other_headers:
    print(line,flush=True)

command = r'''/bin/bash -c "sort -k 1,1 -k 4,4 -k 2,2n -k 5,5n -k 8,8 --field-separator=$'\v' --merge'''
for path in sys.argv[1:]:
    if path.endswith('.gz'):
        command += r' <(zcat {} | sed -n -e "/^[^#]/,$p")'.format(path)
    else:
        command += r' <(sed -n -e "/^[^#]/,$p" {})'.format(path)
command += r'"'
subprocess.call(command, shell=True)

