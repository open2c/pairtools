import sys, argparse

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

parser = argparse.ArgumentParser(
    'Splits .sam entries into different '
    'read pair categories')
parser.add_argument('infile', nargs='?', 
        type=argparse.FileType('r'), 
        default=sys.stdin)
parser.add_argument("--header", type=str, default=None)
parser.add_argument("--out-pairs", type=str, required=True)
parser.add_argument("--out-sam", type=str, required=True)

args = parser.parse_args()

pairs_file = open(args.out_pairs, 'w')
sam_file = open_bamsam(args.out_sam, 'w')

IN_STREAM = args.infile

header_sent = not(args.header)

for line in IN_STREAM.readlines():
    if not header_sent:
        header_file = open_bamsam(args.header, 'r')
        for header_line in header_file.readlines():
            sam_file.write(header_line)
        header_file.close()
        header_sent = True

    if line.startswith('#'):
        pairs_file.write(line)
        continue

    cols = line[:-1].split('\v')
    pairs_file.write('\t'.join(cols[:6]))
    pairs_file.write('\n')
    
    for col in cols[6:]:
        sam_file.write(col)
        sam_file.write('\n')

    
