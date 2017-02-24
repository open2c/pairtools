import sys, argparse

parser = argparse.ArgumentParser('Splits .sam entries into different'
'read pair categories')
parser.add_argument('infile', nargs='?', 
        type=argparse.FileType('rb'), 
        default=sys.stdin.buffer)
parser.add_argument("--header", type=str, default=None)
parser.add_argument("--out-pairs", type=str, required=True)
parser.add_argument("--out-sam", type=str, required=True)

args = parser.parse_args()

pairs_file = open(args.out_pairs, 'wb')
sam_file = open(args.out_sam, 'wb')

IN_STREAM = args.infile

header_sent = not(args.header)

for line in IN_STREAM.readlines():
    if not header_sent:
        with open(args.header, 'rb') as header_file:
            for header_line in header_file.readlines():
                sam_file.write(header_line)
        header_sent = True

    if line.startswith(b'#'):
        pairs_file.write(line)
        continue

    cols = line[:-1].split(b'\v')
    pairs_file.write(b'\t'.join(cols[:6]))
    pairs_file.write(b'\n')
    
    for col in cols:
        sam_file.write(col)
        sam_file.write(b'\n')

    
