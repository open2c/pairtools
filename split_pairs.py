import argparse

parser = argparse.ArgumentParser('Splits .sam entries into different'
'read pair categories')
parser.add_argument("--input", type=str, default=None)
parser.add_argument("--header", type=str, default=None)
parser.add_argument("--out-pairs", type=str, required=True)
parser.add_argument("--out-sam", type=str, required=True)

args = parser.parse_args()

pairs_file = open(args.out_pairs, 'wb')
sam_file = open(args.out_sam, 'wb')

if args.header:
    with open(args.header, 'wb') as header_file:
        for line in header_file.readlines():
            pairs_file.write(line)
            sam_file.write(line)

if args.input:
    IN_STREAM = open(args.input, 'rb')
else:
    IN_STREAM = fileinput.input(mode='rb')

for line in IN_STREAM.readlines():
    if line.startswith(b'#'):
        pairs_file.write(line)
        continue

    cols = line.split(b'\v')
    pairs_file.write(b'\t'.join(cols[:6]))
    pairs_file.write(b'\n')
    
    sam_file.write(cols[6])
    sam_file.write(b'\n')
    sam_file.write(cols[7])

    
