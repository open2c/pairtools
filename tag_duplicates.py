import sys, argparse, pipes

parser = argparse.ArgumentParser(
    'Tags every line of a pairsam with a duplicate tag'
    )
parser.add_argument('infile', nargs='?', 
        type=argparse.FileType('r'), 
        default=sys.stdin)
parser.add_argument("--comment-char", type=str, default="#", help="The first character of comment lines")

args = parser.parse_args()

COMMENT_CHAR = args.comment_char
INSTREAM = args.infile
OUTSTREAM = sys.stdout

for line in INSTREAM.readlines():
    if line.startswith(COMMENT_CHAR):
        OUTSTREAM.write(line)
        continue
    else:
        cols = line[:-1].split('\v')
        cols[0] = 'DD'
        
        for i in range(7, len(cols)):
            sam = cols[i]
            samcols = sam.split('\t')
            samcols[1] = str(int(samcols[1]) | 1024)
            for j in range(11, len(samcols)):
                if samcols[j].startswith('Yt:Z:'):
                    samcols[j] = 'Yt:Z:DD'
            cols[i] = '\t'.join(samcols)


        OUTSTREAM.write('\v'.join(cols))
        OUTSTREAM.write('\n')
