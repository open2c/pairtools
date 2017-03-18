import sys
import argparse
from _distiller_common import open_bgzip

def main():
    parser = argparse.ArgumentParser(
      'Read a pairsam file and print only the pairs of a certain type(s).'
    )
    parser.add_argument(
        'pair_types', 
        nargs='+',
        type=str, 
        help='The list of requested pair types.')
    parser.add_argument(
        '--input',
        type=str, 
        default="",
        help='input file. By default, the input is read from stdin.')
    parser.add_argument(
        "--output", 
        type=str, 
        default="", 
        help="output file. By default, the output is printed into stdout.")
    parser.add_argument(
        "--output-rest", 
        type=str, 
        default="", 
        help="output file for pairs of other types. By default, such pairs are dropped."
        )
    
    args = vars(parser.parse_args())
    
    instream = (open_bgzip(args['input'], mode='r') 
                if args['input'] else sys.stdin)
    outstream = (open_bgzip(args['outstream'], mode='w') 
                 if args['output'] else sys.stdout)
    outstream_rest = (open_bgzip(args['outstream_rest'], mode='w') 
                      if args['output_rest'] else None)

    pair_types = args['pair_types']

    for line in instream.readlines():
        if line.startswith('#'):

            outstream.write(line)
            if outstream_rest:
                outstream_rest.write(line)
        else:
            cols = line.split('\v')
            if any((pair_type == cols[7] for pair_type in pair_types)):
                outstream.write(line)
            elif outstream_rest:
                outstream_rest.write(line)

if __name__ == '__main__':
    main()
