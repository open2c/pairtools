import sys
import argparse
from _distiller_common import open_bgzip

def main():
    parser = argparse.ArgumentParser(
        description='Read a pairsam file and print only the pairs of a certain type(s).'
    )
    parser.add_argument(
        'field',
        nargs=1,
        type=str, 
        choices=['pair_type', 'chrom1', 'chrom2', 'read_id'],
        metavar='field',
        help='The field to filter pairs by. Possible choices are: '
        'pair_type,chrom1,chrom2,read_id')
    parser.add_argument(
        'value', 
        nargs=1,
        type=str, 
        help='Select reads with `field` matching `value`. '
            'Depending on `--match-method`, this argument can be interpreted as '
            'a single value, a comma separated list, a wildcard or a regexp.')
    parser.add_argument(
        '--match-method', 
        choices=['comma_list', 'single_value', 'wildcard', 'regexp'],
        default='comma_list',
        metavar='',
        help='The method of matching of a field to `value`. Possible choices '
            'are: comma_list (`value` is interpreted as a comma-separated list of accepted values), '
                'single_value, wildcard, regexp (a Python-flavor regular expression).'
                )
    
    parser.add_argument(
        '--input',
        type=str, 
        default="",
        help='input pairsam file.'
            ' If the path ends with .gz, the input is gzip-decompressed.'
            ' By default, the input is read from stdin.')
    parser.add_argument(
        "--output", 
        type=str, 
        default="", 
        help='output file.'
            ' If the path ends with .gz, the output is bgzip-compressed.'
            ' By default, the output is printed into stdout.')
    parser.add_argument(
        "--output-rest", 
        type=str, 
        default="", 
        help='output file for pairs of other types. '
            ' If the path ends with .gz, the output is bgzip-compressed.'
            ' By default, such pairs are dropped.'
    )
    parser.add_argument(
        "--send-comments-to", 
        type=str, 
        default="both", 
        choices=['selected', 'rest', 'both', 'none'], 
        help="Which of the outputs should receive header and comment lines")
    
    args = vars(parser.parse_args())
    
    instream = (open_bgzip(args['input'], mode='r') 
                if args['input'] else sys.stdin)
    outstream = (open_bgzip(args['output'], mode='w') 
                 if args['output'] else sys.stdout)
    outstream_rest = (open_bgzip(args['output_rest'], mode='w') 
                      if args['output_rest'] else None)

    colidx = {
        'chrom1':0,
        'chrom1':3,
        'read_id':6,
        'pair_type':7,
        }[args['field'][0]]

    val = args['value'][0]

    if args['match_method'] == 'single_value':
        do_match = lambda x: x==val
    elif args['match_method'] == 'comma_list':
        vals = val.split(',')
        do_match = lambda x: any((i==x for i in vals))
    elif args['match_method'] == 'wildcard':
        import fnmatch, re
        regex = fnmatch.translate(val)
        reobj = re.compile(regex)
        do_match = lambda x: bool(reobj.match(x))
    elif args['match_method'] == 'regexp':
        import re
        reobj = re.compile(val)
        do_match = lambda x: bool(reobj.match(x))
    else:
        raise Exception('An unknown matching method: {}'.format(args['match_method']))

    for line in instream.readlines():
        if line.startswith('#'):

            outstream.write(line)
            if outstream_rest:
                outstream_rest.write(line)
        else:
            cols = line.split('\v')
            if do_match(cols[colidx]):
                outstream.write(line)
            elif outstream_rest:
                outstream_rest.write(line)

    if hasattr(instream, 'close'):
        instream.close()
    if hasattr(outstream, 'close'):
        outstream.close()
    if outstream_rest:
        outstream_rest.close()

if __name__ == '__main__':
    main()
