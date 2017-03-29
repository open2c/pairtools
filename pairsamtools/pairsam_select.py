import sys
import click

from . import _common, cli, __version__

UTIL_NAME = 'pairsam_select'

@cli.command()
@click.argument(
    'field',
    metavar='FIELD',
    type=click.Choice(['pair_type', 'chrom1', 'chrom2', 'read_id']),
)

@click.argument(
    'value', 
    metavar='VALUE',
)

@click.argument(
    'pairsam_path', 
    type=str,
    required=False)
@click.option(
    '--match-method', 
    type=click.Choice(['comma_list', 'single_value', 'wildcard', 'regexp']),
    default='comma_list',
    help='The method of matching of a field to `value`. Possible choices '
        'are: comma_list (`value` is interpreted as a comma-separated list of accepted values), '
        'single_value, wildcard, regexp (a Python-flavor regular expression).',
    show_default=True,
)


@click.option(
    '-o', "--output", 
    type=str, 
    default="", 
    help='output file.'
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, the output is printed into stdout.')

@click.option(
    "--output-rest", 
    type=str, 
    default="", 
    help='output file for pairs of other types. '
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, such pairs are dropped.')

@click.option(
    "--send-comments-to", 
    type=click.Choice(['selected', 'rest', 'both', 'none']),
    default="both", 
    help="Which of the outputs should receive header and comment lines",
    show_default=True)

def select(
    field, value, pairsam_path, match_method,output, output_rest, send_comments_to
    ):
    '''select pairsam entries.

    FIELD : The field to select pairs by. Possible choices are: pair_type,
    chrom1, chrom2, read_id.

    VALUE : Select reads with FIELD matching VALUE. Depending on 
    --match-method, this argument can be interpreted as a single value, 
    a comma separated list, a wildcard or a regexp.

    PAIRSAM_PATH : input .pairsam file. If the path ends with .gz, the input is
    gzip-decompressed. By default, the input is read from stdin.
    '''
    
    instream = (_common.open_bgzip(pairsam_path, mode='r') 
                if pairsam_path else sys.stdin)
    outstream = (_common.open_bgzip(output, mode='w') 
                 if output else sys.stdout)
    outstream_rest = (_common.open_bgzip(output_rest, mode='w') 
                      if output_rest else None)

    colidx = {
        'chrom1':_common.COL_C1,
        'chrom2':_common.COL_C2,
        'read_id':_common.COL_READID,
        'pair_type':_common.COL_PTYPE,
        }[field]

    if match_method == 'single_value':
        do_match = lambda x: x==value
    elif match_method == 'comma_list':
        vals = value.split(',')
        do_match = lambda x: any((i==x for i in vals))
    elif match_method == 'wildcard':
        import fnmatch, re
        regex = fnmatch.translate(value)
        reobj = re.compile(regex)
        do_match = lambda x: bool(reobj.match(x))
    elif match_method == 'regexp':
        import re
        reobj = re.compile(value)
        do_match = lambda x: bool(reobj.match(x))
    else:
        raise Exception('An unknown matching method: {}'.format(match_method))

    header, pairsam_body_stream = _common.get_header(instream)
    header = _common.append_pg_to_sam_header(
        header,
        {'ID': UTIL_NAME,
         'PN': UTIL_NAME,
         'VN': __version__,
         'CL': ' '.join(sys.argv)
         })

    outstream.writelines(header)
    if outstream_rest:
        outstream.writelines(header)

    for line in pairsam_body_stream:
        cols = line.split(_common.PAIRSAM_SEP)
        if do_match(cols[colidx]):
            outstream.write(line)
        elif outstream_rest:
            outstream_rest.write(line)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()

    if outstream_rest:
        outstream_rest.close()

if __name__ == '__main__':
    select()
