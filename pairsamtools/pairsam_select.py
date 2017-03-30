import sys
import click
import re, fnmatch

from . import _common, cli, _headerops

UTIL_NAME = 'pairsam_select'

@cli.command()
@click.argument(
    'condition',
    type=str
)

@click.argument(
    'pairsam_path', 
    type=str,
    required=False)


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
    condition, pairsam_path, output, output_rest, send_comments_to
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

    wildcard_library = {}
    def wildcard_match(x, wildcard):
        if wildcard not in wildcard_library:
            regex = fnmatch.translate(wildcard)
            reobj = re.compile(regex)
            wildcard_library[wildcard] = reobj
        return wildcard_library[wildcard].match(x)

    csv_library = {}
    def csv_match(x, csv):
        if csv not in csv_library:
            csv_library[csv] = set(csv.split(','))
        return x in csv_library[csv]

    regex_library = {}
    def regex_match(x, regex):
        if regex not in regex_library:
            reobj = re.compile(regex)
            regex_library[regex] = reobj
        return regex_library[regex].match(x)

    match_func = compile(condition, '<string>', 'eval')

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l+'\n' for l in header))
    if outstream_rest:
        outstream_rest.writelines((l+'\n' for l in header))


    for line in body_stream:
        cols = line.split(_common.PAIRSAM_SEP)
        chrom1 = cols[_common.COL_C1]
        chrom2 = cols[_common.COL_C2]
        pos1 = int(cols[_common.COL_P1])
        pos2 = int(cols[_common.COL_P2])
        strand1 = cols[_common.COL_S1]
        strand2 = cols[_common.COL_S2]
        read_id = cols[_common.COL_READID]
        pair_type = cols[_common.COL_PTYPE]
        if eval(match_func):
#        if True:
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
