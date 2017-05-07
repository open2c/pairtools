import sys
import click
import re, fnmatch

from . import _io, _pairsam_format, cli, _headerops

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

    CONDITION : A Python expression; if it returns True, select the read pair.
    The column values are available as the following variables: READ_ID, 
    CHROM_1, CHROM_2, POS_1, POS_2, STRAND_1, STRAND_2 (and COLS array to access
    the string values of columns by index). In Bash, quote CONDITION with single
    quotes, and use double quotes for string variables inside CONDITION.

    PAIRSAM_PATH : input .pairsam file. If the path ends with .gz, the input is
    gzip-decompressed. By default, the input is read from stdin.

    The following functions can be used in CONDITION besides the standard Python functions:

    - csv_match(x, csv) - True if variable x is contained in a list of
    comma-separated values, e.g. csv_match(CHROM_1, 'chr1,chr2')

    - wildcard_match(x, wildcard) - True if variable x matches a wildcard,
    e.g. wildcard_match(PAIR_TYPE, 'C*')

    - regex_match(x, regex) - True if variable x matches a Python-flavor regex,
    e.g. regex_match(CHROM_1, 'chr\d')

    \b
    Examples:
    pairsam select '(PAIR_TYPE=="LL") or (PAIR_TYPE=="CC")'
    pairsam select 'CHROM_1==CHROM_2'
    pairsam select 'COLS[1]==COLS[3]'
    pairsam select '(CHROM_1==CHROM_2) and (abs(POS_1 - POS_2) < 1e6)'
    pairsam select '(CHROM_1=="!") and (CHROM_2!="!")'
    pairsam select 'regex_match(CHROM_1, "chr\d+") and regex_match(CHROM_2, "chr\d+")'

    '''
    select_py(
        condition, pairsam_path, output, output_rest, send_comments_to
    )
    
def select_py(
    condition, pairsam_path, output, output_rest, send_comments_to
    ):
    instream = (_io.open_bgzip(pairsam_path, mode='r') 
                if pairsam_path else sys.stdin)
    outstream = (_io.open_bgzip(output, mode='w') 
                 if output else sys.stdout)
    outstream_rest = (_io.open_bgzip(output_rest, mode='w') 
                      if output_rest else None)

    wildcard_library = {}
    def wildcard_match(x, wildcard):
        if wildcard not in wildcard_library:
            regex = fnmatch.translate(wildcard)
            reobj = re.compile(regex)
            wildcard_library[wildcard] = reobj
        return wildcard_library[wildcard].fullmatch(x)

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
        return regex_library[regex].fullmatch(x)
    
    condition = condition.strip()
    condition = condition.replace('PAIR_TYPE', 'COLS[_pairsam_format.COL_PTYPE]')
    condition = condition.replace('READ_ID', 'COLS[_pairsam_format.COL_READID]')
    condition = condition.replace('CHROM_1', 'COLS[_pairsam_format.COL_C1]')
    condition = condition.replace('CHROM_2', 'COLS[_pairsam_format.COL_C2]')
    condition = condition.replace('POS_1', 'int(COLS[_pairsam_format.COL_P1])')
    condition = condition.replace('POS_2', 'int(COLS[_pairsam_format.COL_P2])')
    condition = condition.replace('STRAND_1', 'COLS[_pairsam_format.COL_P1]')
    condition = condition.replace('STRAND_2', 'COLS[_pairsam_format.COL_P2]')
    match_func = compile(condition, '<string>', 'eval')

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l+'\n' for l in header))
    if outstream_rest:
        outstream_rest.writelines((l+'\n' for l in header))

    for line in body_stream:
        COLS = line[:-1].split(_pairsam_format.PAIRSAM_SEP)
        if eval(match_func):
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
