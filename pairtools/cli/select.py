import sys
import click
import re, fnmatch
import warnings

from ..lib import fileio, pairsam_format, headerops
from ..lib.select import evaluate_stream
from . import cli, common_io_options

UTIL_NAME = "pairtools_select"


@cli.command()
@click.argument("condition", type=str)
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output file."
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, the output is printed into stdout.",
)
@click.option(
    "--output-rest",
    type=str,
    default="",
    help="output file for pairs of other types. "
    " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
    " By default, such pairs are dropped.",
)

# Deprecated option to be removed in the future:
# @click.option(
#     "--send-comments-to",
#     type=click.Choice(['selected', 'rest', 'both', 'none']),
#     default="both",
#     help="Which of the outputs should receive header and comment lines",
#     show_default=True)


@click.option(
    "--chrom-subset",
    type=str,
    default=None,
    help="A path to a chromosomes file (tab-separated, 1st column contains "
    "chromosome names) containing a chromosome subset of interest. "
    "If provided, additionally filter pairs with both sides originating from "
    "the provided subset of chromosomes. This operation modifies the #chromosomes: "
    "and #chromsize: header fields accordingly.",
)
@click.option(
    "--startup-code",
    type=str,
    default=None,
    help="An auxiliary code to execute before filtering. "
    "Use to define functions that can be evaluated in the CONDITION statement",
)
@click.option(
    "-t",
    "--type-cast",
    type=(str, str),
    default=(),
    multiple=True,
    help="Cast a given column to a given type. By default, only pos and mapq "
    "are cast to int, other columns are kept as str. Provide as "
    "-t <column_name> <type>, e.g. -t read_len1 int. Multiple entries are allowed.",
)
@click.option(
    "--remove-columns",
    "-r",
    help=f"Comma-separated list of columns to be removed, e.g.: {','.join(pairsam_format.COLUMNS)}",
    type=str,
    default="",
    required=False,
)
@common_io_options
def select(
    condition,
    pairs_path,
    output,
    output_rest,  # send_comments_to,
    chrom_subset,
    startup_code,
    type_cast,
    remove_columns,
    **kwargs,
):
    """Select pairs according to some condition.

    CONDITION : A Python expression; if it returns True, select the read pair.
    Any column declared in the #columns line of the pairs header can be
    accessed by its name. If the header lacks the #columns line, the columns
    are assumed to follow the .pairs/.pairsam standard (readID, chrom1, chrom2,
    pos1, pos2, strand1, strand2, pair_type). Finally, CONDITION has access to
    COLS list which contains the string values of columns. In Bash, quote
    CONDITION with single quotes, and use double quotes for string variables
    inside CONDITION.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip/lz4c. By default, the input is read from stdin.

    The following functions can be used in CONDITION besides the standard Python functions:

    - csv_match(x, csv) - True if variable x is contained in a list of
    comma-separated values, e.g. csv_match(chrom1, 'chr1,chr2')

    - wildcard_match(x, wildcard) - True if variable x matches a wildcard,
    e.g. wildcard_match(pair_type, 'C*')

    - regex_match(x, regex) - True if variable x matches a Python-flavor regex,
    e.g. regex_match(chrom1, 'chr\d')

    \b
    Examples:
    pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")'
    pairtools select 'chrom1==chrom2'
    pairtools select 'COLS[1]==COLS[3]'
    pairtools select '(chrom1==chrom2) and (abs(pos1 - pos2) < 1e6)'
    pairtools select '(chrom1=="!") and (chrom2!="!")'
    pairtools select 'regex_match(chrom1, "chr\d+") and regex_match(chrom2, "chr\d+")'

    pairtools select 'True' --chrom-subset mm9.reduced.chromsizes

    """
    select_py(
        condition,
        pairs_path,
        output,
        output_rest,  # send_comments_to,
        chrom_subset,
        startup_code,
        type_cast,
        remove_columns,
        **kwargs,
    )


def select_py(
    condition,
    pairs_path,
    output,
    output_rest,  # send_comments_to,
    chrom_subset,
    startup_code,
    type_cast,
    remove_columns,
    **kwargs,
):

    instream = fileio.auto_open(
        pairs_path,
        mode="r",
        nproc=kwargs.get("nproc_in"),
        command=kwargs.get("cmd_in", None),
    )
    outstream = fileio.auto_open(
        output,
        mode="w",
        nproc=kwargs.get("nproc_out"),
        command=kwargs.get("cmd_out", None),
    )

    # Optional output created only if requested:
    outstream_rest = None
    if output_rest:
        outstream_rest = fileio.auto_open(
            output_rest,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )

    # Parse the input stream:
    header, body_stream = headerops.get_header(instream)

    # Modify the header:
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)

    # Filter out unwanted columns:
    if remove_columns:
        input_columns = headerops.extract_column_names(header)
        remove_columns = remove_columns.split(",")
        for col in remove_columns:
            if col in pairsam_format.COLUMNS_PAIRS:
                warnings.warn(
                    f"Removing required {col} column for .pairs format. Output is not .pairs anymore"
                )
            elif col in pairsam_format.COLUMNS_PAIRSAM:
                warnings.warn(
                    f"Removing required {col} column for .pairsam format. Output is not .pairsam anymore"
                )
        updated_columns = [x for x in input_columns if x not in remove_columns]

        if len(updated_columns) == len(input_columns):
            warnings.warn(
                f"Some column(s) {','.join(remove_columns)} not in the file, the operation has no effect"
            )
        else:
            header = headerops.set_columns(header, updated_columns)

    # Update the chromosomes:
    new_chroms = None
    if chrom_subset is not None:
        new_chroms = [l.strip().split("\t")[0] for l in open(chrom_subset, "r")]

    if new_chroms is not None:
        header = headerops.subset_chroms_in_pairsheader(header, new_chroms)
    outstream.writelines((l + "\n" for l in header))
    if output_rest:
        outstream_rest.writelines((l + "\n" for l in header))

    column_names = headerops.extract_column_names(header)
    if len(column_names) == 0:
        column_names = pairsam_format.COLUMNS

    # Columns filtration rule:
    if remove_columns:
        column_scheme = [input_columns.index(COL) for COL in updated_columns]

    # Format the condition:
    condition = condition.strip()
    if new_chroms is not None:
        condition = (
            f"({condition}) and (chrom1 in {new_chroms}) and (chrom2 in {new_chroms})"
        )

    for filter_passed, line in evaluate_stream(
        body_stream, condition, column_names, type_cast, startup_code
    ):
        COLS = line.rstrip('\n').split(pairsam_format.PAIRSAM_SEP)

        if remove_columns:
            COLS = [
                COLS[idx] for idx in column_scheme
            ]  # re-order the columns according to the scheme:
            line = pairsam_format.PAIRSAM_SEP.join(COLS) + "\n"  # form the line

        if filter_passed:
            outstream.write(line)
        elif outstream_rest:
            outstream_rest.write(line)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()

    if output_rest and outstream_rest != sys.stdout:
        outstream_rest.close()


if __name__ == "__main__":
    select()
