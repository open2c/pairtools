import sys
import click
import warnings
import subprocess

from ..lib import fileio, pairsam_format, headerops
from ..lib.parse_pysam import AlignmentFilePairtoolized
from . import cli, common_io_options


UTIL_NAME = "pairtools_header"


@cli.group()
def header():
    """
    Manipulate the .pairs/.pairsam header
    """
    pass


# Common options for all header tools:
def register_subcommand(func):
    return header.command()(
        click.argument("pairs_path", type=str, required=False)(
            click.option(
                "-o",
                "--output",
                type=str,
                default="",
                help="output file."
                " If the path ends with .gz or .lz4, the output is bgzip-/lz4c-compressed."
                " By default, the output is printed into stdout.",
            )(
                click.option(
                    "--nproc-in",
                    type=int,
                    default=1,
                    show_default=True,
                    help="Number of processes used by the auto-guessed input decompressing command.",
                )(
                    click.option(
                        "--nproc-out",
                        type=int,
                        default=8,
                        show_default=True,
                        help="Number of processes used by the auto-guessed output compressing command.",
                    )(
                        click.option(
                            "--cmd-in",
                            type=str,
                            default=None,
                            help="A command to decompress the input. "
                            "If provided, fully overrides the auto-guessed command. "
                            "Does not work with stdin. "
                            "Must read input from stdin and print output into stdout. "
                            "EXAMPLE: pbgzip -dc -n 3",
                        )(
                            click.option(
                                "--cmd-out",
                                type=str,
                                default=None,
                                help="A command to compress the output. "
                                "If provided, fully overrides the auto-guessed command. "
                                "Does not work with stdout. "
                                "Must read input from stdin and print output into stdout. "
                                "EXAMPLE: pbgzip -c -n 8",
                            )(func)
                        )
                    )
                )
            )
        )
    )


def add_arg_help(func):
    func.__doc__ = func.__doc__.format(
        """
    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by bgzip/lz4c. By default, the input is read from stdin.
    """
    )
    return func


@register_subcommand
@add_arg_help
@click.option(
    "--chroms-path",
    type=str,
    default=None,
    required=False,
    help="Chromosome order used to flip interchromosomal mates: "
    "path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose "
    "first column lists scaffold names. Any scaffolds not listed will be "
    "ordered lexicographically following the names provided.",
)
@click.option(
    "--sam-path",
    type=str,
    default=None,
    required=False,
    help="Input sam file to inherit the header."
    " Either --sam or --chroms-path should be provided to store the chromosome sizes in the header.",
)
@click.option(
    "--columns",
    type=click.STRING,
    default="",
    help="Report columns describing alignments "
    "Can take multiple values as a comma-separated list."
    f"By default, assign standard .pairs columns: {','.join(pairsam_format.COLUMNS)}",
)
@click.option(
    "--extra-columns",
    type=click.STRING,
    default="",
    help="Report extra columns describing alignments "
    "Can take multiple values as a comma-separated list.",
)
@click.option(
    "--assembly",
    type=str,
    default="",
    help="Name of genome assembly (e.g. hg19, mm10) to store in the pairs header.",
)
@click.option(
    "--no-flip",
    is_flag=True,
    help="If specified, assume that the pairs are not filpped in genomic order and instead preserve "
    "the order in which they were sequenced.",
)
@click.option(
    "--pairs/--pairsam",
    is_flag=True,
    default=True,
    help=f"If pairs, then the defult columns will be set to: {','.join(pairsam_format.COLUMNS_PAIRS)}"
    f"\nif pairsam, then to: {','.join(pairsam_format.COLUMNS_PAIRSAM)}",
)
def generate(pairs_path, output, chroms_path, sam_path, columns, assembly, **kwargs):
    """
    Generate the header

    """
    generate_py(pairs_path, output, chroms_path, sam_path, columns, assembly, **kwargs)


def generate_py(pairs_path, output, chroms_path, sam_path, columns, assembly, **kwargs):

    instream = fileio.auto_open(
        pairs_path,
        mode="r",
        nproc=kwargs.get("nproc_in"),
        command=kwargs.get("cmd_in", None),
    )
    header, body_stream = headerops.get_header(instream, ignore_warning=True)

    # Parse chromosome sizes present in the input chromosomes:
    if chroms_path and not sam_path:
        chromsizes = headerops.get_chromsizes_from_file(chroms_path)
        # chromosomes = headerops.get_chromsizes_from_file(chroms_path)

    # Parse chromosome sizes present in sam input:
    if sam_path:  # open input sam file with pysam
        input_sam = AlignmentFilePairtoolized(
            sam_path, "r", threads=kwargs.get("nproc_in")
        )
        samheader = input_sam.header
        chromsizes = headerops.get_chromsizes_from_pysam_header(samheader)
        # if chroms_path:
        #     chromosomes = headerops.get_chrom_order(chroms_path, list(chromsizes.keys()))
        # else:
        #     chromosomes = chromsizes.keys()

    # Read the input columns:
    if columns:
        columns = columns.split(",")
    else:
        if kwargs.get("pairs", True):
            columns = pairsam_format.COLUMNS_PAIRS
        else:
            columns = pairsam_format.COLUMNS_PAIRSAM

    extra_columns = kwargs.get("extra_columns", "")
    if extra_columns:
        columns += extra_columns.split(",")

    # Write new header to the pairsam file
    new_header = headerops.make_standard_pairsheader(
        assembly=assembly,
        chromsizes=chromsizes,
        columns=columns,
        shape="whole matrix" if kwargs["no_flip"] else "upper triangle",
    )

    if sam_path:
        new_header = headerops.insert_samheader_pysam(new_header, samheader)

    new_header = headerops.append_new_pg(new_header, ID=UTIL_NAME, PN=UTIL_NAME)

    # Check that the number of columns in the body corresponds to the header:
    if not headerops.validate_cols(instream, columns):
        raise ValueError(
            f"Number of columns mismatch:\n\t#columns: {headerops.SEP_COLS.join(columns)}\n\t{body_stream.readline()}"
        )

    ########
    # Write the output after successful checks:
    outstream = (
        fileio.auto_open(
            output,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output
        else sys.stdout
    )

    outstream.writelines((l + "\n" for l in new_header))
    outstream.flush()

    if body_stream == sys.stdin:
        for line in body_stream:
            outstream.write(line)
    else:
        command = r"""
                   /bin/bash -c 'export LC_COLLATE=C; export LANG=C; cat """

        if kwargs.get("cmd_in", None):
            command += r""" <(cat {} | {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path, kwargs["cmd_in"]
            )
        elif pairs_path.endswith(".gz"):
            command += (
                r""" <(bgzip -dc -@ {} {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                    kwargs["nproc_in"], pairs_path
                )
            )
        elif pairs_path.endswith(".lz4"):
            command += r""" <(lz4c -dc {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path
            )
        else:
            command += r""" <(sed -n -e '\''/^[^#]/,$p'\'' {})""".format(pairs_path)
        command += "'"

        subprocess.check_call(command, shell=True, stdout=outstream)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


@register_subcommand
@add_arg_help
@click.option(
    "--reference-file", "-r", help="Header file for transfer", type=str, required=True
)
def transfer(pairs_path, output, reference_file, **kwargs):
    """
    Transfer the header from one pairs file to another

    """
    transfer_py(pairs_path, output, reference_file, **kwargs)


def transfer_py(pairs_path, output, reference_file, **kwargs):

    instream = fileio.auto_open(
        pairs_path,
        mode="r",
        nproc=kwargs.get("nproc_in"),
        command=kwargs.get("cmd_in", None),
    )
    header, body_stream = headerops.get_header(instream, ignore_warning=True)

    # Read the header from reference file
    instream_header = fileio.auto_open(
        reference_file,
        mode="r",
        nproc=kwargs.get("nproc_in"),
        command=kwargs.get("cmd_in", None),
    )
    reference_header, _ = headerops.get_header(instream_header)
    # Close the reference stream after extraction of the header:
    if instream_header != sys.stdin:
        instream_header.close()

    reference_columns = headerops.extract_column_names(reference_header)

    # Check that the number of columns in the body corresponds to the header:
    if not headerops.validate_cols(instream, reference_columns):
        raise ValueError(
            f"Number of columns mismatch:\n\t#columns: {headerops.SEP_COLS.join(reference_columns)}\n\t{body_stream.readline()}"
        )

    ########
    # Write the output after successful checks:
    outstream = (
        fileio.auto_open(
            output,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output
        else sys.stdout
    )

    reference_header = headerops.append_new_pg(
        reference_header, ID=UTIL_NAME, PN=UTIL_NAME
    )
    outstream.writelines((l + "\n" for l in reference_header))
    outstream.flush()

    if body_stream == sys.stdin:
        for line in body_stream:
            outstream.write(line)
    else:
        command = r"""
                   /bin/bash -c 'export LC_COLLATE=C; export LANG=C; cat """

        if kwargs.get("cmd_in", None):
            command += r""" <(cat {} | {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path, kwargs["cmd_in"]
            )
        elif pairs_path.endswith(".gz"):
            command += (
                r""" <(bgzip -dc -@ {} {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                    kwargs["nproc_in"], pairs_path
                )
            )
        elif pairs_path.endswith(".lz4"):
            command += r""" <(lz4c -dc {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path
            )
        else:
            command += r""" <(sed -n -e '\''/^[^#]/,$p'\'' {})""".format(pairs_path)
        command += "'"

        subprocess.check_call(command, shell=True, stdout=outstream)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


@register_subcommand
@add_arg_help
@click.option(
    "--columns",
    "-c",
    help=f"Comma-separated list of columns to be added, e.g.: {','.join(pairsam_format.COLUMNS)}",
    type=str,
    required=True,
)
def set_columns(pairs_path, output, columns, **kwargs):
    """
    Add the columns to the .pairs/pairsam file
    """
    set_columns_py(pairs_path, output, columns, **kwargs)


def set_columns_py(pairs_path, output, columns, **kwargs):
    instream = (
        fileio.auto_open(
            pairs_path,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        if pairs_path
        else sys.stdin
    )
    outstream = (
        fileio.auto_open(
            output,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output
        else sys.stdout
    )

    header, body_stream = headerops.get_header(instream)
    header = headerops.set_columns(header, columns.split(","))
    outstream.writelines((l + "\n" for l in header))
    outstream.flush()

    if body_stream == sys.stdin:
        for line in body_stream:
            outstream.write(line)
    else:
        command = r"""
                   /bin/bash -c 'export LC_COLLATE=C; export LANG=C; cat """

        if kwargs.get("cmd_in", None):
            command += r""" <(cat {} | {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path, kwargs["cmd_in"]
            )
        elif pairs_path.endswith(".gz"):
            command += (
                r""" <(bgzip -dc -@ {} {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                    kwargs["nproc_in"], pairs_path
                )
            )
        elif pairs_path.endswith(".lz4"):
            command += r""" <(lz4c -dc {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path
            )
        else:
            command += r""" <(sed -n -e '\''/^[^#]/,$p'\'' {})""".format(pairs_path)
        command += "'"

        subprocess.check_call(command, shell=True, stdout=outstream)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


@register_subcommand
@add_arg_help
@click.option(
    "--reference-file",
    "-r",
    help="Header file for comparison (optional)",
    type=str,
    required=False,
    default="",
)
@click.option(
    "--reference-columns",
    "-c",
    help=f"Comma-separated list of columns fro check (optional), e.g.: {','.join(pairsam_format.COLUMNS)}",
    type=str,
    required=False,
    default="",
)
def validate_columns(pairs_path, output, reference_file, reference_columns, **kwargs):
    """
    Validate the columns of the .pairs/pairsam file against reference or within file.
    If the checks pass, then returns full pairs file. Otherwise throws an exception.

    If reference_file is provided, check:
        1) columns are the same between pairs and reference_file
        2) number of columns in the pairs body is the same as the number of columns

    If reference_columns are provided, check:
        1) pairs columns are the same as provided
        2) number of columns in the pairs body is the same as the number of columns

    If no reference_file or columns, then check only the number of columns in the pairs body.
    Checks only the first line in the pairs stream!

    """
    validate_columns_py(pairs_path, output, reference_file, reference_columns, **kwargs)


def validate_columns_py(
    pairs_path, output, reference_file, reference_columns, **kwargs
):

    instream = fileio.auto_open(
        pairs_path,
        mode="r",
        nproc=kwargs.get("nproc_in"),
        command=kwargs.get("cmd_in", None),
    )
    header, body_stream = headerops.get_header(instream)
    pairs_columns = headerops.extract_column_names(header)

    # Convert reference columns string into list, if provided
    if reference_columns:
        reference_columns = reference_columns.split(",")

    # Read the header from reference file
    if reference_file:
        instream_header = fileio.auto_open(
            reference_file,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        reference_header, _ = headerops.get_header(instream_header)
        # Close the reference stream after extraction of the header:
        if instream_header != sys.stdin:
            instream_header.close()

        if reference_columns:
            warnings.warn(
                "--reference-columns are ignored, as --reference-file is provided"
            )

        reference_columns = headerops.extract_column_names(reference_header)

    if reference_columns:
        if pairs_columns != reference_columns:
            raise ValueError(
                f"Pairs columns differ from reference columns:\n\t{pairs_columns}\n\t{reference_columns}"
            )

    # Check that the number of columns in the body corresponds to the header:
    if not headerops.validate_cols(instream, pairs_columns):
        raise ValueError(
            f"Number of columns mismatch:\n\t#columns: {headerops.SEP_COLS.join(pairs_columns)}\n\t{body_stream.readline()}"
        )

    ########
    # Write the output after successful checks:
    outstream = (
        fileio.auto_open(
            output,
            mode="w",
            nproc=kwargs.get("nproc_out"),
            command=kwargs.get("cmd_out", None),
        )
        if output
        else sys.stdout
    )
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l + "\n" for l in header))
    outstream.flush()

    if body_stream == sys.stdin:
        for line in body_stream:
            outstream.write(line)
    else:
        command = r"""
                   /bin/bash -c 'export LC_COLLATE=C; export LANG=C; cat """

        if kwargs.get("cmd_in", None):
            command += r""" <(cat {} | {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path, kwargs["cmd_in"]
            )
        elif pairs_path.endswith(".gz"):
            command += (
                r""" <(bgzip -dc -@ {} {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                    kwargs["nproc_in"], pairs_path
                )
            )
        elif pairs_path.endswith(".lz4"):
            command += r""" <(lz4c -dc {} | sed -n -e '\''/^[^#]/,$p'\'')""".format(
                pairs_path
            )
        else:
            command += r""" <(sed -n -e '\''/^[^#]/,$p'\'' {})""".format(pairs_path)
        command += "'"

        subprocess.check_call(command, shell=True, stdout=outstream)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()


if __name__ == "__main__":
    header()
