from collections import defaultdict
import sys
import copy
import itertools
import warnings

import numpy as np
import pandas as pd

from .. import __version__
from . import pairsam_format
from .fileio import ParseError, get_stream_handlers

from .._logging import get_logger

logger = get_logger()

PAIRS_FORMAT_VERSION = "1.0.0"
SEP_COLS = " "
SEP_CHROMS = " "
COMMENT_CHAR = "#"



def get_header(instream, comment_char=COMMENT_CHAR, ignore_warning=False):
    """Returns a header from the stream and an the reaminder of the stream
    with the actual data.
    Parameters
    ----------
    instream : a file object
        An input stream.
    comment_char : str
        The character prepended to header lines (use '@' when parsing sams,
        '#' when parsing pairsams).
    ignore_warning : bool
        If True, then no warning will be generated if header of pairs file is empty.
    Returns
    -------
    header : list
        The header lines, stripped of terminal spaces and newline characters.
    remainder_stream : stream/file-like object
        Stream with the remaining lines.

    """
    header = []
    if not comment_char:
        raise ValueError("Please, provide a comment char!")
    comment_byte = comment_char.encode()

    readline_f, peek_f = get_stream_handlers(instream)
    current_peek = peek_f(1)
    while current_peek.startswith(comment_byte):
        # consuming a line from buffer guarantees
        # that the remainder of the buffer starts
        # with the beginning of the line.
        line = readline_f()
        if isinstance(line, bytes):
            line = line.decode()
        # append line to header, since it does start with header
        header.append(line.rstrip('\n'))
        # peek into the remainder of the instream
        current_peek = peek_f(1)
    # apparently, next line does not start with the comment
    # return header and the instream, advanced to the beginning of the data

    if len(header) == 0 and not ignore_warning:
        logger.warning(
            "Headerless input, please, add the header by `pairtools header generate` or `pairtools header transfer`"
        )

    return header, instream


def extract_fields(header, field_name, save_rest=False):
    """
    Extract the specified fields from the pairs header and return
    a list of corresponding values, even if a single field was found.
    Additionally, can return the list of intact non-matching entries.
    """

    fields = []
    rest = []
    for l in header:
        if l.lstrip(COMMENT_CHAR).startswith(field_name + ":"):
            fields.append(l.split(":", 1)[1].rstrip('\n').lstrip())
        elif save_rest:
            rest.append(l)

    if save_rest:
        return fields, rest
    else:
        return fields


def extract_column_names(header):
    """
    Extract column names from header lines.
    """
    columns = extract_fields(header, "columns")

    if len(columns) != 0:
        return columns[0].split(SEP_COLS)
    else:
        return []


def validate_cols(stream, columns):
    """
    Validate that the number of columns coincides between stream and columns.
    Checks only the first line in the pairs stream!

    Note that it irreversibly removes the header from the stream.

    Parameters
    ----------
    stream: input stream, body or full .pairs file
    columns: columns to validate against

    Returns
    -------
    True if the number of columns is identical between file and columns
    """

    comment_byte = COMMENT_CHAR.encode()
    readline_f, peek_f = get_stream_handlers(stream)

    current_peek = peek_f(1)
    while current_peek.startswith(comment_byte):
        # consuming a line from buffer guarantees
        # that the remainder of the buffer starts
        # with the beginning of the line.
        line = readline_f()
        # peek into the remainder of the instream
        current_peek = peek_f(1)

    line = readline_f()
    if isinstance(line, bytes):
        line = line.decode()

    ncols_body = len(line.split(pairsam_format.PAIRSAM_SEP))
    ncols_reference = (
        len(columns) if isinstance(columns, list) else columns.split(SEP_COLS)
    )

    return ncols_body == ncols_reference


def validate_header_cols(stream, header):
    """Validate that the number of columns corresponds between the stream and header"""

    columns = extract_column_names(header)
    return validate_cols(stream, header)


def is_empty_header(header):
    if len(header) == 0:
        return True
    if not header[0].startswith("##"):
        return True
    else:
        return False


def extract_chromsizes(header):
    """
    Extract chromosome sizes from header lines.
    """

    chromsizes_str = extract_fields(header, "chromsize")
    chromsizes_str = list(zip(*[s.split(SEP_CHROMS) for s in chromsizes_str]))
    chromsizes = pd.Series(data=chromsizes_str[1], index=chromsizes_str[0]).astype(
        np.int64
    )

    return chromsizes


def get_chromsizes_from_pysam_header(samheader):
    """Convert pysam header to pairtools chromosomes dict (ordered by Python default since 3.7).

    Example of pysam header converted to dict:
    dict([
        ('SQ', [{'SN': 'chr1', 'LN': 248956422},
         {'SN': 'chr10', 'LN': 133797422},
         {'SN': 'chr11', 'LN': 135086622},
         {'SN': 'chr12', 'LN': 133275309}]),
        ('PG', [{'ID': 'bwa', 'PN': 'bwa', 'VN': '0.7.17-r1188', 'CL': 'bwa mem -t 8 -SP -v1 hg38.fa test_1.1.fastq.gz test_2.1.fastq.gz'}])
    ])
    """
    SQs = samheader.to_dict()["SQ"]
    chromsizes = [(sq["SN"], int(sq["LN"])) for sq in SQs]
    return dict(chromsizes)


def get_chromsizes_from_file(chroms_file):
    """
    Produce an "enumeration" of chromosomes based on the list
    of chromosomes
    """
    chrom_sizes = dict()
    with open(chroms_file, "rt") as f:
        for line in f:
            chrom, size = line.strip().split("\t")
            chrom_sizes[chrom] = int(size)

    return chrom_sizes


def get_chromsizes_from_pysam_header(samheader):
    """Convert pysam header to pairtools chromosomes (ordered dict).

    Example of pysam header converted to dict:
    dict([
        ('SQ', [{'SN': 'chr1', 'LN': 248956422},
         {'SN': 'chr10', 'LN': 133797422},
         {'SN': 'chr11', 'LN': 135086622},
         {'SN': 'chr12', 'LN': 133275309}]),
        ('PG', [{'ID': 'bwa', 'PN': 'bwa', 'VN': '0.7.17-r1188', 'CL': 'bwa mem -t 8 -SP -v1 hg38.fa test_1.1.fastq.gz test_2.1.fastq.gz'}])
    ])
    """
    SQs = samheader.to_dict()["SQ"]
    chromsizes = [(sq["SN"], int(sq["LN"])) for sq in SQs]
    return dict(chromsizes)


def get_chrom_order(chroms_file, sam_chroms=None):
    """
    Produce an "enumeration" of chromosomes based on the list
    of chromosomes
    """
    chrom_enum = dict()
    i = 1
    with open(chroms_file, "rt") as f:
        for line in f:
            chrom = line.strip().split("\t")[0]
            if chrom and ((not sam_chroms) or (chrom in sam_chroms)):
                chrom_enum[chrom] = i
                i += 1

    if sam_chroms:
        remaining = sorted(
            chrom for chrom in sam_chroms if chrom not in chrom_enum.keys()
        )

        for chrom in remaining:
            chrom_enum[chrom] = i
            i += 1

    return chrom_enum


def make_standard_pairsheader(
    assembly=None,
    chromsizes=None,
    columns=pairsam_format.COLUMNS,
    shape="upper triangle",
):
    header = []
    header.append("## pairs format v{}".format(PAIRS_FORMAT_VERSION))
    header.append("#shape: {}".format(shape))

    header.append(
        "#genome_assembly: {}".format(assembly if assembly is not None else "unknown")
    )

    if chromsizes is not None:
        try:
            chromsizes = chromsizes.items()
        except AttributeError:
            pass
        for chrom, length in chromsizes:
            header.append("#chromsize: {} {}".format(chrom, length))

    header.append("#columns: " + SEP_COLS.join(columns))

    return header


def subset_chroms_in_pairsheader(header, chrom_subset):
    new_header = []
    for line in header:
        if line.startswith("#chromsize:"):
            if line.strip().split()[1] in chrom_subset:
                new_header.append(line)
        elif line.startswith("#chromosomes:"):
            line = SEP_CHROMS.join(
                ["#chromosomes:"]
                + [c for c in line.strip().split()[1:] if c in chrom_subset]
            )
            new_header.append(line)
        else:
            new_header.append(line)
    return new_header


def insert_samheader(header, samheader):
    """Insert samheader into header."""
    new_header = [l for l in header if not l.startswith("#columns")]
    if samheader:
        new_header += ["#samheader: " + l for l in samheader]
    new_header += [l for l in header if l.startswith("#columns")]
    return new_header


def insert_samheader_pysam(header, samheader):
    """Insert samheader into header,pysam version."""
    new_header = [l for l in header if not l.startswith("#columns")]
    if samheader:
        new_header += ["#samheader: " + l for l in str(samheader).strip().split("\n")]
    new_header += [l for l in header if l.startswith("#columns")]
    return new_header


def mark_header_as_sorted(header):
    header = copy.deepcopy(header)
    if is_empty_header(header):
        raise Exception("Input file is not valid .pairs, has no header or is empty.")
    if not any([l.startswith("#sorted") for l in header]):
        if header[0].startswith("##"):
            header.insert(1, "#sorted: chr1-chr2-pos1-pos2")
        else:
            header.insert(0, "#sorted: chr1-chr2-pos1-pos2")
    for i in range(len(header)):
        if header[i].startswith("#chromosomes"):
            chroms = header[i][12:].strip().split(SEP_CHROMS)
            header[i] = "#chromosomes: {}".format(SEP_CHROMS.join(sorted(chroms)))
    return header


def append_new_pg(header, ID="", PN="", VN=None, CL=None, force=False):
    header = copy.deepcopy(header)
    if is_empty_header(header):
        raise Exception("Input file is not valid .pairs, has no header or is empty.")
    samheader, other_header = extract_fields(header, "samheader", save_rest=True)
    new_samheader = _add_pg_to_samheader(samheader, ID, PN, VN, CL, force)
    new_header = insert_samheader(other_header, new_samheader)
    return new_header


def _update_header_entry(header, field, new_value):
    header = copy.deepcopy(header)
    found = False
    newline = "#{}: {}".format(field, new_value)
    for i in range(len(header)):
        if header[i].startswith(COMMENT_CHAR + field):
            header[i] = newline
            found = True
    if not found:
        if header[-1].startswith("#columns"):
            header.insert(-1, newline)
        else:
            header.append(newline)
    return header


def _add_pg_to_samheader(samheader, ID="", PN="", VN=None, CL=None, force=False):
    """Append a @PG record to an existing sam header. If the header comes
    from a merged file and thus has multiple chains of @PG, append the
    provided PG to all of the chains, adding the numerical suffix of the
    branch to the ID.

    Parameters
    ----------

    header : list of str
    ID, PN, VN, CL : std
        The keys of a new @PG record. If absent, VN is the version of
        pairtools and CL is taken from sys.argv.
    force : bool
        If True, ignore the inconsistencies among @PG records of the existing
        header.

    Returns
    -------
    new_header : list of str
        A list of new headers lines, stripped of newline characters.
    """

    if VN is None:
        VN = __version__
    if CL is None:
        CL = " ".join(sys.argv)

    pre_pg_header = [
        line.strip()
        for line in samheader
        if line.startswith("@HD") or line.startswith("@SQ") or line.startswith("@RG")
    ]

    post_pg_header = [
        line.strip()
        for line in samheader
        if not line.startswith("@HD")
        and (not line.startswith("@SQ"))
        and (not line.startswith("@RG"))
        and (not line.startswith("@PG"))
    ]

    pg_chains = _parse_pg_chains(samheader, force=force)

    for i, br in enumerate(pg_chains):
        new_pg = {"ID": ID, "PN": PN, "VN": VN, "CL": CL}
        new_pg["PP"] = br[-1]["ID"]
        if len(pg_chains) > 1:
            new_pg["ID"] = new_pg["ID"] + "-" + str(i + 1) + "." + str(len(br) + 1)
        new_pg["raw"] = _format_pg(**new_pg)
        br.append(new_pg)

    new_header = (
        pre_pg_header + [pg["raw"] for br in pg_chains for pg in br] + post_pg_header
    )

    return new_header


def _format_pg(**kwargs):
    out = ["@PG"] + [
        "{}:{}".format(field, kwargs[field])
        for field in ["ID", "PN", "CL", "PP", "DS", "VN"]
        if field in kwargs
    ]
    return "\t".join(out)


def _parse_pg_chains(header, force=False):
    pg_chains = []

    parsed_pgs = []
    for l in header:
        if l.startswith("@PG"):
            tag_value_pairs = l.strip().split("\t")[1:]
            if not all(":" in tvp for tvp in tag_value_pairs):
                warnings.warn(
                    f"Skipping the following @PG line, as it does not follow the SAM header standard of TAG:VALUE: {l}"
                )
                continue
            parsed_tvp = dict(
                [tvp.split(":", maxsplit=1) for tvp in tag_value_pairs if ":" in tvp]
            )
            if parsed_tvp:
                parsed_tvp["raw"] = l.strip()
                parsed_pgs.append(parsed_tvp)

    while True:
        if len(parsed_pgs) == 0:
            break

        for i in range(len(parsed_pgs)):
            pg = parsed_pgs[i]
            if "PP" not in pg:
                pg_chains.append([pg])
                parsed_pgs.pop(i)
                break
            else:
                matching_chains = [
                    branch for branch in pg_chains if branch[-1]["ID"] == pg["PP"]
                ]
                if len(matching_chains) > 1:
                    if force:
                        matching_chains[0].append(pg)
                        parsed_pgs.pop(i)
                        break

                    else:
                        raise ParseError(
                            "Multiple @PG records with the IDs identical to the PP field of another record:\n"
                            + "\n".join([br[-1]["raw"] for br in matching_chains])
                            + "\nvs\n"
                            + pg["raw"]
                        )

                if len(matching_chains) == 1:
                    matching_chains[0].append(pg)
                    parsed_pgs.pop(i)
                    break

            if force:
                pg_chains.append([pg])
                parsed_pgs.pop(i)
                break
            else:
                raise ParseError(
                    "Cannot find the parental @PG record for the @PG records:\n"
                    + "\n".join([pg["raw"] for pg in parsed_pgs])
                )

    return pg_chains


def _toposort(dag, tie_breaker):
    """
    Topological sort on a directed acyclic graph

    Uses Kahn's algorithm with a custom tie-breaking option. The
    dictionary ``dag`` can be interpreted in two ways:

    1. A dependency graph (i.e. arcs point from values to keys),
       and the generator yields items with no dependences followed
       by items that depend on previous ones.

    2. Arcs point from keys to values, in which case the generator produces
       a **reverse** topological ordering of the nodes.

    Parameters
    ----------
    dag: dict of nodes to sets of nodes
        Directed acyclic graph encoded as a dictionary.
    tie_breaker: callable
        Function that picks a tie breaker from a set of nodes with no
        unprocessed dependences.

    Returns
    -------
    Generator

    Notes
    -----
    See <https://en.wikipedia.org/wiki/Topological_sorting for more info>.
    Based in part on activestate recipe:
    <http://code.activestate.com/recipes/578272-topological-sort/> by Sam
    Denton (MIT licensed).

    """
    # Drop self-edges.
    for k, v in dag.items():
        v.discard(k)

    # Find all nodes that don't depend on anything
    # and include them with empty dependencies.
    indep_nodes = set.union(*dag.values()) - set(dag.keys())
    dag.update({node: set() for node in indep_nodes})

    while True:
        if not indep_nodes:
            break
        out = tie_breaker(indep_nodes)
        indep_nodes.discard(out)
        del dag[out]

        yield out

        for node, deps in dag.items():
            deps.discard(out)
            if len(deps) == 0:
                indep_nodes.add(node)

    if len(dag) != 0:
        raise ValueError("Circular dependencies exist: {} ".format(list(dag.items())))


def merge_chrom_lists(*lsts):
    sentinel = "!NONE!"

    g = defaultdict(set)
    for lst in lsts:
        if len(lst) == 1:
            g[lst[0]].add(sentinel)
        for a, b in zip(lst[:-1], lst[1:]):
            g[b].add(a)
    if len(g) == 0:
        return []

    chrom_list = list(_toposort(g.copy(), tie_breaker=min))
    if sentinel in chrom_list:
        chrom_list.remove(sentinel)
    chrom_list = sorted(chrom_list)
    return chrom_list


def _merge_samheaders(samheaders, force=False):
    # first, append an HD line if it is present in any files
    # if different lines are present, raise an error
    HDs = set.union(
        *[
            set(line for line in samheader if line.startswith("@HD"))
            for samheader in samheaders
        ]
    )
    if len(HDs) > 1 and not force:
        raise ParseError("More than one unique @HD line is found in samheaders!")
    HDs = [list(HDs)[0]] if HDs else []

    # second, confirm that all files had the same SQ lines
    # add SQs from the first file, keeping its order
    SQs = [
        set(line for line in samheader if line.startswith("@SQ"))
        for samheader in samheaders
    ]
    common_SQs = set.intersection(*SQs)

    SQs_same = all([len(samheader) == len(common_SQs) for samheader in SQs])

    if not SQs_same and not (force):
        raise ParseError("The SQ (sequence) lines of the sam headers are not identical")
    SQs = [line for line in samheaders[0] if line.startswith("@SQ")]

    # third, append _all_ PG chains, adding a unique index according to the
    # provided merging order
    PGs = []
    for i, samheader in enumerate(samheaders):
        for line in samheader:
            if line.startswith("@PG"):
                split_line = line.split("\t")
                for j in range(len(split_line)):
                    if split_line[j].startswith("ID:") or split_line[j].startswith(
                        "PP:"
                    ):
                        split_line[j] = split_line[j] + "-" + str(i + 1)
                PGs.append("\t".join(split_line))

    # finally, add all residual unique lines
    rest = sum(
        [
            list(
                set(
                    line
                    for line in samheader
                    if (not line.startswith("@HD"))
                    and (not line.startswith("@SQ"))
                    and (not line.startswith("@PG"))
                )
            )
            for samheader in samheaders
        ],
        [],
    )

    new_header = []
    new_header += HDs
    new_header += SQs
    new_header += PGs
    new_header += rest

    return new_header


def _merge_pairheaders(pairheaders, force=False):
    new_header = []

    # first, add all keys that are expected to be the same among all headers
    keys_expected_identical = [
        "## pairs format",
        "#sorted:",
        "#shape:",
        "#genome_assembly:",
        "#columns:",
    ]

    keys_orginal = [l.split()[0] for header in pairheaders for l in header]

    for k in keys_expected_identical:
        lines = [[l for l in header if l.startswith(k)] for header in pairheaders]
        same = all([l == lines[0] for l in lines])
        if not (same or force):
            raise ParseError(
                "The following header entries must be the same "
                "the merged files: {}".format(k)
            )
        new_header += lines[0]

    # second, merge and add the chromsizes fields.
    chrom_lists = []
    chromsizes = {}
    for header in pairheaders:
        chromlist = []
        for line in header:
            if line.startswith("#chromsize:"):
                chrom, length = line.strip("#chromsize:").split()
                chromsizes[chrom] = length
                chromlist.append(chrom)
        chrom_lists.append(chromlist)

    chroms_merged = merge_chrom_lists(*chrom_lists)
    if "#chromosomes:" in keys_orginal:
        chrom_line = "#chromosomes: {}".format(" ".join(chroms_merged))
        new_header.extend([chrom_line])

    chromsize_lines = [
        "#chromsize: {} {}".format(chrom, chromsizes[chrom]) for chrom in chroms_merged
    ]
    new_header.extend(chromsize_lines)

    # finally, add a sorted list of other unique fields
    other_lines = sorted(
        set(
            l
            for h in pairheaders
            for l in h
            if not any(
                l.startswith(k)
                for k in keys_expected_identical + ["#chromosomes", "#chromsize"]
            )
        )
    )
    if other_lines:
        if new_header[-1].startswith("#columns"):
            new_header = new_header[:-1] + other_lines + [new_header[-1]]
        else:
            new_header = new_header + other_lines

    return new_header


def all_same_columns(pairheaders):
    key_target = "#columns:"
    lines = [[l for l in header if l.startswith(key_target)] for header in pairheaders]
    all_same = all([l == lines[0] for l in lines])
    return all_same


def merge_headers(headers, force=False):
    samheaders, pairheaders = zip(
        *[extract_fields(h, "samheader", save_rest=True) for h in headers]
    )
    # HD headers contain information that becomes invalid after processing
    # with distiller. Do not print into the output.
    new_pairheader = _merge_pairheaders(pairheaders, force=False)
    new_samheader = _merge_samheaders(samheaders, force=force)

    new_header = insert_samheader(new_pairheader, new_samheader)

    return new_header


def append_columns(header, columns):
    """
    Appends columns to the header, separated by SEP_COLS

    Parameters
    ----------
    header: Previous header
    columns: List of column names to append

    Returns
    -------
    Modified header (appended columns to the field "#columns")
    """
    for i in range(len(header)):
        if header[i].startswith("#columns: "):
            header[i] += SEP_COLS + SEP_COLS.join(columns)
    return header


def get_colnames(header):
    """
    Get column names of the header, separated by SEP_COLS

    Parameters
    ----------
    header: Previous header

    Returns
    -------
    List of column names
    """
    for i in range(len(header)):
        if header[i].startswith("#columns: "):
            columns = header[i].split(SEP_COLS)[1:]
            return columns
    return []


def set_columns(header, columns):
    """
    Set columns to the header, separated by SEP_COLS

    Parameters
    ----------
    header: Previous header
    columns: List of column names to append

    Returns
    -------
    Modified header (appended columns to the field "#columns")
    """
    for i in range(len(header)):
        if header[i].startswith("#columns:"):
            header[i] = "#columns:" + SEP_COLS + SEP_COLS.join(columns)
    return header


# def _guess_genome_assembly(samheader):
#    PG = [l for l in samheader if l.startswith('@PG') and '\tID:bwa' in l][0]
#    CL = [field for field in PG.split('\t') if field.startswith('CL:')]
#
#    return ga
