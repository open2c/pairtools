#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import sys
import click
import subprocess
import shutil
import warnings
import logging

from ..lib import fileio, pairsam_format, headerops
from . import cli, common_io_options

UTIL_NAME = "pairtools_sort"

logging.basicConfig(level=logging.DEBUG)

@cli.command()
@click.argument("pairs_path", type=str, required=False)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="output pairs file. If not specified, output to stdout.",
)
@click.option(
    "--c1", type=str, default="1", help="Chrom 1 column (1-based index or name, e.g., 'chrom1')",
)
@click.option(
    "--c2", type=str, default="3", help="Chrom 2 column (1-based index or name, e.g., 'chrom2')",
)
@click.option(
    "--p1", type=str, default="2", help="Position 1 column (1-based index or name, e.g., 'pos1')",
)
@click.option(
    "--p2", type=str, default="4", help="Position 2 column (1-based index or name, e.g., 'pos2')",
)
@click.option(
    "--pt", type=str, default="7", help="Pair type column (1-based index or name, e.g., 'pair_type'), optional",
    required=False,
)
@click.option(
    "--extra-col", nargs=1, type=str, multiple=True,
    help="Extra column (name or index) for sorting. Can be repeated.",
)
@click.option(
    "--nproc", type=int, default=1, show_default=True, help="Number of processes for sorting.",
)
@click.option(
    "--tmpdir", type=str, default="", help="Temporary folder for sorting intermediates.",
)
@click.option(
    "--memory", type=str, default="2G", show_default=True, help="Memory for sorting.",
)
@click.option(
    "--compress-program", type=str, default="auto", show_default=True,
    help="Compression binary (e.g., gzip, lz4c). 'auto' uses lz4c if available, else gzip.",
)
@common_io_options
def sort(pairs_path, output, c1, c2, p1, p2, pt, extra_col, nproc, tmpdir, memory, compress_program, **kwargs):
    sort_py(pairs_path, output, c1, c2, p1, p2, pt, extra_col, nproc, tmpdir, memory, compress_program, **kwargs)

def sort_py(pairs_path, output, c1, c2, p1, p2, pt, extra_col, nproc, tmpdir, memory, compress_program, **kwargs):
    logging.debug("Starting sort_py")
    instream = fileio.auto_open(pairs_path, mode="r", nproc=kwargs.get("nproc_in"), command=kwargs.get("cmd_in"))
    outstream = fileio.auto_open(output, mode="w", nproc=kwargs.get("nproc_out"), command=kwargs.get("cmd_out"))

    header, body_stream = headerops.get_header(instream)
    logging.debug(f"Header: {header}")
    header = headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    header = headerops.mark_header_as_sorted(header)

    column_names = headerops.extract_column_names(header)
    logging.debug(f"Column names: {column_names}")
    col_indices = headerops.map_columns(
        column_names, {"c1": c1, "c2": c2, "p1": p1, "p2": p2, "pt": pt or ""}
    )
    logging.debug(f"Column indices: {col_indices}")

    canonical_order = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
    if "pair_type" in col_indices and col_indices["pair_type"] is not None:
        canonical_order.append("pair_type")
    extra_cols = [headerops.parse_column(col, column_names) for col in extra_col]
    canonical_order.extend([column_names[i] for i in extra_cols if column_names[i] not in canonical_order])
    logging.debug(f"Canonical order: {canonical_order}")

    header = headerops.set_columns(header, canonical_order)
    outstream.writelines((l + "\n" for l in header))
    outstream.flush()
    logging.debug("Header written")

    body_lines = list(body_stream)
    logging.debug(f"Body lines count: {len(body_lines)}")
    if not body_lines:
        logging.warning("No body lines to sort")
        if instream != sys.stdin:
            instream.close()
        if outstream != sys.stdout:
            outstream.close()
        return

    compress_program = "lz4c" if (compress_program == "auto" and shutil.which("lz4c")) else "gzip"
    logging.debug(f"Compress program: {compress_program}")

    sort_keys = [
        (col_indices["chrom1"] + 1, ""),
        (col_indices["pos1"] + 1, "n"),
        (col_indices["chrom2"] + 1, ""),
        (col_indices["pos2"] + 1, "n"),
    ]
    if "pair_type" in col_indices and col_indices["pair_type"] is not None:
        sort_keys.append((col_indices["pair_type"] + 1, ""))
    for col in extra_cols:
        colname = column_names[col]
        sort_keys.append((col + 1, "n" if issubclass(pairsam_format.DTYPES_PAIRSAM.get(colname, str), int) else ""))

    cols = " ".join(f"-k {i},{i}{mode}" for i, mode in sort_keys)
    command = (
        f"sort {cols} --stable --field-separator='\t' "  # Correct tab separator
        f"--parallel={nproc} {f'--temporary-directory={tmpdir}' if tmpdir else ''} "
        f"-S {memory} --compress-program={compress_program}"
    ).strip()
    logging.debug(f"Sort command: {command}")

    try:
        process = subprocess.run(
            command,
            input="\n".join(body_lines),
            text=True,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10
        )
        logging.debug(f"Subprocess return code: {process.returncode}")
        if process.returncode != 0:
            logging.error(f"Sort failed: {process.stderr}")
            raise subprocess.CalledProcessError(process.returncode, command, process.stdout, process.stderr)

        sorted_lines = [line for line in process.stdout.splitlines() if line.strip()]  # Filter empty lines
        logging.debug(f"Sorted lines count: {len(sorted_lines)}")
        for line in sorted_lines:
            fields = line.split(pairsam_format.PAIRSAM_SEP)
            if len(fields) < len(canonical_order):
                logging.error(f"Line too short: {line}, expected {len(canonical_order)} fields, got {len(fields)}")
                raise ValueError(f"Invalid line: {line}")
            reordered = [fields[col_indices.get(col, column_names.index(col))] for col in canonical_order]
            outstream.write(pairsam_format.PAIRSAM_SEP.join(reordered) + "\n")
        outstream.flush()
        logging.debug("Sorted output written")
    except subprocess.TimeoutExpired as e:
        logging.error(f"Subprocess timed out: {e.stderr}")
        raise

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
    logging.debug("sort_py completed")

if __name__ == "__main__":
    sort()