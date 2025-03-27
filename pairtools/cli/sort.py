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
@click.option("-o", "--output", type=str, default="", help="output pairs file. If the path ends with .gz or .lz4, the output is compressed by bgzip or lz4, correspondingly. By default, the output is printed into stdout.")
@click.option("--c1", type=str, default=pairsam_format.COLUMNS_PAIRS[1], help=f"Chrom 1 column; default {pairsam_format.COLUMNS_PAIRS[1]} (index or name, e.g., '1' or 'chrom1') [input format option]")
@click.option("--c2", type=str, default=pairsam_format.COLUMNS_PAIRS[3], help=f"Chrom 2 column; default {pairsam_format.COLUMNS_PAIRS[3]} (index or name, e.g., '3' or 'chrom2') [input format option]")
@click.option("--p1", type=str, default=pairsam_format.COLUMNS_PAIRS[2], help=f"Position 1 column; default {pairsam_format.COLUMNS_PAIRS[2]} (index or name, e.g., '2' or 'pos1') [input format option]")
@click.option("--p2", type=str, default=pairsam_format.COLUMNS_PAIRS[4], help=f"Position 2 column; default {pairsam_format.COLUMNS_PAIRS[4]} (index or name, e.g., '4' or 'pos2') [input format option]")
@click.option("--pt", type=str, default=pairsam_format.COLUMNS_PAIRS[7], help=f"Pair type column; default {pairsam_format.COLUMNS_PAIRS[7]} (index or name, e.g., '7' or 'pair_type'), optional [input format option]", required=False)
@click.option("--extra-col", nargs=1, type=str, multiple=True, help="Extra column (index or name) to sort by. Can be repeated. Example: --extra-col '5' --extra-col 'phase1'. [output format option]")
@click.option("--nproc", type=int, default=8, show_default=True, help="Number of processes to split the sorting work between (subprocess mode only).")
@click.option("--tmpdir", type=str, default="", help="Custom temporary folder for sorting intermediates (subprocess mode only).")
@click.option("--memory", type=str, default="2G", show_default=True, help="The amount of memory used by default (subprocess mode only).")
@click.option("--compress-program", type=str, default="auto", show_default=True, help="A binary to compress temporary sorted chunks in subprocess mode. Must decompress input with -d flag. Suggested: gzip, lzop, lz4c, snzip. 'auto' uses lz4c if available, else gzip.")
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

    col_map = {"c1": c1, "c2": c2, "p1": p1, "p2": p2, "pt": pt or ""}
    aliases = {"chr1": "chrom1", "chr2": "chrom2"}
    col_indices = {}
    for key, value in col_map.items():
        if value:
            if value.isnumeric():
                col_indices[key] = int(value) - 1
            else:
                try:
                    col_indices[key] = column_names.index(value)
                except ValueError:
                    alias = aliases.get(value)
                    if alias and alias in column_names:
                        col_indices[key] = column_names.index(alias)
                    else:
                        raise ValueError(f"Column '{value}' not found in header columns: {column_names}")
        elif key != "pt":
            raise ValueError(f"Required column '{key}' not specified and not found in header")

    logging.debug(f"Column indices: {col_indices}")
    extra_cols = [headerops.parse_column(col, column_names) for col in extra_col]

    outstream.writelines((l + "\n" for l in header))
    outstream.flush()
    logging.debug("Header written")

    body_lines = [line for line in body_stream if line.strip()]
    logging.debug(f"Body lines count: {len(body_lines)}")
    if not body_lines:
        logging.warning("No body lines to sort")
        if instream != sys.stdin:
            instream.close()
        if outstream != sys.stdout:
            outstream.close()
        return

    def sort_key(line):
        fields = line.split(pairsam_format.PAIRSAM_SEP)
        key = [
            fields[col_indices["c1"]],  # chrom1
            int(fields[col_indices["p1"]]),  # pos1
            fields[col_indices["c2"]],  # chrom2
            int(fields[col_indices["p2"]]),  # pos2
        ]
        if "pt" in col_indices and col_indices["pt"] is not None:
            key.append(fields[col_indices["pt"]])
        for col in extra_cols:
            key.append(int(fields[col]) if issubclass(pairsam_format.DTYPES_PAIRSAM.get(column_names[col], str), int) else fields[col])
        return tuple(key)

    if len(body_lines) < 1000:
        logging.debug("Using in-memory sorting")
        sorted_lines = sorted(body_lines, key=sort_key)
        for line in sorted_lines:
            outstream.write(line + "\n")
    else:
        logging.debug("Using subprocess sorting")
        compress_program = "lz4c" if (compress_program == "auto" and shutil.which("lz4c")) else "gzip"
        logging.debug(f"Compress program: {compress_program}")
        columns = [c1, p1, c2, p2]  # Updated order
        if pt:
            columns.append(pt)
        columns.extend(extra_col)
        cols = []
        for col in columns:
            colindex = int(col) - 1 if col.isnumeric() else column_names.index(col)
            colindex += 1
            cols.append(f"-k {colindex},{colindex}{'n' if issubclass(pairsam_format.DTYPES_PAIRSAM.get(column_names[colindex-1], str), int) else ''}")
        cols = " ".join(cols)
        command = f"/bin/bash -c 'export LC_COLLATE=C; export LANG=C; sort {cols} --stable --field-separator=$'\\t' --parallel={nproc} {f'--temporary-directory={tmpdir}' if tmpdir else ''} -S {memory} --compress-program={compress_program}'".strip()
        logging.debug(f"Sort command: {command}")
        with subprocess.Popen(command, stdin=subprocess.PIPE, bufsize=-1, shell=True, stdout=outstream) as process:
            stdin_wrapper = io.TextIOWrapper(process.stdin, "utf-8")
            for line in body_lines:
                stdin_wrapper.write(line)
            stdin_wrapper.flush()
            process.communicate()

    outstream.flush()
    logging.debug("Sorted output written")

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
    logging.debug("sort_py completed")

if __name__ == "__main__":
    sort()