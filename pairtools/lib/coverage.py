#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coverage analysis API for genomic pairs data.

This module provides functions to calculate coverage from pairs files,
supporting various filtering and binning options.
"""

import multiprocessing
from functools import partial
from typing import Iterator, Optional, Tuple, Union

import bioframe
import numpy as np
import pandas as pd
from pairtools.lib import fileio, headerops

pd.set_option("mode.copy_on_write", True)


def read_pairs(
    pairs: Union[str, object], chunksize: int, threads: int = 1
) -> Tuple[Iterator[pd.DataFrame], Union[dict, pd.Series]]:
    """
    Read pairs file and extract chromosome sizes.

    Args:
        pairs: Path to pairs file or file-like object
        chunksize: Number of pairs to process at once
        threads: Number of threads for reading

    Returns:
        Tuple of (pairs DataFrame iterator, chromsizes dict)
    """
    pairs_stream = (
        fileio.auto_open(
            pairs,
            mode="r",
            nproc=threads,
        )
        if isinstance(pairs, str)
        else pairs
    )

    header, pairs_body = headerops.get_header(pairs_stream)

    cols = headerops.extract_column_names(header)
    chromsizes = headerops.extract_chromsizes(header)

    pairs_df = pd.read_csv(
        pairs_body,
        header=None,
        names=cols,
        chunksize=chunksize,
        sep="\t",
        dtype={"chrom1": str, "chrom2": str},
        low_memory=False,
    )  # type: ignore
    return pairs_df, chromsizes


def intervals_to_increments(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert genomic intervals to increment/decrement events.

    Args:
        df: DataFrame with columns 'chrom', 'start', 'end'

    Returns:
        DataFrame with columns 'chrom', 'pos', 'inc'
    """
    inc_df = pd.concat(
        [
            df[["chrom", "start"]].rename(columns={"start": "pos"}).eval("inc=1"),
            df[["chrom", "end"]].rename(columns={"end": "pos"}).eval("inc=-1"),
        ]
    )
    return inc_df


def aggregate_increments(inc_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate increment events by position.

    Args:
        inc_df: DataFrame with increment events

    Returns:
        Aggregated DataFrame sorted by chromosome and position
    """
    inc_df = (
        inc_df.groupby(["chrom", "pos"])
        .sum()
        .reset_index()
        .sort_values(["chrom", "pos"], ignore_index=True)
    )
    return inc_df


def coverage_chunk(df_chunk: pd.DataFrame, binsize: int) -> pd.DataFrame:
    """
    Calculate coverage for a chunk of intervals.

    Args:
        df_chunk: DataFrame with genomic intervals
        binsize: Size of bins for aggregation

    Returns:
        DataFrame with coverage counts per bin
    """
    count_df = aggregate_increments(intervals_to_increments(df_chunk))

    count_df.insert(2, "count", count_df["inc"].cumsum())
    count_df["bin"] = count_df["pos"] // binsize
    if binsize > 1:
        count_df = count_df.groupby(["chrom", "bin"], as_index=False).agg(
            {"count": "sum"}
        )
        count_df["bin"] = count_df["bin"].astype(int)
    return count_df


def produce_chunks(
    pairs_stream: Iterator[pd.DataFrame], side: str, end: str
) -> Iterator[pd.DataFrame]:
    """
    Generate chunks of genomic intervals from pairs stream.

    Args:
        pairs_stream: Iterator of pairs DataFrames
        side: Which side to use (0=both, 1=first, 2=second)
        end: Which end to use (0=default, 3=3', 5=5')

    Yields:
        DataFrames with genomic intervals
    """
    if end == "0":
        end = ""
    for pairs_df in pairs_stream:
        if pairs_df.shape[0] == 0:
            continue
        if side != "0":
            pairs_df["start"] = pairs_df[f"pos{end}{side}"] - 1
            pairs_df["end"] = pairs_df[f"pos{end}{side}"]
            pairs_df["chrom"] = pairs_df[f"chrom{side}"]
            pairs_df["strand"] = pairs_df[f"strand{side}"]
            yield pairs_df[["chrom", "start", "end", "strand"]]
        else:
            for side in ["1", "2"]:
                pairs_df["start"] = pairs_df[f"pos{end}{side}"] - 1
                pairs_df["end"] = pairs_df[f"pos{end}{side}"]
                pairs_df["chrom"] = pairs_df[f"chrom{side}"]
                pairs_df["strand"] = pairs_df[f"strand{side}"]
                yield pairs_df[["chrom", "start", "end", "strand"]]


def process_chunk(
    pairs_chunk: pd.DataFrame,
    shift: int,
    binsize: int,
    chromsizes: Union[dict, pd.Series],
) -> pd.DataFrame:
    """
    Process a single chunk of pairs data.

    Args:
        pairs_chunk: DataFrame with pairs data
        shift: Shift value for strand-specific adjustment
        binsize: Bin size for coverage calculation
        chromsizes: Dictionary of chromosome sizes

    Returns:
        DataFrame with coverage data
    """
    shift_values = np.where(pairs_chunk["strand"] == "+", shift, -shift)
    pairs_chunk[["start", "end"]] = pairs_chunk[["start", "end"]].add(
        shift_values, axis=0
    )
    pairs_chunk = pairs_chunk[
        (pairs_chunk["start"] >= 0)
        & (pairs_chunk["end"] <= (pairs_chunk["chrom"].map(chromsizes)))
    ]
    pairs_chunk = pairs_chunk[["chrom", "start", "end"]]
    pairs_chunk.sort_values(["chrom", "start", "end"], inplace=True)
    pairs_chunk.reset_index(drop=True, inplace=True)
    return coverage_chunk(pairs_chunk, binsize)


def calculate_coverage(
    pairs_chunks: Iterator[pd.DataFrame],
    chromsizes: Union[dict, pd.Series],
    side: str = "1",
    end: str = "0",
    shift: int = 0,
    window_size: int = 10,
    threads: int = 1,
) -> pd.DataFrame:
    """
    Calculate coverage from pairs chunks.

    Args:
        pairs_chunks: Iterator of DataFrames with pairs data
        chromsizes: Dictionary or pandas Series with chromosome sizes
        side: Which side to use (0=both, 1=first, 2=second)
        end: Which end to use (0=default, 3=3', 5=5')
        shift: Shift value for strand-specific adjustment
        window_size: Window size for binning
        threads: Number of threads to use

    Returns:
        DataFrame with coverage data
    """
    chunks = produce_chunks(pairs_chunks, side, end)
    func = partial(
        process_chunk, shift=shift, binsize=window_size, chromsizes=chromsizes
    )

    if threads > 1:
        with multiprocessing.Pool(threads) as pool:
            coverage_chunks = pool.map(func, chunks)
    else:
        coverage_chunks = map(func, chunks)

    coverage_df = (
        pd.concat(coverage_chunks, ignore_index=True)
        .groupby(["chrom", "bin"], as_index=False)
        .agg({"count": "sum"})
        .reset_index(drop=True)
    )
    coverage_df["start"] = (coverage_df["bin"] * window_size).astype(int)
    coverage_df["end"] = coverage_df["start"] + window_size
    coverage_df["end"] = np.where(
        coverage_df["end"] <= coverage_df["chrom"].map(chromsizes),
        coverage_df["end"],
        coverage_df["chrom"].map(chromsizes),
    )

    # Filter out empty intervals
    coverage_df = coverage_df[coverage_df["end"] - coverage_df["start"] > 0]
    coverage_df = coverage_df[["chrom", "start", "end", "count"]]

    return coverage_df


def save_coverage(
    coverage_df: pd.DataFrame,
    output_file: str,
    chromsizes: dict,
    output_bigwig: Optional[str] = None,
) -> None:
    """
    Save coverage data to file(s).

    Args:
        coverage_df: DataFrame with coverage data
        output_file: Path to output tab-separated file
        chromsizes: Dictionary of chromosome sizes
        output_bigwig: Optional path to output BigWig file
    """
