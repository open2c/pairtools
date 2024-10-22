import numpy as np
import pandas as pd

from .regions import assign_regs_c
from . import pairsio
import bioframe


def geomprog(factor, start=1):
    yield start
    while True:
        start *= factor
        yield start


def _geomrange(start, end, factor, endpoint):
    prev = np.nan
    for i in geomprog(factor, start):
        x = int(round(i))

        if x > end:
            break

        if x == prev:
            continue

        prev = x
        yield x

    if endpoint and prev != end:
        yield end


def geomrange(start, end, factor, endpoint=False):
    return np.fromiter(_geomrange(start, end, factor, endpoint), dtype=int)


def geomspace(start, end, num=50, endpoint=True):
    factor = (end / start) ** (1 / num)
    return geomrange(start, end, factor, endpoint=endpoint)


def _to_float(arr_or_scalar):
    if np.isscalar(arr_or_scalar):
        return float(arr_or_scalar)
    else:
        return np.asarray(arr_or_scalar).astype(float)


def assign_regs(chroms, pos, regs):
    gb_regs = regs.sort_values(["chrom", "start", "end"]).groupby(["chrom"])

    regs_dict = {
        chrom.encode(): regs_per_chrom[["start", "end"]]
        .values.flatten()
        .astype(np.int64)
        for chrom, regs_per_chrom in gb_regs
    }

    return assign_regs_c(np.asarray(chroms).astype("bytes"), np.asarray(pos), regs_dict)


def cartesian_df_product(df1, df2, suffixes=["1", "2"]):
    return pd.merge(
        left=df1.assign(cartesian_product_dummy=1),
        right=df2.assign(cartesian_product_dummy=1),
        on=["cartesian_product_dummy"],
        how="outer",
        suffixes=suffixes,
    ).drop("cartesian_product_dummy", axis="columns")


def make_empty_scaling(regions, dist_bins, multiindex=True):

    if dist_bins[0] != 0:
        dist_bins = np.r_[0, dist_bins]
    if dist_bins[-1] != np.iinfo(np.int64).max:
        dist_bins = np.r_[dist_bins, np.iinfo(np.int64).max]

    strands_table = pd.DataFrame(
        {"strand1": ["+", "+", "-", "-"], "strand2": ["+", "-", "+", "-"]}
    )
    dists_table = pd.DataFrame(
        list(zip(dist_bins[:-1], dist_bins[1:])), columns=["min_dist", "max_dist"]
    )

    out = regions.join(regions, on=None, lsuffix="1", rsuffix="2")
    out = cartesian_df_product(out, strands_table)
    out = cartesian_df_product(out, dists_table)

    if multiindex:
        index_by = [
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "strand1",
            "strand2",
            "min_dist",
            "max_dist",
        ]
        out.set_index(index_by, inplace=True)

    return out


def make_empty_cross_region_table(
    regions, drop_same_reg=True, split_by_strand=True, multiindex=True
):
    out = cartesian_df_product(regions, regions)
    if split_by_strand:
        strands_table = pd.DataFrame(
            {"strand1": ["+", "+", "-", "-"], "strand2": ["+", "-", "+", "-"]}
        )
        out = cartesian_df_product(out, strands_table)

    if drop_same_reg:
        out = out[
            (out["chrom1"] != out["chrom2"])
            | (out["start1"] != out["start2"])
            | (out["end1"] != out["end2"])
        ]

    if multiindex:
        index_by = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
        if split_by_strand:
            index_by += ["strand1", "strand2"]

        out.set_index(index_by, inplace=True)

    return out


def bins_pairs_by_distance(
    pairs_df,
    dist_bins,
    regions=None,
    chromsizes=None,
    ignore_trans=False,
    keep_unassigned=False,
):

    dist_bins = np.r_[dist_bins, np.iinfo(np.int64).max]
    if regions is None:
        if chromsizes is None:
            chroms = sorted(
                set.union(set(pairs_df.chrom1.unique()), set(pairs_df.chrom2.unique()))
            )
            regions = pd.DataFrame({"chrom": chroms, "start": 0, "end": -1})
            regions = regions[["chrom", "start", "end"]]

            region_starts1, region_starts2 = 0, 0
            region_ends1, region_ends2 = -1, -1

        else:
            region_ends1 = pairs_df.chrom1.map(chromsizes).fillna(-1).astype(np.int64)
            region_ends2 = pairs_df.chrom2.map(chromsizes).fillna(-1).astype(np.int64)
            region_starts1 = np.where(region_ends1 > 0, 0, -1)
            region_starts2 = np.where(region_ends2 > 0, 0, -1)
            regions = pd.DataFrame(
                [
                    {"chrom": chrom, "start": 0, "end": length}
                    for chrom, length in chromsizes.items()
                ]
            )
            regions = regions[["chrom", "start", "end"]]

        try:
            regions = bioframe.from_any(regions)
        except Exception as e:
            raise ValueError(f"Invalid viewframe created from pairs file, {e}")

    else:

        if not bioframe.is_viewframe(regions):
            try:
                regions = bioframe.from_any(regions)
            except Exception as e:
                raise ValueError(
                    f"Provided regions cannot be converted to viewframe, {e}"
                )

        regions = regions[["chrom", "start", "end"]]

        _, region_starts1, region_ends1 = assign_regs(
            pairs_df.chrom1.values, pairs_df.pos1.values, regions
        ).T
        _, region_starts2, region_ends2 = assign_regs(
            pairs_df.chrom2.values, pairs_df.pos2.values, regions
        ).T

    pairs_reduced_df = pd.DataFrame(
        {
            "chrom1": pairs_df.chrom1.values,
            "start1": region_starts1,
            "end1": region_ends1,
            "chrom2": pairs_df.chrom2.values,
            "start2": region_starts2,
            "end2": region_ends2,
            "strand1": pairs_df.strand1.values,
            "strand2": pairs_df.strand2.values,
            "dist_bin_idx": np.searchsorted(
                dist_bins, np.abs(pairs_df.pos1 - pairs_df.pos2), side="right"
            ),
            "n_pairs": 1,
        },
        copy=False,
    )

    if not keep_unassigned:
        pairs_reduced_df = (
            pairs_reduced_df.query("(start1 >= 0) and (start2 >= 0)")
            # do not test for end1 and end2, as they can be -1 if regions and not specified
            .reset_index(drop=True)
        )

    pairs_reduced_df["min_dist"] = np.where(
        pairs_reduced_df["dist_bin_idx"] > 0,
        dist_bins[pairs_reduced_df["dist_bin_idx"] - 1],
        0,
    )

    pairs_reduced_df["max_dist"] = np.where(
        pairs_reduced_df["dist_bin_idx"] < len(dist_bins) - 1,
        dist_bins[pairs_reduced_df["dist_bin_idx"]],
        np.iinfo(np.int64).max,
    )

    # importantly, in the future, we may want to extend the function to plot scalings
    # for pairs from different regions!

    cis_region_pairs = (
        (pairs_reduced_df.chrom1 == pairs_reduced_df.chrom2)
        & (pairs_reduced_df.start1 == pairs_reduced_df.start2)
        & (pairs_reduced_df.end1 == pairs_reduced_df.end2)
    )

    pairs_for_scaling_mask = (
        cis_region_pairs
        & (pairs_reduced_df.min_dist > 0)
        & (pairs_reduced_df.max_dist < np.iinfo(np.int64).max)
    )

    pairs_for_scaling_df = pairs_reduced_df.loc[pairs_for_scaling_mask]

    pairs_for_scaling_counts = pairs_for_scaling_df.groupby(
        by=[
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "strand1",
            "strand2",
            "min_dist",
            "max_dist",
        ]
    ).agg({"n_pairs": "sum"})

    pairs_for_scaling_counts = (
        make_empty_scaling(regions, dist_bins)
        .assign(n_pairs=0)
        .add(pairs_for_scaling_counts, fill_value=0)
    )
    pairs_for_scaling_counts["n_pairs"] = pairs_for_scaling_counts["n_pairs"].astype(
        np.int64
    )

    if ignore_trans:
        pairs_no_scaling_counts = None
    else:
        pairs_no_scaling_df = pairs_reduced_df.loc[~cis_region_pairs]

        pairs_no_scaling_counts = pairs_no_scaling_df.groupby(
            by=[
                "chrom1",
                "start1",
                "end1",
                "chrom2",
                "start2",
                "end2",
                "strand1",
                "strand2",
            ]
        ).agg({"n_pairs": "sum"})

        pairs_no_scaling_counts = (
            make_empty_cross_region_table(regions)
            .assign(n_pairs=0)
            .add(pairs_no_scaling_counts, fill_value=0)
        )
        pairs_no_scaling_counts["n_pairs"] = pairs_no_scaling_counts["n_pairs"].astype(
            np.int64
        )

    return pairs_for_scaling_counts, pairs_no_scaling_counts


def contact_areas_same_reg(min_dist, max_dist, region_length):

    min_dist = _to_float(min_dist)
    max_dist = _to_float(max_dist)
    scaffold_length = _to_float(region_length)
    outer_areas = np.maximum(region_length - min_dist, 0) ** 2
    inner_areas = np.maximum(region_length - max_dist, 0) ** 2
    return 0.5 * (outer_areas - inner_areas)


def _contact_areas_diff_reg(
    min_dist, max_dist, region_start1, region_end1, region_start2, region_end2
):

    return (
        contact_areas_same_reg(min_dist, max_dist, np.abs(region_end2 - region_start1))
        + contact_areas_same_reg(
            min_dist, max_dist, np.abs(region_end1 - region_start2)
        )
        - contact_areas_same_reg(
            min_dist, max_dist, np.abs(region_start1 - region_start2)
        )
        - contact_areas_same_reg(min_dist, max_dist, np.abs(region_end1 - region_end2))
    )


def _contact_areas_trans(min_dist, max_dist, region_length1, region_length2):

    return (
        contact_areas_same_reg(min_dist, max_dist, region_length1 + region_length2)
        - contact_areas_same_reg(min_dist, max_dist, region_length1)
        - contact_areas_same_reg(min_dist, max_dist, region_length2)
    )


def compute_scaling(
    pairs,
    regions=None,
    chromsizes=None,
    dist_range=(int(1e0), int(1e9)),
    n_dist_bins_decade=8,
    chunksize=int(1e7),
    ignore_trans=False,
    keep_unassigned=False,
    filter_f=None,
    nproc_in=4,
):
    """
    Compute the contact-frequency-vs-distance (aka "scaling") curve from a table of contacts.

    Parameters
    ----------
    pairs : pd.DataFrame or str or file-like object
        A table with pairs of genomic coordinates representing contacts.
        It can be a pandas DataFrame, a path to a pairs file, or a file-like object.
    regions : bioframe viewframe or None, optional
        Genomic regions of interest. It can be anything that can serve as input to bioframe.from_any,
        or None if not applicable.
    chromsizes : pd.DataFrame or None, optional
        Additional dataframe with chromosome sizes, if different from regions.
    dist_range : tuple of int, optional
        The range of distances to calculate the scaling curve. Default is (10, 1000000000).
    n_dist_bins : int, optional
        The number of distance bins per order of magnitude in a log10-space. Default is 8.
    chunksize : int, optional
        Size of chunks for calculations. Default is 10000000.
    ignore_trans : bool, optional
        Ignore trans interactions or not. Default is False.
    keep_unassigned : bool, optional
        Keep pairs that are not assigned to any region or not. Default is False.
    filter_f : function or None, optional
        A function that to filter contacts. Default is None.
    nproc_in : int, optional
        Number of processes to use for reading pairs file. Default is 1.

    Returns
    -------
    sc : pd.DataFrame
        Scaling information for each distance bin.
    trans_counts : pd.DataFrame or None
        Trans interaction counts for each distance bin. None if ignore_trans is True.
    """

    dist_bins = geomspace(
        dist_range[0],
        dist_range[1],
        int(np.round(np.log10(dist_range[1] / dist_range[0]) * n_dist_bins_decade)),
    )

    if isinstance(pairs, pd.DataFrame):
        pairs_df = pairs

    elif isinstance(pairs, str) or hasattr(pairs, "buffer") or hasattr(pairs, "peek"):
        pairs_df, _, chromsizes = pairsio.read_pairs(
            pairs, nproc=nproc_in, chunksize=chunksize
        )
    else:
        raise ValueError(
            "pairs must be either a path to a pairs file or a pd.DataFrame"
        )

    sc, trans_counts = None, None
    for pairs_chunk in [pairs_df] if isinstance(pairs_df, pd.DataFrame) else pairs_df:
        if filter_f:
            pairs_chunk = filter_f(pairs_chunk)
        sc_chunk, trans_counts_chunk = bins_pairs_by_distance(
            pairs_chunk,
            dist_bins,
            regions=regions,
            chromsizes=chromsizes,
            ignore_trans=ignore_trans,
            keep_unassigned=keep_unassigned,
        )

        sc = sc_chunk if sc is None else sc.add(sc_chunk, fill_value=0)

        trans_counts = (
            trans_counts_chunk
            if trans_counts is None
            else trans_counts.add(trans_counts_chunk, fill_value=0)
        )

    #         if not (isinstance(regions, pd.DataFrame) and
    #                  (set(regions.columns) == set(['chrom', 'start','end']))):
    #             raise ValueError('regions must be provided as a dict or chrom-indexed Series of chromsizes or as a bedframe.')

    sc.reset_index(inplace=True)
    sc["n_bp2"] = contact_areas_same_reg(
        sc["min_dist"], sc["max_dist"], sc["end1"] - sc["start1"]
    )

    if not ignore_trans:
        trans_counts.reset_index(inplace=True)
        trans_counts["n_bp2"] = (trans_counts["end1"] - trans_counts["start1"]) * (
            trans_counts["end2"] - trans_counts["start2"]
        )

    return sc, trans_counts


def norm_scaling_factor(bins, cfreqs, norm_window):
    """
    Calculate the normalization factor for a contact-frequency-vs-distance curve,
    by setting the average contact frequency in a specified range of distances to 1.0.

    Args:
        bins (array-like): The distance bins.
        cfreqs (array-like): The contact frequencies.
        norm_window (tuple of float): A tuple with the range of distances to use for normalization.

    Returns:
        float: The normalization scaling factor.
    """

    lo, hi = np.searchsorted(bins, norm_window)
    return cfreqs[lo : hi + 1].mean()


def norm_scaling(bins, cfreqs, norm_window, log_input=False):
    """
    Normalize a contact-frequency-vs-distance curve, by setting the average contact frequency
    in a given window to 1.0.

    Args:
        bins (array-like): The distance bins.
        cfreqs (array-like): The contact frequencies.
        norm_window (tuple of float): A tuple with the range of distances to use for normalization.
        log_input (bool, optional): Whether the input contact frequencies were log-transformed. Defaults to False.

    Returns:
        float or array-like: The normalized contact frequencies.
    """

    norm = norm_scaling_factor(bins, cfreqs, norm_window)
    if log_input:
        return cfreqs - norm
    else:
        return cfreqs / norm


def unity_norm_scaling(bins, cfreqs, norm_range=(1e4, 1e9)):
    bin_lens = np.diff(bins)
    bin_mids = np.sqrt(bins[1:] * bins[:-1])

    if norm_range is None:
        norm_cfreqs = cfreqs / np.sum(1.0 * (bin_lens * cfreqs)[np.isfinite(cfreqs)])
    else:
        norm_cfreqs = cfreqs / np.sum(
            1.0
            * (bin_lens * cfreqs)[
                np.isfinite(cfreqs)
                & (bin_mids > norm_range[0])
                & (bin_mids < norm_range[1])
            ]
        )

    return norm_cfreqs
