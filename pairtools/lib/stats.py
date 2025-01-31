import numpy as np
import pandas as pd
from scipy import special
from collections.abc import Mapping
import sys
import yaml
from . import fileio
from .select import evaluate_df

from .._logging import get_logger

logger = get_logger()


def parse_number(s):
    if s.isdigit():
        return int(s)
    elif s.replace(".", "", 1).isdigit():
        return float(s)
    else:
        return s


def flat_dict_to_nested(input_dict, sep="/"):
    output_dict = {}

    for key, value in input_dict.items():
        if type(key) == tuple:
            key_parts = key
        elif type(key) == str:
            key_parts = key.split(sep)
        else:
            raise ValueError(
                f"Key type can be either str or tuple. Found key {key} of type {type(key)}."
            )

        current_dict = output_dict
        for key_part in key_parts[:-1]:
            current_dict = current_dict.setdefault(key_part, {})
        current_dict[key_parts[-1]] = value

    return output_dict


def nested_dict_to_flat(d, tuple_keys=False, sep="/"):
    """Flatten a nested dictionary to a flat dictionary.

    Parameters
    ----------
    d: dict
        A nested dictionary to flatten.
    tuple_keys: bool
        If True, keys will be joined into tuples. Otherwise, they will be joined into strings.
    sep: str
        The separator to use between the parent key and the key if tuple_keys==False.
    Returns
    -------
    dict
        A flat dictionary.
    """

    if tuple_keys:
        join_keys = lambda k1, k2: (k1,) + k2
    else:
        join_keys = lambda k1, k2: (k1 + sep + k2) if k2 else k1

    out = {}

    for k1, v1 in d.items():
        if isinstance(v1, dict):
            out.update(
                {
                    join_keys(k1, k2): v2
                    for k2, v2 in nested_dict_to_flat(v1, tuple_keys, sep).items()
                }
            )
        else:
            if tuple_keys:
                out[(k1,)] = v1
            else:
                out[k1] = v1

    return out


def is_nested_dict(d):
    """Check if a dictionary is nested.

    Parameters
    ----------
    d: dict
        A dictionary to check.
    Returns
    -------
    bool
        True if the dictionary is nested, False otherwise.
    """

    if not isinstance(d, dict):
        return False

    for v in d.values():
        if isinstance(v, dict):
            return True

    return False


def is_tuple_keyed_dict(d):
    """Check if a dictionary is tuple-keyed.

    Parameters
    ----------
    d: dict
        A dictionary to check.
    Returns
    -------
    bool
        True if the dictionary is tuple-keyed, False otherwise.
    """

    if not isinstance(d, dict):
        return False

    for k, v in d.items():
        if not isinstance(k, tuple):
            return False
        if isinstance(v, dict):
            return False

    return True


def is_str_keyed_dict(d):
    """Check if a dictionary is string-keyed.

    Parameters
    ----------
    d: dict
        A dictionary to check.
    Returns
    -------
    bool
        True if the dictionary is string-keyed, False otherwise.
    """

    if not isinstance(d, dict):
        return False

    for k, v in d.keys():
        if not isinstance(k, str):
            return False
        if isinstance(v, dict):
            return False

    return True


def swap_levels_nested_dict(nested_dict, level1, level2, sep="/"):
    """Swap the order of two levels in a nested dictionary.

    Parameters
    ----------
    nested_dict: dict
        A nested dictionary.
    level1: int
        The index of the first level to swap.
    level2: int
        The index of the second level to swap.
    Returns
    -------
    dict
        A nested dictionary with the levels swapped.
    """

    if is_tuple_keyed_dict(nested_dict):
        out = {}
        for k1, v1 in nested_dict.items():
            k1_list = list(k1)
            k1_list[level1], k1_list[level2] = k1_list[level2], k1_list[level1]
            out[tuple(k1_list)] = v1
        return out

    elif is_nested_dict(nested_dict):
        out = nested_dict_to_flat(nested_dict, tuple_keys=True)
        out = swap_levels_nested_dict(out, level1, level2)
        out = flat_dict_to_nested(out)
        return out

    elif is_str_keyed_dict(nested_dict):
        out = nested_dict_to_flat(nested_dict, sep=sep)
        out = swap_levels_nested_dict(out, level1, level2)
        out = {sep.join(k): v for k, v in out.items()}
        return out

    else:
        raise ValueError(
            "Input dictionary must be either nested, string-keyed or tuple-keyed"
        )


class PairCounter(Mapping):
    """
    A Counter for Hi-C pairs that accumulates various statistics.
    PairCounter implements two interfaces to access multi-level statistics:
    1. as a nested dict, e.g. pairCounter['pair_types']['LL']
    2. as a flat dict, with the level keys separated by '/', e.g. pairCounter['pair_types/LL']
    Other features:
    -- PairCounters can be saved into/loaded from a file
    -- multiple PairCounters can be merged via addition.
    """

    _SEP = "\t"
    _KEY_SEP = "/"
    DIST_FREQ_REL_DIFF_THRESHOLD = 0.05
    N_DIST_BINS_DECADE_DEFAULT = 8
    MIN_LOG10_DIST_DEFAULT = 0
    MAX_LOG10_DIST_DEFAULT = 9

    def __init__(
        self,
        min_log10_dist=MIN_LOG10_DIST_DEFAULT,
        max_log10_dist=MAX_LOG10_DIST_DEFAULT,
        n_dist_bins_decade=N_DIST_BINS_DECADE_DEFAULT,
        bytile_dups=False,
        filters=None,
        **kwargs,
    ):
        # Define filters and parameters for filters evaluation:
        if filters is not None:
            self.filters = filters
        else:
            self.filters = {"no_filter": ""}
        self.startup_code = kwargs.get("startup_code", "")
        self.type_cast = kwargs.get("type_cast", ())
        self.engine = kwargs.get("engine", "pandas")

        # Define default filter:
        if "no_filter" not in self.filters:
            self.filters["no_filter"] = ""
        self._stat = {key: {} for key in self.filters}

        # some variables used for initialization:
        # genomic distance bining for the ++/--/-+/+- distribution
        log10_dist_bin_step = 1.0 / n_dist_bins_decade
        self._dist_bins = np.unique(
            np.r_[
                0,
                np.round(
                    10
                    ** np.arange(
                        min_log10_dist, max_log10_dist + 0.001, log10_dist_bin_step
                    )
                ).astype(np.int_),
            ]
        )

        # establish structure of an empty _stat:
        for key in self.filters:
            self._stat[key]["filter_expression"] = self.filters[key]
            self._stat[key]["total"] = 0
            self._stat[key]["total_unmapped"] = 0
            self._stat[key]["total_single_sided_mapped"] = 0
            # total_mapped = total_dups + total_nodups
            self._stat[key]["total_mapped"] = 0
            self._stat[key]["total_dups"] = 0
            self._stat[key]["total_nodups"] = 0
            ########################################
            # the rest of stats are based on nodups:
            ########################################
            self._stat[key]["cis"] = 0
            self._stat[key]["trans"] = 0
            self._stat[key]["pair_types"] = {}

            # to be removed:
            # self._stat[key]["dedup"] = {}

            self._stat[key]["cis_1kb+"] = 0
            self._stat[key]["cis_2kb+"] = 0
            self._stat[key]["cis_4kb+"] = 0
            self._stat[key]["cis_10kb+"] = 0
            self._stat[key]["cis_20kb+"] = 0
            self._stat[key]["cis_40kb+"] = 0
            self._stat[key]["summary"] = dict(
                [
                    ("frac_cis", 0),
                    ("frac_cis_1kb+", 0),
                    ("frac_cis_2kb+", 0),
                    ("frac_cis_4kb+", 0),
                    ("frac_cis_10kb+", 0),
                    ("frac_cis_20kb+", 0),
                    ("frac_cis_40kb+", 0),
                    ("frac_dups", 0),
                    ("complexity_naive", 0),
                ]
            )

            self._stat[key]["chrom_freq"] = {}

            self._stat[key]["dist_freq"] = {
                "+-": {bin.item(): 0 for bin in self._dist_bins},
                "-+": {bin.item(): 0 for bin in self._dist_bins},
                "--": {bin.item(): 0 for bin in self._dist_bins},
                "++": {bin.item(): 0 for bin in self._dist_bins},
            }

            self._stat[key]["chromsizes"] = {}

            # Summaries are derived from other stats and are recalculated on merge

        self._save_bytile_dups = bytile_dups
        if self._save_bytile_dups:
            self._bytile_dups = pd.DataFrame(
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["tile", "parent_tile"]
                )
            )
        self._summaries_calculated = False

    def __getitem__(self, key, filter="no_filter"):
        if isinstance(key, str):
            # let's strip any unintentional '/'
            # from either side of the key
            key = key.strip("/")
            if self._KEY_SEP in key:
                # multi-key to access nested elements
                k_fields = key.split(self._KEY_SEP)
            else:
                # single-key access flat part of PairCounter
                # or to access highest level of hierarchy
                return self._stat[filter][key]
        else:
            # clearly an error:
            raise ValueError("{} is not a valid key: must be str".format(key))

        # K_FIELDS:
        # process multi-key case:
        # in this case key must be in ['pair_types','chrom_freq','dist_freq','dedup']
        # get the first 'k' and keep the remainders in 'k_fields'
        k = k_fields.pop(0)
        if k in ["pair_types", "dedup"]:
            # assert there is only one element in key_fields left:
            # 'pair_types' and 'dedup' treated the same
            if len(k_fields) == 1:
                return self._stat[filter][k][k_fields[0]]
            else:
                raise ValueError(
                    "{} is not a valid key: {} section implies 1 identifier".format(
                        key, k
                    )
                )

        elif k == "chrom_freq":
            # assert remaining key_fields == [chr1, chr2]:
            if len(k_fields) == 2:
                return self._stat[filter][k][tuple(k_fields)]
            else:
                raise ValueError(
                    "{} is not a valid key: {} section implies 2 identifiers".format(
                        key, k
                    )
                )

        elif k == "dist_freq":
            # assert that last element of key_fields is the 'directions'
            # THIS IS DONE FOR CONSISTENCY WITH .stats FILE
            # SHOULD THAT BE CHANGED IN .stats AND HERE AS WELL?
            if len(k_fields) == 2:
                # assert 'dirs' in ['++','--','+-','-+']
                dirs = k_fields.pop()
                # there is only genomic distance range of the bin that's left:
                (bin_range,) = k_fields
                # extract left border of the bin "1000000+" or "1500-6000":
                dist_bin_left = int(
                    bin_range.strip("+")
                    if bin_range.endswith("+")
                    else bin_range.split("-")[0]
                )
                # store corresponding value:
                return self._stat[filter]["dist_freq"][dirs][dist_bin_left]
            else:
                raise ValueError(
                    "{} is not a valid key: {} section implies 2 identifiers".format(
                        key, k
                    )
                )
        else:
            raise ValueError("{} is not a valid key".format(k))

    def __iter__(self):
        return iter(self._stat)

    def __len__(self):
        return len(self._stat)

    def find_dist_freq_convergence_distance(self, rel_threshold):
        """Finds the largest distance at which the frequency of pairs of reads
        with different strands deviates from their average by the specified
        relative threshold."""

        out = {}
        all_strands = ["++", "--", "-+", "+-"]

        for filter in self.filters:
            out[filter] = {}

            dist_freqs_by_strands = {
                strands: np.array(
                    list(self._stat[filter]["dist_freq"][strands].values())
                )
                for strands in all_strands
            }

            # Calculate the average frequency of pairs with different strands
            avg_freq_all_strands = np.mean(
                np.vstack(list(dist_freqs_by_strands.values())), axis=0
            )

            # Calculate the largest distance at which the frequency of pairs of at least one strand combination deviates from the average by the given threshold
            rel_deviations = {
                strands: np.nan_to_num(
                    np.abs(dist_freqs_by_strands[strands] - avg_freq_all_strands)
                    / avg_freq_all_strands
                )
                for strands in all_strands
            }

            idx_maxs = {strand: 0 for strand in all_strands}
            for strands in all_strands:
                bin_exceeds = rel_deviations[strands] > rel_threshold
                if np.any(bin_exceeds):
                    idx_maxs[strands] = np.max(np.nonzero(bin_exceeds))

            # Find the largest distance and the strand combination where frequency of pairs deviates from the average by the given threshold:
            convergence_bin_idx = 0
            convergence_strands = "??"
            convergence_dist = "0"

            for strands in all_strands:
                if idx_maxs[strands] > convergence_bin_idx:
                    convergence_bin_idx = idx_maxs[strands]
                    convergence_strands = strands

                    if idx_maxs[strands] < len(self._dist_bins):
                        convergence_dist = self._dist_bins[convergence_bin_idx + 1]
                    else:
                        convergence_dist = np.iinfo(np.int64)

            out[filter]["convergence_dist"] = convergence_dist
            out[filter]["strands_w_max_convergence_dist"] = convergence_strands
            out[filter]["convergence_rel_diff_threshold"] = rel_threshold

            out[filter]["n_cis_pairs_below_convergence_dist"] = {
                strands: dist_freqs_by_strands[strands][: convergence_bin_idx + 1].sum()
                for strands in all_strands
                for strands in all_strands
            }

            out[filter]["n_cis_pairs_below_convergence_dist_all_strands"] = sum(
                list(out[filter]["n_cis_pairs_below_convergence_dist"].values())
            )

            n_cis_pairs_above_convergence_dist = {
                strands: dist_freqs_by_strands[strands][convergence_bin_idx + 1 :].sum()
                for strands in all_strands
                for strands in all_strands
            }

            out[filter]["n_cis_pairs_above_convergence_dist_all_strands"] = sum(
                list(n_cis_pairs_above_convergence_dist.values())
            )

            norms = dict(
                cis=self._stat[filter]["cis"],
                total_mapped=self._stat[filter]["total_mapped"],
            )

            if "total_nodups" in self._stat[filter]:
                norms["total_nodups"] = self._stat[filter]["total_nodups"]

            for key, norm_factor in norms.items():
                out[filter][f"frac_{key}_in_cis_below_convergence_dist"] = {
                    strands: n_cis_pairs / norm_factor
                    for strands, n_cis_pairs in out[filter][
                        "n_cis_pairs_below_convergence_dist"
                    ].items()
                }

                out[filter][f"frac_{key}_in_cis_below_convergence_dist_all_strands"] = (
                    sum(
                        list(
                            out[filter][
                                f"frac_{key}_in_cis_below_convergence_dist"
                            ].values()
                        )
                    )
                )

                out[filter][f"frac_{key}_in_cis_above_convergence_dist_all_strands"] = (
                    sum(list(n_cis_pairs_above_convergence_dist.values())) / norm_factor
                )

        return out

    def calculate_summaries(self):
        """calculate summary statistics (fraction of cis pairs at different cutoffs,
        complexity estimate) based on accumulated counts. Results are saved into
        self._stat["filter_name"]['summary"]
        """
        convergence_stats = self.find_dist_freq_convergence_distance(
            self.DIST_FREQ_REL_DIFF_THRESHOLD
        )

        for filter_name in self.filters.keys():
            for cis_count in (
                "cis",
                "cis_1kb+",
                "cis_2kb+",
                "cis_4kb+",
                "cis_10kb+",
                "cis_20kb+",
                "cis_40kb+",
            ):
                self._stat[filter_name]["summary"][f"frac_{cis_count}"] = (
                    (
                        self._stat[filter_name][cis_count]
                        / self._stat[filter_name]["total_nodups"]
                    )
                    if self._stat[filter_name]["total_nodups"] > 0
                    else 0
                )

            self._stat[filter_name]["summary"]["dist_freq_convergence"] = (
                convergence_stats[filter_name]
            )

            self._stat[filter_name]["summary"]["frac_dups"] = (
                (
                    self._stat[filter_name]["total_dups"]
                    / self._stat[filter_name]["total_mapped"]
                )
                if self._stat[filter_name]["total_mapped"] > 0
                else 0
            )

            self._stat[filter_name]["summary"]["complexity_naive"] = (
                estimate_library_complexity(
                    self._stat[filter_name]["total_mapped"],
                    self._stat[filter_name]["total_dups"],
                    0,
                )
            )

            if filter_name == "no_filter" and self._save_bytile_dups:
                # Estimate library complexity with information by tile, if provided:
                if self._bytile_dups.shape[0] > 0:
                    self._stat[filter_name]["dups_by_tile_median"] = int(
                        round(
                            self._bytile_dups["dup_count"].median()
                            * self._bytile_dups.shape[0]
                        )
                    )
                if "dups_by_tile_median" in self._stat[filter_name].keys():
                    self._stat[filter_name]["summary"][
                        "complexity_dups_by_tile_median"
                    ] = estimate_library_complexity(
                        self._stat[filter_name]["total_mapped"],
                        self._stat[filter_name]["total_dups"],
                        self._stat[filter_name]["total_dups"]
                        - self._stat[filter_name]["dups_by_tile_median"],
                    )

            self._summaries_calculated = True

    @classmethod
    def from_file(cls, file_handle, n_dist_bins_decade=N_DIST_BINS_DECADE_DEFAULT):
        """create instance of PairCounter from file
        Parameters
        ----------
        file_handle: file handle
        Returns
        -------
        PairCounter
            new PairCounter filled with the contents of the input file
        """
        # fill in from file - file_handle:
        default_filter = "no_filter"
        stat_from_file = cls(n_dist_bins_decade=n_dist_bins_decade)
        raw_stat = {}
        for l in file_handle:
            key_val_pair = l.strip().split(cls._SEP)
            if len(key_val_pair) == 0:
                # skip empty lines:
                continue
            if len(key_val_pair) != 2:
                # expect two _SEP separated values per line:
                raise fileio.ParseError(
                    "{} is not a valid stats file".format(file_handle.name)
                )
            raw_stat[key_val_pair[0]] = parse_number(key_val_pair[1])

        ## TODO: check if raw_stat does not contain any unknown keys

        # Convert flat dict to nested dict
        stat_from_file._stat[default_filter].update(
            flat_dict_to_nested(raw_stat, sep=cls._KEY_SEP)
        )

        stat_from_file._stat[default_filter]["chrom_freq"] = nested_dict_to_flat(
            stat_from_file._stat[default_filter]["chrom_freq"], tuple_keys=True
        )

        bin_to_left_val = lambda bin: int(
            bin.rstrip("+") if ("+" in bin) else bin.split("-")[0]
        )

        stat_from_file._stat[default_filter]["dist_freq"] = {
            bin_to_left_val(k): v
            for k, v in stat_from_file._stat[default_filter]["dist_freq"].items()
        }

        stat_from_file._stat[default_filter]["dist_freq"] = swap_levels_nested_dict(
            stat_from_file._stat[default_filter]["dist_freq"], 0, 1
        )

        return stat_from_file

    @classmethod
    def from_yaml(cls, file_handle, n_dist_bins_decade=N_DIST_BINS_DECADE_DEFAULT):
        """create instance of PairCounter from file
        Parameters
        ----------
        file_handle: file handle
        Returns
        -------
        PairCounter
            new PairCounter filled with the contents of the input file
        """
        # fill in from file - file_handle:
        stat = yaml.safe_load(file_handle)
        stat_from_file = cls(
            n_dist_bins_decade=n_dist_bins_decade,
            filters={
                key: val.get("filter_expression", "") for key, val in stat.items()
            },
        )

        for key, filter in stat.items():
            chromdict = {}
            for chroms in stat[key]["chrom_freq"].keys():
                chromdict[tuple(chroms.split(cls._KEY_SEP))] = stat[key]["chrom_freq"][
                    chroms
                ]
            stat[key]["chrom_freq"] = chromdict
        stat_from_file._stat = stat
        return stat_from_file

    def add_pair(
        self,
        chrom1,
        pos1,
        strand1,
        chrom2,
        pos2,
        strand2,
        pair_type,
        unmapped_chrom="!",
        filter="no_filter",
    ):
        """Gather statistics for a Hi-C pair and add to the PairCounter.
        Parameters
        ----------
        chrom1: str
            chromosome of the first read
        pos1: int
            position of the first read
        strand1: str
            strand of the first read
        chrom2: str
            chromosome of the first read
        pos2: int
            position of the first read
        strand2: str
            strand of the first read
        pair_type: str
            type of the mapped pair of reads
        unmapped_chrom: str
            what string denotes chromosomes in unmapped pairs (default: "!")
        filter: str
            name of the filter toward which the pair should count (default: "no_filter")
        """

        self._stat[filter]["total"] += 1
        # collect pair type stats including DD:
        self._stat[filter]["pair_types"][pair_type] = (
            self._stat[filter]["pair_types"].get(pair_type, 0) + 1
        )
        if chrom1 == unmapped_chrom and chrom2 == unmapped_chrom:
            self._stat[filter]["total_unmapped"] += 1
        elif chrom1 != unmapped_chrom and chrom2 != unmapped_chrom:
            self._stat[filter]["total_mapped"] += 1
            # only mapped ones can be duplicates:
            if pair_type == "DD":
                self._stat[filter]["total_dups"] += 1
            else:
                self._stat[filter]["total_nodups"] += 1
                self._stat[filter]["chrom_freq"][(chrom1, chrom2)] = (
                    self._stat[filter]["chrom_freq"].get((chrom1, chrom2), 0) + 1
                )

                if chrom1 == chrom2:
                    self._stat[filter]["cis"] += 1
                    dist = np.abs(pos2 - pos1)
                    dist_bin = self._dist_bins[
                        np.searchsorted(self._dist_bins, dist, "right") - 1
                    ]
                    self._stat[filter]["dist_freq"][strand1 + strand2][dist_bin] += 1

                    for dist_kb in [1, 2, 4, 10, 20, 40]:
                        if dist >= dist_kb * 1000:
                            self._stat[filter][f"cis_{dist_kb}kb+"] += 1

                else:
                    self._stat[filter]["trans"] += 1
        else:
            self._stat[filter]["total_single_sided_mapped"] += 1

    def add_pairs_from_dataframe(self, df, unmapped_chrom="!"):
        """Gather statistics for Hi-C pairs in a dataframe and add to the PairCounter.

        Parameters
        ----------
        df: pd.DataFrame
            DataFrame with pairs. Needs to have columns:
                'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'pair_type'
        """
        for key in self.filters.keys():
            if key == "no_filter":
                df_filtered = df.copy()
            else:
                condition = self.filters[key]
                filter_passed = evaluate_df(
                    df,
                    condition,
                    type_cast=self.type_cast,
                    startup_code=self.startup_code,
                    engine=self.engine,
                )
                df_filtered = df.loc[filter_passed, :].reset_index(drop=True)
            total_count = df_filtered.shape[0]
            self._stat[key]["total"] += total_count

            # collect pair type stats including DD:
            for pair_type, type_count in (
                df_filtered["pair_type"].value_counts().items()
            ):
                self._stat[key]["pair_types"][pair_type] = (
                    self._stat[key]["pair_types"].get(pair_type, 0) + type_count
                )

            # Count the unmapped by the "unmapped" chromosomes (debatable, as WW are also marked as ! and they might be mapped):
            unmapped_count = np.logical_and(
                df_filtered["chrom1"] == unmapped_chrom,
                df_filtered["chrom2"] == unmapped_chrom,
            ).sum()
            self._stat[key]["total_unmapped"] += int(unmapped_count)

            # Count the mapped:
            df_mapped = df_filtered.loc[
                (df_filtered["chrom1"] != unmapped_chrom)
                & (df_filtered["chrom2"] != unmapped_chrom),
                :,
            ]
            mapped_count = df_mapped.shape[0]

            self._stat[key]["total_mapped"] += mapped_count
            self._stat[key]["total_single_sided_mapped"] += int(
                total_count - (mapped_count + unmapped_count)
            )

            # Count the duplicates:
            if "duplicate" in df_mapped.columns:
                mask_dups = df_mapped["duplicate"]
            else:
                mask_dups = df_mapped["pair_type"] == "DD"
            df_dups = df_mapped[mask_dups]
            dups_count = df_dups.shape[0]
            self._stat[key]["total_dups"] += int(dups_count)
            self._stat[key]["total_nodups"] += int(mapped_count - dups_count)

            df_nodups = df_mapped.loc[~mask_dups, :]
            mask_cis = df_nodups["chrom1"] == df_nodups["chrom2"]
            df_cis = df_nodups.loc[mask_cis, :].copy()

            # Count pairs per chromosome:
            for (chrom1, chrom2), chrom_count in (
                df_nodups[["chrom1", "chrom2"]].value_counts().items()
            ):
                self._stat[key]["chrom_freq"][(chrom1, chrom2)] = (
                    self._stat[key]["chrom_freq"].get((chrom1, chrom2), 0) + chrom_count
                )

            # Count cis-trans by pairs:

            self._stat[key]["cis"] += df_cis.shape[0]
            self._stat[key]["trans"] += df_nodups.shape[0] - df_cis.shape[0]

            # Count cis distance frequencies:
            dist = np.abs(df_cis["pos2"].values - df_cis["pos1"].values)

            df_cis.loc[:, "bin_idx"] = (
                np.searchsorted(self._dist_bins, dist, "right") - 1
            )

            for (strand1, strand2, bin_id), strand_bin_count in (
                df_cis[["strand1", "strand2", "bin_idx"]].value_counts().items()
            ):
                self._stat[key]["dist_freq"][strand1 + strand2][
                    self._dist_bins[bin_id].item()
                ] += strand_bin_count

            self._stat[key]["cis_1kb+"] += int(np.sum(dist >= 1000))
            self._stat[key]["cis_2kb+"] += int(np.sum(dist >= 2000))
            self._stat[key]["cis_4kb+"] += int(np.sum(dist >= 4000))
            self._stat[key]["cis_10kb+"] += int(np.sum(dist >= 10000))
            self._stat[key]["cis_20kb+"] += int(np.sum(dist >= 20000))
            self._stat[key]["cis_40kb+"] += int(np.sum(dist >= 40000))

            ### Add by-tile dups
            if key == "no_filter" and self._save_bytile_dups and (df_dups.shape[0] > 0):
                bytile_dups = analyse_bytile_duplicate_stats(df_dups)
                self._bytile_dups = self._bytile_dups.add(
                    bytile_dups, fill_value=0
                ).astype(int)

    def add_chromsizes(self, chromsizes):
        """Add chromsizes field to the output stats
        Parameters
        ----------
        chromsizes: Dataframe with chromsizes, read by headerops.chromsizes
        """

        chromsizes = chromsizes.to_dict()
        for filter in self._stat.keys():
            self._stat[filter]["chromsizes"] = chromsizes
        return

    def __add__(self, other, filter="no_filter"):
        # both PairCounter are implied to have a list of common fields:
        #
        # 'total', 'total_unmapped', 'total_single_sided_mapped', 'total_mapped',
        # 'cis', 'trans', 'pair_types', 'cis_1kb+', 'cis_2kb+',
        # 'cis_10kb+', 'cis_20kb+', 'chrom_freq', 'dist_freq', 'dedup'
        #
        # If 'chromsizes' are present, they must be identical
        #
        # initialize empty PairCounter for the result of summation:
        sum_stat = PairCounter()
        # use the empty PairCounter to iterate over:
        for k, v in sum_stat._stat[filter].items():
            if k != "chromsizes" and (
                (k not in self._stat[filter]) or (k not in other._stat[filter])
            ):
                # Skip any missing fields and warn
                logger.warning(
                    f"{k} not found in at least one of the input stats, skipping"
                )
                continue

            # not nested fields are summed trivially:
            if isinstance(v, int):
                sum_stat._stat[filter][k] = (
                    self._stat[filter][k] + other._stat[filter][k]
                )

            # sum nested dicts/arrays in a context dependet manner:
            else:
                if k in ["pair_types", "dedup"]:
                    # handy function for summation of a pair of dicts:
                    # https://stackoverflow.com/questions/10461531/merge-and-sum-of-two-dictionaries
                    sum_dicts = lambda dict_x, dict_y: {
                        key: dict_x.get(key, 0) + dict_y.get(key, 0)
                        for key in set(dict_x) | set(dict_y)
                    }
                    # sum a pair of corresponding dicts:
                    sum_stat._stat[filter][k] = sum_dicts(
                        self._stat[filter][k], other._stat[filter][k]
                    )

                elif k == "chrom_freq":
                    # union list of keys (chr1,chr2) with potential duplicates:
                    union_keys_with_dups = list(self._stat[filter][k].keys()) + list(
                        other._stat[filter][k].keys()
                    )
                    # dict.fromkeys will take care of keys' order and duplicates in a consistent manner:
                    # https://stackoverflow.com/questions/1720421/how-to-concatenate-two-lists-in-python
                    # last comment to the 3rd Answer
                    sum_stat._stat[filter][k] = dict.fromkeys(union_keys_with_dups)
                    # perform a summation:
                    for union_key in sum_stat._stat[filter][k]:
                        sum_stat._stat[filter][k][union_key] = self._stat[filter][
                            k
                        ].get(union_key, 0) + other._stat[filter][k].get(union_key, 0)

                elif k == "dist_freq":
                    for dirs in sum_stat[k]:
                        from functools import reduce

                        def reducer(accumulator, element):
                            for key, value in element.items():
                                accumulator[key] = accumulator.get(key, 0) + value
                            return accumulator

                        sum_stat[k][dirs] = reduce(
                            reducer,
                            [self._stat[filter][k][dirs], other._stat[filter][k][dirs]],
                            {},
                        )
                        # sum_stat[k][dirs] = self._stat[filter][k][dirs] + other._stat[filter][k][dirs]

                elif k == "chromsizes":
                    if k in self._stat[filter] and k in other._stat[filter]:
                        if self._stat[filter][k] == other._stat[filter][k]:
                            sum_stat._stat[filter][k] = self._stat[filter][k]
                        elif (
                            len(self._stat[filter][k]) == 0
                            or len(other._stat[filter][k]) == 0
                        ):
                            logger.warning(
                                "One of the stats has no chromsizes recorded,"
                                "writing the one that is present to the output"
                            )
                            if len(self._stat[filter][k]) > 0:
                                sum_stat._stat[filter][k] = self._stat[filter][k]
                            else:
                                sum_stat._stat[filter][k] = other._stat[filter][k]
                        else:
                            raise ValueError(
                                "Can't merge stats with different chromsizes"
                            )
                    else:
                        logger.warning(
                            "One or both stats don't have chromsizes recorded"
                        )

        sum_stat.calculate_summaries()

        return sum_stat

    # we need this to be able to sum(list_of_PairCounters)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def flatten(self, filter="no_filter"):
        """return a flattened dict (formatted same way as .stats file)
        Performed for a single filter."""
        # dict for flat store:
        flat_stat = {}

        # Storing statistics
        for k, v in self._stat[filter].items():
            if isinstance(v, int):
                flat_stat[k] = v
            # store nested dicts/arrays in a context dependent manner:
            # nested categories are stored only if they are non-trivial
            else:
                if (k == "dist_freq") and v:
                    for i in range(len(self._dist_bins)):
                        for dirs, freqs in v.items():
                            dist = self._dist_bins[i]

                            # last bin is treated differently: "100000+" vs "1200-3000":
                            if i < len(self._dist_bins) - 1:
                                dist_next = self._dist_bins[i + 1]
                                formatted_key = self._KEY_SEP.join(
                                    ["{}", "{}-{}", "{}"]
                                ).format(k, dist, dist_next, dirs)
                            elif i == len(self._dist_bins) - 1:
                                formatted_key = self._KEY_SEP.join(
                                    ["{}", "{}+", "{}"]
                                ).format(k, dist, dirs)
                            else:
                                raise ValueError(
                                    "There is a mismatch between dist_freq bins in the instance"
                                )

                            # store key,value pair:
                            try:
                                flat_stat[formatted_key] = freqs[dist]
                            except:
                                # in some previous versions of stats, last bin was not reported, so we need to skip it now:
                                if (dist not in freqs) and (
                                    i == len(self._dist_bins) - 1
                                ):
                                    flat_stat[formatted_key] = 0
                                else:
                                    raise ValueError(
                                        f"Error in {k} {dirs} {dist} {dist_next} {freqs}: source and destination bins do not match"
                                    )

                elif (k in ["pair_types", "dedup", "chromsizes", "summary"]) and v:
                    # 'pair_types' and 'dedup' are simple dicts inside,
                    # treat them the exact same way:
                    flat_stat.update(
                        {
                            k + self._KEY_SEP + k2: v2
                            for k2, v2 in nested_dict_to_flat(
                                v, sep=self._KEY_SEP
                            ).items()
                        }
                    )

                elif (k == "chrom_freq") and v:
                    for (chrom1, chrom2), freq in v.items():
                        formatted_key = self._KEY_SEP.join(["{}", "{}", "{}"]).format(
                            k, chrom1, chrom2
                        )
                        # store key,value pair:
                        flat_stat[formatted_key] = freq

        # return flattened dict
        return flat_stat

    def format_yaml(self, filter="no_filter"):
        """return a formatted dict (for the yaml output)
        Performed for all filters at once."""

        from copy import deepcopy

        formatted_stat = {filter_name: {} for filter_name in self.filters.keys()}

        # Storing statistics for each filter
        for filter_name in self.filters.keys():
            for k, v in self._stat[filter_name].items():
                if k == "chrom_freq":
                    v = {self._KEY_SEP.join(k2): v2 for k2, v2 in v.items()}
                if v:
                    formatted_stat[filter_name][k] = deepcopy(v)
            # return formatted dict
        formatted_stat = nested_dict_to_flat(formatted_stat, tuple_keys=True)
        for k in formatted_stat:
            v = formatted_stat[k]
            if isinstance(v, np.generic):
                formatted_stat[k] = v.item()
        formatted_stat = flat_dict_to_nested(formatted_stat)

        return formatted_stat

    def save(self, outstream, yaml=False, filter="no_filter"):
        """save PairCounter to tab-delimited text file.
        Flattened version of PairCounter is stored in the file.
        Parameters
        ----------
        outstream: file handle
        yaml: is output in yaml format or table
        filter: filter to output in tsv mode
        Note
        ----
        The order of the keys is not guaranteed
        Merging several .stats is not associative with respect to key order:
        merge(A,merge(B,C)) != merge(merge(A,B),C).
        Theys shou5ld match exactly, however, when soprted:
        sort(merge(A,merge(B,C))) == sort(merge(merge(A,B),C))
        """

        if not self._summaries_calculated:
            self.calculate_summaries()

        # write flattened version of the PairCounter to outstream,
        # will output all the filters
        if yaml:
            import yaml

            data = self.format_yaml()
            yaml.dump(data, outstream, default_flow_style=False, sort_keys=False)
        else:  # will output a single filter
            data = self.flatten(filter=filter)
            for k, v in data.items():
                outstream.write("{}{}{}\n".format(k, self._SEP, v))

    def save_bytile_dups(self, outstream):
        """save bytile duplication counts to a tab-delimited text file.
        Parameters
        ----------
        outstream: file handle
        """
        if self._save_bytile_dups:
            self._bytile_dups.reset_index().to_csv(outstream, sep="\t", index=False)
        else:
            logger.error("Bytile dups are not calculated, cannot save.")

    def __repr__(self):
        return str(self._stat)


##################
# Other functions:


def do_merge(output, files_to_merge, **kwargs):
    # Parse all stats files.
    stats = []
    for stat_file in files_to_merge:
        f = fileio.auto_open(
            stat_file,
            mode="r",
            nproc=kwargs.get("nproc_in"),
            command=kwargs.get("cmd_in", None),
        )
        # use a factory method to instanciate PairCounter
        if kwargs.get("yaml", False):
            stat = PairCounter.from_yaml(
                f,
                n_dist_bins_decade=kwargs.get(
                    "n_dist_bins_decade", PairCounter.N_DIST_BINS_DECADE_DEFAULT
                ),
            )
        else:
            stat = PairCounter.from_file(
                f,
                n_dist_bins_decade=kwargs.get(
                    "n_dist_bins_decade", PairCounter.N_DIST_BINS_DECADE_DEFAULT
                ),
            )
        stats.append(stat)
        f.close()

    # combine stats from several files (files_to_merge):
    out_stat = sum(stats)

    # Save merged stats.
    outstream = fileio.auto_open(
        output,
        mode="w",
        nproc=kwargs.get("nproc_out"),
        command=kwargs.get("cmd_out", None),
    )

    # save statistics to file ...
    out_stat.save(outstream)

    if outstream != sys.stdout:
        outstream.close()


def estimate_library_complexity(nseq, ndup, nopticaldup=0):
    """Estimate library complexity accounting for optical/clustering duplicates
    Parameters
    ----------
    nseq : int
        Total number of sequences
    ndup : int
        Total number of duplicates
    nopticaldup : int, optional
        Number of non-PCR duplicates, by default 0
    Returns
    -------
    float
        Estimated complexity
    """
    nseq = nseq - nopticaldup
    if nseq == 0:
        logger.warning("Empty of fully duplicated library, can't estimate complexity")
        return 0
    ndup = ndup - nopticaldup
    u = (nseq - ndup) / nseq
    if u == 0:
        logger.warning(
            "All the sequences are duplicates. Do you run complexity estimation on duplicates file?"
        )
        return 0
    seq_to_complexity = special.lambertw(-np.exp(-1 / u) / u).real + 1 / u
    complexity = float(nseq / seq_to_complexity)  # clean np.int64 data type
    return complexity


def analyse_bytile_duplicate_stats(df_dups, tile_dup_regex=False):
    """Count by-tile duplicates
    Parameters
    ----------
    dups : pd.DataFrame
        Dataframe with duplicates that contains pared read IDs
    tile_dup_regex : bool, optional
        See extract_tile_info for details, by default False
    Returns
    -------
    pd.DataFrame
        Grouped multi-indexed dataframe of pairwise by-tile duplication counts
    """

    df_dups = df_dups.copy()

    df_dups["tile"] = extract_tile_info(df_dups["readID"], regex=tile_dup_regex)
    df_dups["parent_tile"] = extract_tile_info(
        df_dups["parent_readID"], regex=tile_dup_regex
    )

    df_dups["same_tile"] = df_dups["tile"] == df_dups["parent_tile"]
    bytile_dups = (
        df_dups.groupby(["tile", "parent_tile"])
        .size()
        .reset_index(name="dup_count")
        .sort_values(["tile", "parent_tile"])
    )
    bytile_dups[["tile", "parent_tile"]] = np.sort(
        bytile_dups[["tile", "parent_tile"]].values, axis=1
    )
    bytile_dups = bytile_dups.groupby(["tile", "parent_tile"]).sum()
    return bytile_dups


def extract_tile_info(series, regex=False):
    """Extract the name of the tile for each read name in the series
    Parameters
    ----------
    series : pd.Series
        Series containing read IDs
    regex : bool, optional
        Regex to extract fields from the read IDs that correspond to tile IDs.
        By default False, uses a faster predefined approach for typical Illumina
        read names
        Example: r"(?:\w+):(?:\w+):(\w+):(\w+):(\w+):(?:\w+):(?:\w+)"
    Returns
    -------
    Series
        Series containing tile IDs as strings
    """
    if regex:
        split = series.str.extractall(regex).unstack().droplevel(1, axis=1)
        if split.shape[1] < 4:
            raise ValueError(
                f"Unable to convert tile names, does your readID have the tile information?\nHint: SRA removes tile information from readID.\nSample of your readIDs:\n{series.head()}"
            )
        return split[0] + ":" + split[1] + ":" + split[2]
    else:
        try:
            split = [":".join(name.split(":")[2:5]) for name in series]
        except:
            raise ValueError(
                f"Unable to convert tile names, does your readID have the tile information?\nHint: SRA removes tile information from readID.\nSample of your readIDs:\n{series.head()}"
            )
        return split


def yaml2pandas(yaml_path):
    """Generate a pandas DataFrame with stats from a yaml file

    Formats the keys within each filter using the PairCounter.flatten() method, to
    achieve same naming as in non-yaml stats files.

    Parameters
    ----------
    yaml_path : str
        Path to a yaml-formatted file with stats

    Returns
    -------
    pd.DataFrame
        Dataframe with filter names in the index and stats in columns
    """
    counter = PairCounter.from_yaml(open(yaml_path, "r"))
    stats = pd.concat(
        [
            pd.DataFrame(counter.flatten(filter=filter), index=[filter])
            for filter in counter.filters
        ]
    )
    return stats
