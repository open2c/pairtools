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

    def __init__(
        self,
        min_log10_dist=0,
        max_log10_dist=9,
        log10_dist_bin_step=0.25,
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
        self._dist_bins = np.r_[
            0,
            np.round(
                10
                ** np.arange(
                    min_log10_dist, max_log10_dist + 0.001, log10_dist_bin_step
                )
            ).astype(np.int),
        ]

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
            self._stat[key]["dedup"] = {}

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
                dist_bin_left = (
                    bin_range.strip("+")
                    if bin_range.endswith("+")
                    else bin_range.split("-")[0]
                )
                # get the index of that bin:
                bin_idx = (
                    np.searchsorted(self._dist_bins, int(dist_bin_left), "right") - 1
                )
                # store corresponding value:
                return self._stat[filter]["dist_freq"][dirs][bin_idx]
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

    def calculate_summaries(self):
        """calculate summary statistics (fraction of cis pairs at different cutoffs,
        complexity estimate) based on accumulated counts. Results are saved into
        self._stat["filter_name"]['summary"]
        """
        for key in self.filters.keys():
            self._stat[key]["summary"]["frac_dups"] = (
                (self._stat[key]["total_dups"] / self._stat[key]["total_mapped"])
                if self._stat[key]["total_mapped"] > 0
                else 0
            )

            for cis_count in (
                "cis",
                "cis_1kb+",
                "cis_2kb+",
                "cis_4kb+",
                "cis_10kb+",
                "cis_20kb+",
                "cis_40kb+",
            ):
                self._stat[key]["summary"][f"frac_{cis_count}"] = (
                    (self._stat[key][cis_count] / self._stat[key]["total_nodups"])
                    if self._stat[key]["total_nodups"] > 0
                    else 0
                )

            self._stat[key]["summary"][
                "complexity_naive"
            ] = estimate_library_complexity(
                self._stat[key]["total_mapped"], self._stat[key]["total_dups"], 0
            )
            if key == "no_filter" and self._save_bytile_dups:
                # Estimate library complexity with information by tile, if provided:
                if self._bytile_dups.shape[0] > 0:
                    self._stat[key]["dups_by_tile_median"] = int(
                        round(
                            self._bytile_dups["dup_count"].median()
                            * self._bytile_dups.shape[0]
                        )
                    )
                if "dups_by_tile_median" in self._stat[key].keys():
                    self._stat[key]["summary"][
                        "complexity_dups_by_tile_median"
                    ] = estimate_library_complexity(
                        self._stat[key]["total_mapped"],
                        self._stat[key]["total_dups"],
                        self._stat[key]["total_dups"]
                        - self._stat[key]["dups_by_tile_median"],
                    )

            self._summaries_calculated = True

    @classmethod
    def from_file(cls, file_handle):
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
        stat_from_file = cls()
        for l in file_handle:
            fields = l.strip().split(cls._SEP)
            if len(fields) == 0:
                # skip empty lines:
                continue
            if len(fields) != 2:
                # expect two _SEP separated values per line:
                raise fileio.ParseError(
                    "{} is not a valid stats file".format(file_handle.name)
                )
            # extract key and value, then split the key:
            putative_key, putative_val = fields[0], fields[1]
            key_fields = putative_key.split(cls._KEY_SEP)
            # we should impose a rigid structure of .stats or redo it:
            if len(key_fields) == 1:
                key = key_fields[0]
                if key in stat_from_file._stat[default_filter]:
                    stat_from_file._stat[default_filter][key] = int(fields[1])
                else:
                    raise fileio.ParseError(
                        "{} is not a valid stats file: unknown field {} detected".format(
                            file_handle.name, key
                        )
                    )
            else:
                # in this case key must be in ['pair_types','chrom_freq','dist_freq','dedup', 'summary']
                # get the first 'key' and keep the remainders in 'key_fields'
                key = key_fields.pop(0)
                if key in ["pair_types", "dedup", "summary", "chromsizes"]:
                    # assert there is only one element in key_fields left:
                    # 'pair_types', 'dedup', 'summary' and 'chromsizes' treated the same
                    if len(key_fields) == 1:
                        try:
                            stat_from_file._stat[default_filter][key][
                                key_fields[0]
                            ] = int(fields[1])
                        except ValueError:
                            stat_from_file._stat[default_filter][key][
                                key_fields[0]
                            ] = float(fields[1])
                    else:
                        raise fileio.ParseError(
                            "{} is not a valid stats file: {} section implies 1 identifier".format(
                                file_handle.name, key
                            )
                        )

                elif key == "chrom_freq":
                    # assert remaining key_fields == [chr1, chr2]:
                    if len(key_fields) == 2:
                        stat_from_file._stat[default_filter][key][
                            tuple(key_fields)
                        ] = int(fields[1])
                    else:
                        raise fileio.ParseError(
                            "{} is not a valid stats file: {} section implies 2 identifiers".format(
                                file_handle.name, key
                            )
                        )

                elif key == "dist_freq":
                    # assert that last element of key_fields is the 'directions'
                    if len(key_fields) == 2:
                        # assert 'dirs' in ['++','--','+-','-+']
                        dirs = key_fields.pop()
                        # there is only genomic distance range of the bin that's left:
                        (bin_range,) = key_fields
                        # extract left border of the bin "1000000+" or "1500-6000":
                        dist_bin_left = (
                            bin_range.strip("+")
                            if bin_range.endswith("+")
                            else bin_range.split("-")[0]
                        )
                        # get the index of that bin:
                        bin_idx = (
                            np.searchsorted(
                                stat_from_file._dist_bins, int(dist_bin_left), "right"
                            )
                            - 1
                        )
                        # store corresponding value:
                        stat_from_file._stat[default_filter][key][dirs][bin_idx] = int(
                            fields[1]
                        )
                    else:
                        raise fileio.ParseError(
                            "{} is not a valid stats file: {} section implies 2 identifiers".format(
                                file_handle.name, key
                            )
                        )
                else:
                    raise fileio.ParseError(
                        "{} is not a valid stats file: unknown field {} detected".format(
                            file_handle.name, key
                        )
                    )
        # return PairCounter from a non-empty dict:
        return stat_from_file

    @classmethod
    def from_yaml(cls, file_handle):
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
        stat_from_file = cls()

        stat = yaml.safe_load(file_handle)
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
        """

        self._stat[filter]["total"] += 1
        # collect pair type stats including DD:
        self._stat[filter]["pair_types"][pair_type] = (
            self._stat[filter]["pair_types"].get(pair_type, 0) + 1
        )
        if chrom1 == "!" and chrom2 == "!":
            self._stat[filter]["total_unmapped"] += 1
        elif chrom1 != "!" and chrom2 != "!":
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
                    bin = self._dist_bins[
                        np.searchsorted(self._dist_bins, dist, "right") - 1
                    ]
                    self._stat[filter]["dist_freq"][strand1 + strand2][bin] += 1
                    if dist >= 1000:
                        self._stat[filter]["cis_1kb+"] += 1
                    if dist >= 2000:
                        self._stat[filter]["cis_2kb+"] += 1
                    if dist >= 4000:
                        self._stat[filter]["cis_4kb+"] += 1
                    if dist >= 10000:
                        self._stat[filter]["cis_10kb+"] += 1
                    if dist >= 20000:
                        self._stat[filter]["cis_20kb+"] += 1
                    if dist >= 40000:
                        self._stat[filter]["cis_40kb+"] += 1

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
                k not in self._stat[filter] or k not in other._stat[filter]
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
                if k in ["pair_types", "dedup", "summary"]:
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
            # store nested dicts/arrays in a context dependet manner:
            # nested categories are stored only if they are non-trivial
            else:
                if (k == "dist_freq") and v:
                    for i in range(len(self._dist_bins)):
                        for dirs, freqs in v.items():
                            # last bin is treated differently: "100000+" vs "1200-3000":
                            if i != len(self._dist_bins) - 1:
                                dist = self._dist_bins[i]
                                dist_next = self._dist_bins[i + 1]
                                formatted_key = self._KEY_SEP.join(
                                    ["{}", "{}-{}", "{}"]
                                ).format(k, dist, dist_next, dirs)
                            else:
                                formatted_key = self._KEY_SEP.join(
                                    ["{}", "{}+", "{}"]
                                ).format(k, dist, dirs)
                            # store key,value pair:
                            flat_stat[formatted_key] = freqs[dist]
                elif (k in ["pair_types", "dedup", "chromsizes"]) and v:
                    # 'pair_types' and 'dedup' are simple dicts inside,
                    # treat them the exact same way:
                    for k_item, freq in v.items():
                        formatted_key = self._KEY_SEP.join(["{}", "{}"]).format(
                            k, k_item
                        )
                        # store key,value pair:
                        flat_stat[formatted_key] = freq
                elif (k == "chrom_freq") and v:
                    for (chrom1, chrom2), freq in v.items():
                        formatted_key = self._KEY_SEP.join(["{}", "{}", "{}"]).format(
                            k, chrom1, chrom2
                        )
                        # store key,value pair:
                        flat_stat[formatted_key] = freq
                elif (k == "summary") and v:
                    for key, frac in v.items():
                        formatted_key = self._KEY_SEP.join(["{}", "{}"]).format(k, key)
                        # store key,value pair:
                        flat_stat[formatted_key] = frac

        # return flattened dict
        return flat_stat

    def format_yaml(self, filter="no_filter"):
        """return a formatted dict (for the yaml output)
        Performed for all filters at once."""

        from copy import deepcopy

        formatted_stat = {key: {} for key in self.filters.keys()}

        # Storing statistics for each filter
        for key in self.filters.keys():
            for k, v in self._stat[key].items():
                if isinstance(v, int):
                    formatted_stat[key][k] = v
                # store nested dicts/arrays in a context dependet manner:
                # nested categories are stored only if they are non-trivial
                else:
                    if (k != "chrom_freq") and v:
                        # simple dicts inside
                        # treat them the exact same way:
                        formatted_stat[key][k] = deepcopy(v)
                    elif (k == "chrom_freq") and v:
                        # need to convert tuples of chromosome names to str
                        freqs = {}
                        for (chrom1, chrom2), freq in sorted(v.items()):
                            freqs[
                                self._KEY_SEP.join(["{}", "{}"]).format(chrom1, chrom2)
                            ] = freq
                            # store key,value pair:
                        formatted_stat[key][k] = deepcopy(freqs)
            # return formatted dict
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
            stat = PairCounter.from_yaml(f)
        else:
            stat = PairCounter.from_file(f)
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
