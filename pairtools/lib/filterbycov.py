import numpy as np
import warnings

from .dedup import mark_split_pair_as_dup
from . import pairsam_format


def fetchadd(key, mydict):
    key = key.strip()
    if key not in mydict:
        mydict[key] = len(mydict)
    return mydict[key]


def ar(mylist, val):
    return np.array(mylist, dtype={8: np.int8, 16: np.int16, 32: np.int32}[val])


def _filterbycov(c1_in, p1_in, c2_in, p2_in, max_dist, method):
    """
    This is a slow version of the filtering code used for testing purposes only
    Use cythonized version in the future!!
    """

    c1 = np.asarray(c1_in, dtype=int)
    p1 = np.asarray(p1_in, dtype=int)
    c2 = np.asarray(c2_in, dtype=int)
    p2 = np.asarray(p2_in, dtype=int)

    M = np.r_[
        np.c_[c1, p1], np.c_[c2, p2]
    ]  # M is a table of (chrom, pos) with 2*N rows

    assert c1.shape[0] == c2.shape[0]
    N = 2 * c1.shape[0]

    ind_sorted = np.lexsort((M[:, 1], M[:, 0]))  # sort by chromosomes, then positions
    # M[ind_sorted]
    # ind_sorted
    # M, M[ind_sorted]

    if method == "sum":
        proximity_count = np.zeros(
            N
        )  # keeps track of how many molecules each framgent end is close to
    elif method == "max":
        proximity_count = np.zeros(N)
    else:
        raise ValueError("Unknown method: {}".format(method))

    low = 0
    high = 1
    while True:

        # boundary case finish
        if low == N:
            break

        # boundary case  - CHECK
        if high == N:
            low += 1
            high = low + 1
            continue

        # check if "high" is proximal enough to "low"

        # first, if chromosomes not equal, we have gone too far, and the positions are not proximal
        if M[ind_sorted[low], 0] != M[ind_sorted[high], 0]:
            low += 1
            high = low + 1  # restart high
            continue

        # next, if positions are not proximal, increase low, and continue
        elif np.abs(M[ind_sorted[high], 1] - M[ind_sorted[low], 1]) > max_dist:
            low += 1
            high = low + 1  # restart high
            continue

        # if on the same chromosome, and the distance is "proximal enough", add to count of both "low" and "high" positions
        else:
            proximity_count[low] += 1
            proximity_count[high] += 1

        high += 1

    # unsort proximity count
    # proximity_count = proximity_count[ind_sorted]
    proximity_count[ind_sorted] = np.copy(proximity_count)
    # print(M)
    # print(proximity_count)

    # if method is sum of pairs
    if method == "sum":
        pcounts = proximity_count[0 : N // 2] + proximity_count[N // 2 :] + 1
    elif method == "max":
        pcounts = np.maximum(
            proximity_count[0 : N // 2] + 1, proximity_count[N // 2 :] + 1
        )
    else:
        raise ValueError("Unknown method: {}".format(method))

    return pcounts


def streaming_filterbycov(
    method,
    max_dist,
    max_cov,
    sep,
    c1ind,
    c2ind,
    p1ind,
    p2ind,
    s1ind,
    s2ind,
    unmapped_chrom,
    instream,
    outstream,
    outstream_high,
    outstream_unmapped,
    out_stat,
    mark_multi,
):

    # doing everything in memory
    maxind = max(c1ind, c2ind, p1ind, p2ind, s1ind, s2ind)

    # if we do stats in the dedup, we need PAIR_TYPE
    # i do not see way around this:
    if out_stat:
        ptind = pairsam_format.COL_PTYPE
        maxind = max(maxind, ptind)

    c1 = []
    c2 = []
    p1 = []
    p2 = []
    s1 = []
    s2 = []
    line_buffer = []
    cols_buffer = []
    chromDict = {}
    strandDict = {}
    n_unmapped = 0
    n_high = 0
    n_low = 0

    instream = iter(instream)
    while True:
        rawline = next(instream, None)
        stripline = rawline.strip() if rawline else None

        # take care of empty lines not at the end of the file separately
        if rawline and (not stripline):
            warnings.warn("Empty line detected not at the end of the file")
            continue

        if stripline:
            cols = stripline.split(sep)
            if len(cols) <= maxind:
                raise ValueError(
                    "Error parsing line {}: ".format(stripline)
                    + " expected {} columns, got {}".format(maxind, len(cols))
                )

            if (cols[c1ind] == unmapped_chrom) or (cols[c2ind] == unmapped_chrom):

                if outstream_unmapped:
                    outstream_unmapped.write(stripline)
                    # don't forget terminal newline
                    outstream_unmapped.write("\n")

                # add a pair to PairCounter if stats output is requested:
                if out_stat:
                    out_stat.add_pair(
                        cols[c1ind],
                        int(cols[p1ind]),
                        cols[s1ind],
                        cols[c2ind],
                        int(cols[p2ind]),
                        cols[s2ind],
                        cols[ptind],
                    )
            else:
                line_buffer.append(stripline)
                cols_buffer.append(cols)

                c1.append(fetchadd(cols[c1ind], chromDict))
                c2.append(fetchadd(cols[c2ind], chromDict))
                p1.append(int(cols[p1ind]))
                p2.append(int(cols[p2ind]))
                s1.append(fetchadd(cols[s1ind], strandDict))
                s2.append(fetchadd(cols[s2ind], strandDict))

        else:  # when everything is loaded in memory...

            res = _filterbycov(c1, p1, c2, p2, max_dist, method)

            for i in range(len(res)):
                # not high-frequency interactor pairs:
                if not res[i] > max_cov:
                    outstream.write(line_buffer[i])
                    # don't forget terminal newline
                    outstream.write("\n")
                    if out_stat:
                        out_stat.add_pair(
                            cols_buffer[i][c1ind],
                            int(cols_buffer[i][p1ind]),
                            cols_buffer[i][s1ind],
                            cols_buffer[i][c2ind],
                            int(cols_buffer[i][p2ind]),
                            cols_buffer[i][s2ind],
                            cols_buffer[i][ptind],
                        )
                # high-frequency interactor pairs:
                else:
                    if out_stat:
                        out_stat.add_pair(
                            cols_buffer[i][c1ind],
                            int(cols_buffer[i][p1ind]),
                            cols_buffer[i][s1ind],
                            cols_buffer[i][c2ind],
                            int(cols_buffer[i][p2ind]),
                            cols_buffer[i][s2ind],
                            "FF",
                        )
                    if outstream_high:
                        outstream_high.write(
                            # DD-marked pair:
                            sep.join(mark_split_pair_as_dup(cols_buffer[i]))
                            if mark_multi
                            # pair as is:
                            else line_buffer[i]
                        )
                        # don't forget terminal newline
                        outstream_high.write("\n")

            # flush buffers and perform necessary checks here:
            c1 = []
            c2 = []
            p1 = []
            p2 = []
            s1 = []
            s2 = []
            line_buffer = line_buffer[len(res) :]
            cols_buffer = cols_buffer[len(res) :]
            if not stripline:
                if len(line_buffer) != 0:
                    raise ValueError(
                        "{} lines left in the buffer, ".format(len(line_buffer))
                        + "should be none;"
                        + "something went terribly wrong"
                    )
                break

            break
