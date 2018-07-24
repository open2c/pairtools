#!/usr/bin/env python
# -*- coding: utf-8  -*-
import sys
import ast 
import warnings
import pathlib

import click

import numpy as np

from . import _dedup, _fileio, _pairsam_format, _headerops, cli, common_io_options
from .pairtools_markasdup import mark_split_pair_as_dup
from .pairtools_stats import PairCounter

UTIL_NAME = 'pairtools_filterbycov'

######################################
## TODO: - output stats after filtering
## edit/update mark as dup to mark as multi
###################################

@cli.command()
@click.argument(
    'pairs_path', 
    type=str,
    required=False)
@click.option(
    "-o", "--output", 
    type=str, 
    default="", 
    help='output file for pairs from low coverage regions.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed.'
        ' By default, the output is printed into stdout.')
@click.option(
    "--output-highcov",
    type=str, 
    default="", 
    help='output file for pairs from high coverage regions.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed.'
        ' If the path is the same as in --output or -, output duplicates together '
        ' with deduped pairs. By default, duplicates are dropped.')
@click.option(
    "--output-unmapped",
    type=str, 
    default="", 
    help='output file for unmapped pairs. '
        'If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed. '
        'If the path is the same as in --output or -, output unmapped pairs together '
        'with deduped pairs. If the path is the same as --output-highcov, '
        'output unmapped reads together. By default, unmapped pairs are dropped.')
@click.option(
    "--output-stats", 
    type=str, 
    default="", 
    help='output file for statistics of multiple interactors. '
        ' If file exists, it will be open in the append mode.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed.'
        ' By default, statistics are not printed.')
@click.option(
    "--max-cov",
    type=int, 
    default=8,
    help='The maximum allowed coverage per region.'
    )
@click.option(
    "--max-dist",
    type=int, 
    default=500,
    help='The resolution for calculating coverage. For each pair, the local '
    'coverage around each end is calculated as (1 + the number of neighbouring '
    'pairs within +/- max_dist bp) ')
@click.option(
    '--method',
    type=click.Choice(['max', 'sum']),
    default="max",  
    help='calculate the number of neighbouring pairs as either the sum or the max'
    ' of the number of neighbours on the two sides',
    show_default=True)
@click.option(
    "--sep",
    type=str, 
    default=_pairsam_format.PAIRSAM_SEP_ESCAPE, 
    help=r"Separator (\t, \v, etc. characters are "
          "supported, pass them in quotes) ")
@click.option(
    "--comment-char", 
    type=str, 
    default="#", 
    help="The first character of comment lines")
@click.option(
    "--send-header-to", 
    type=click.Choice(['lowcov', 'highcov', 'both', 'none']),
    default="both", 
    help="Which of the outputs should receive header and comment lines")
@click.option(
    "--c1", 
    type=int, 
    default=_pairsam_format.COL_C1,  
    help='Chrom 1 column; default {}'.format(_pairsam_format.COL_C1))
@click.option(
    "--c2", 
    type=int, 
    default=_pairsam_format.COL_C2,  
    help='Chrom 2 column; default {}'.format(_pairsam_format.COL_C2))
@click.option(
    "--p1", 
    type=int, 
    default=_pairsam_format.COL_P1,  
    help='Position 1 column; default {}'.format(_pairsam_format.COL_P1))
@click.option(
    "--p2", 
    type=int, 
    default=_pairsam_format.COL_P2,  
    help='Position 2 column; default {}'.format(_pairsam_format.COL_P2))
@click.option(
    "--s1", 
    type=int, 
    default=_pairsam_format.COL_S1,  
    help='Strand 1 column; default {}'.format(_pairsam_format.COL_S1))
@click.option(
    "--s2", 
    type=int, 
    default=_pairsam_format.COL_S2,  
    help='Strand 2 column; default {}'.format(_pairsam_format.COL_S2))
@click.option(
    "--unmapped-chrom", 
    type=str, 
    default=_pairsam_format.UNMAPPED_CHROM,  
    help='Placeholder for a chromosome on an unmapped side; default {}'.format(_pairsam_format.UNMAPPED_CHROM))
@click.option(
    "--mark-multi", 
    is_flag=True,
    help='If specified, duplicate pairs are marked as FF in "pair_type" and '
         'as a duplicate in the sam entries.')

@common_io_options

def filterbycov(
    pairs_path, output, output_highcov,
    output_unmapped, output_stats,
    max_dist,max_cov, method, 
    sep, comment_char, send_header_to,
    c1, c2, p1, p2, s1, s2, unmapped_chrom, mark_multi, **kwargs
    ):
    '''Remove pairs from regions of high coverage. 
    
    Find and remove pairs with >(MAX_COV-1) neighbouring pairs
    within a +/- MAX_DIST bp window around either side. Useful for single-cell 
    Hi-C experiments, where coverage is naturally limited by the chromosome 
    copy number.

    PAIRS_PATH : input triu-flipped sorted .pairs or .pairsam file.  If the
    path ends with .gz/.lz4, the input is decompressed by pbgzip/lz4c. 
    By default, the input is read from stdin.
    '''
    filterbycov_py(
        pairs_path, output, output_highcov,
        output_unmapped,output_stats,
        max_dist,max_cov, method, 
        sep, comment_char, send_header_to,
        c1, c2, p1, p2, s1, s2, unmapped_chrom, mark_multi,
        **kwargs
        )
    
    
def filterbycov_py(
    pairs_path, output, output_highcov,
    output_unmapped, output_stats,
    max_dist,max_cov, method, 
    sep, comment_char, send_header_to,
    c1, c2, p1, p2, s1, s2, unmapped_chrom, mark_multi,
    **kwargs
    ):
    
    
    ## Prepare input, output streams based on selected outputs
    ## Default ouput stream is low-frequency interactors
    sep = ast.literal_eval('"""' + sep + '"""')
    send_header_to_lowcov = send_header_to in ['both', 'lowcov']
    send_header_to_highcov = send_header_to in ['both', 'highcov']

    instream = (_fileio.auto_open(pairs_path, mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if pairs_path else sys.stdin)
    outstream = (_fileio.auto_open(output, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)
    out_stats_stream = (_fileio.auto_open(output_stats, mode='w', 
                           nproc=kwargs.get('nproc_out'),
                           command=kwargs.get('cmd_out', None)) 
             if output_stats else None)

    # generate empty PairCounter if stats output is requested:
    out_stat = PairCounter() if output_stats else None    
    
    # output the high-frequency interacting pairs
    if not output_highcov:
        outstream_high = None
    elif (output_highcov == '-' or 
          (pathlib.Path(output_highcov).absolute() == pathlib.Path(output).absolute())):
        outstream_high = outstream
    else:
        outstream_high = _fileio.auto_open(output_highcov, mode='w', 
                                            nproc=kwargs.get('nproc_out'),
                                            command=kwargs.get('cmd_out', None)) 
    
    # output unmapped pairs
    if not output_unmapped:
        outstream_unmapped = None
    elif (output_unmapped == '-' or 
        (pathlib.Path(output_unmapped).absolute() == pathlib.Path(output).absolute())):
        outstream_unmapped = outstream
    elif (pathlib.Path(output_unmapped).absolute() == pathlib.Path(output_highcov).absolute()):
        outstream_unmapped = outstream_high
    else:
        outstream_unmapped = _fileio.auto_open(output_unmapped, mode='w', 
                                            nproc=kwargs.get('nproc_out'),
                                            command=kwargs.get('cmd_out', None)) 
        
    # prepare file headers
    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)

    # header for low-frequency interactors
    if send_header_to_lowcov:
        outstream.writelines((l+'\n' for l in header))
        
    # header for high-frequency interactors
    if send_header_to_highcov and outstream_high and (outstream_high != outstream):
        outstream_high.writelines((l+'\n' for l in header))
        
    # header for unmapped pairs
    if (outstream_unmapped and (outstream_unmapped != outstream) 
            and (outstream_unmapped != outstream_high)):
        outstream_unmapped.writelines((l+'\n' for l in header))
  
    # perform filtering of pairs based on low/high-frequency of interaction
    streaming_filterbycov( method, max_dist,max_cov, sep,
        c1, c2, p1, p2, s1, s2, unmapped_chrom,
        body_stream, outstream, outstream_high,
        outstream_unmapped, out_stat, mark_multi)
    
    ## FINISHED!
    # save statistics to a file if it was requested: TO BE TESTED
    if out_stat:
        out_stat.save(out_stats_stream)

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()

    if outstream_high and (outstream_high != outstream):
        outstream_high.close()

    if (outstream_unmapped and (outstream_unmapped != outstream) 
            and (outstream_unmapped != outstream_high)):
        outstream_unmapped.close()

    if out_stats_stream:
        out_stats_stream.close()
    
    
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

    c1 = np.asarray(c1_in,dtype=int)
    p1 = np.asarray(p1_in,dtype=int)
    c2 = np.asarray(c2_in,dtype=int)
    p2 = np.asarray(p2_in,dtype=int)
    
    M = np.r_[np.c_[c1,p1],np.c_[c2,p2]] # M is a table of (chrom, pos) with 2*N rows

    assert (c1.shape[0] == c2.shape[0])
    N = 2*c1.shape[0]

    ind_sorted = np.lexsort((M[:,1],M[:,0])) # sort by chromosomes, then positions
    # M[ind_sorted]
    # ind_sorted
    # M, M[ind_sorted]


    if (method == 'sum'):
        proximity_count = np.zeros(N) # keeps track of how many molecules each framgent end is close to
    elif (method == 'max'):
        proximity_count = np.zeros(N) 
    else:
        raise ValueError('Unknown method: {}'.format(method))

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
        if M[ind_sorted[low],0] != M[ind_sorted[high],0]: 
            low += 1
            high = low + 1  # restart high 
            continue

        # next, if positions are not proximal, increase low, and continue
        elif np.abs(M[ind_sorted[high],1] - M[ind_sorted[low],1]) > max_dist:
            low += 1
            high = low + 1  # restart high 
            continue

        # if on the same chromosome, and the distance is "proximal enough", add to count of both "low" and "high" positions
        else:
            proximity_count[low] += 1
            proximity_count[high] += 1

        high += 1

    # unsort proximity count
    #proximity_count = proximity_count[ind_sorted]
    proximity_count[ind_sorted] = np.copy(proximity_count)
    #print(M)
    #print(proximity_count)

    # if method is sum of pairs
    if method == 'sum':
        pcounts = proximity_count[0:N//2] + proximity_count[N//2:] + 1
    elif method == 'max':
        pcounts = np.maximum(proximity_count[0:N//2]+1,
                             proximity_count[N//2:]+1)
    else:
        raise ValueError('Unknown method: {}'.format(method))
       
    return pcounts


def streaming_filterbycov( method, max_dist, max_cov, sep,
        c1ind, c2ind, p1ind, p2ind, s1ind, s2ind, unmapped_chrom,
        instream, outstream, outstream_high,
        outstream_unmapped, out_stat, mark_multi):
    
    # doing everything in memory
    maxind = max(c1ind, c2ind, p1ind, p2ind, s1ind, s2ind)

    # if we do stats in the dedup, we need PAIR_TYPE
    # i do not see way around this:
    if out_stat:
        ptind = _pairsam_format.COL_PTYPE
        maxind = max(maxind, ptind)

    c1 = []; c2 = []; p1 = []; p2 = []; s1 = []; s2 = []
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
                    + " expected {} columns, got {}".format(maxind, len(cols)))
                
            if ((cols[c1ind] == unmapped_chrom)
                or (cols[c2ind] == unmapped_chrom)):

                if outstream_unmapped:
                    outstream_unmapped.write(stripline)
                    # don't forget terminal newline
                    outstream_unmapped.write("\n")

                # add a pair to PairCounter if stats output is requested:
                if out_stat:
                    out_stat.add_pair(cols[c1ind],  int(cols[p1ind]),  cols[s1ind],
                                      cols[c2ind],  int(cols[p2ind]),  cols[s2ind],
                                      cols[ptind])                    
            else:
                line_buffer.append(stripline)
                cols_buffer.append(cols)

                c1.append(fetchadd(cols[c1ind], chromDict))
                c2.append(fetchadd(cols[c2ind], chromDict))
                p1.append(int(cols[p1ind]))
                p2.append(int(cols[p2ind]))
                s1.append(fetchadd(cols[s1ind], strandDict))
                s2.append(fetchadd(cols[s2ind], strandDict))    
    
        else: # when everything is loaded in memory... 

            res = _filterbycov(c1, p1, c2, p2, max_dist, method)

            for i in range(len(res)):
                # not high-frequency interactor pairs:
                if not res[i] > max_cov:
                    outstream.write(line_buffer[i])
                    # don't forget terminal newline
                    outstream.write("\n")
                    if out_stat:
                        out_stat.add_pair(cols_buffer[i][c1ind],
                                          int(cols_buffer[i][p1ind]),
                                          cols_buffer[i][s1ind],
                                          cols_buffer[i][c2ind],
                                          int(cols_buffer[i][p2ind]),
                                          cols_buffer[i][s2ind],
                                          cols_buffer[i][ptind])
                # high-frequency interactor pairs:
                else:
                    if out_stat:
                        out_stat.add_pair(cols_buffer[i][c1ind],
                                          int(cols_buffer[i][p1ind]),
                                          cols_buffer[i][s1ind],
                                          cols_buffer[i][c2ind],
                                          int(cols_buffer[i][p2ind]),
                                          cols_buffer[i][s2ind],
                                          'FF' )
                    if outstream_high:
                        outstream_high.write(
                          # FF-marked pair:
                          sep.join(mark_split_pair_as_dup(cols_buffer[i])) 
                          if mark_multi
                          # pair as is:
                          else line_buffer[i] )
                        # don't forget terminal newline
                        outstream_high.write('\n')
                    
            # flush buffers and perform necessary checks here:
            c1 = []; c2 = []; p1 = []; p2 = []; s1 = []; s2 = []
            line_buffer = line_buffer[len(res):]
            cols_buffer = cols_buffer[len(res):]
            if not stripline:
                if(len(line_buffer) != 0):                
                    raise ValueError(
                        "{} lines left in the buffer, ".format(len(line_buffer))
                        + "should be none;"
                        + "something went terribly wrong")
                break
                      
            break

    
    
if __name__ == '__main__':
    filterbycov()
