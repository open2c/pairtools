#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import OrderedDict
import subprocess
import fileinput
import itertools
import click
import pipes
import sys
import os
import io

from . import _fileio, _pairsam_format, _headerops, cli, common_io_options
from .pairtools_stats import PairCounter


UTIL_NAME = 'pairtools_parse'

EXTRA_COLUMNS = [
    'mapq',
    'pos5',
    'pos3',
    'cigar',
    'read_len',
    'matched_bp',
    'algn_ref_span',
    'algn_read_span',
    'dist_to_5',
    'dist_to_3',
    'seq'
]

@cli.command()
@click.argument(
    'sam_path', 
    type=str,
    required=False)
@click.option(
    "-c", "--chroms-path",
    type=str,
    required=True,
    help='Chromosome order used to flip interchromosomal mates: '
         'path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose '
         'first column lists scaffold names. Any scaffolds not listed will be '
         'ordered lexicographically following the names provided.')
@click.option(
    "-o", "--output", 
    type=str, 
    default="", 
    help='output file. '
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4-compressed.'
         'By default, the output is printed into stdout. ')
@click.option(
    "--assembly", 
    type=str, 
    help='Name of genome assembly (e.g. hg19, mm10) to store in the pairs header.')
@click.option(
    "--min-mapq", 
    type=int, 
    default=1,
    show_default=True,
    help='The minimal MAPQ score to consider a read as uniquely mapped')
@click.option(
    "--max-molecule-size", 
    type=int, 
    default=2000,
    show_default=True,
    help='The maximal size of a Hi-C molecule; used to rescue single ligations'
         'from molecules with three alignments.')
@click.option(
    "--drop-readid", 
    is_flag=True,
    help='If specified, do not add read ids to the output')
@click.option(
    "--drop-seq", 
    is_flag=True,
    help='If specified, remove sequences and PHREDs from the sam fields')
@click.option(
    "--drop-sam", 
    is_flag=True,
    help='If specified, do not add sams to the output')
@click.option(
    "--add-columns", 
    type=click.STRING,
    default='',
    help='Report extra columns describing alignments '
         'Possible values (can take multiple values as a comma-separated '
         'list): a SAM tag (any pair of uppercase letters) or {}.'.format(
             ', '.join(EXTRA_COLUMNS)))
@click.option(
    "--output-parsed-alignments", 
    type=str, 
    default="", 
    help='output file for all parsed alignments, including walks.'
        ' Useful for debugging and rnalysis of walks.'
        ' If file exists, it will be open in the append mode.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4-compressed.'
        ' By default, not used.'
        )
@click.option(
    "--output-stats", 
    type=str, 
    default="", 
    help='output file for various statistics of pairs file. '
        ' By default, statistics is not generated.')
@click.option(
    '--report-alignment-end', 
    type=click.Choice(['5', '3']),
    default='5',
    help='specifies whether the 5\' or 3\' end of the alignment is reported as'
    ' the position of the Hi-C read.')
@click.option(
    '--max-inter-align-gap',
    type=int,
    default=20,
    show_default=True,
    help='read segments that are not covered by any alignment and'
         ' longer than the specified value are treated as "null" alignments.'
         ' These null alignments convert otherwise linear alignments into walks,' 
         ' and affect how they get reported as a Hi-C pair (see --walks-policy).' 
    )
@click.option(
    "--walks-policy", 
    type=click.Choice(['mask', 'all', '5any', '5unique', '3any', '3unique']),
    default='mask',
    help='the policy for reporting unrescuable walks (reads containing more'
    ' than one alignment on one or both sides, that can not be explained by a'
    ' single ligation between two mappable DNA fragments).'
    ' "mask" - mask walks (chrom="!", pos=0, strand="-"); '
    ' "all" - report all pairs of consecutive alignments [NOT IMPLEMENTED]; '
    ' "5any" - report the 5\'-most alignment on each side;'
    ' "5unique" - report the 5\'-most unique alignment on each side, if present;'
    ' "3any" - report the 3\'-most alignment on each side;'
    ' "3unique" - report the 3\'-most unique alignment on each side, if present.',
    show_default=True
    )
@click.option(
    "--no-flip", 
    is_flag=True,
    help='If specified, do not flip pairs in genomic order and instead preserve ' 
         'the order in which they were sequenced.')

@common_io_options

def parse(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
          drop_readid, drop_seq, drop_sam, add_columns,
          output_parsed_alignments, output_stats, **kwargs):
    '''Find ligation junctions in .sam, make .pairs.

    SAM_PATH : an input .sam/.bam file with paired-end sequence alignments of
    Hi-C molecules. If the path ends with .bam, the input is decompressed from
    bam with samtools. By default, the input is read from stdin.
    '''
    parse_py(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
             drop_readid, drop_seq, drop_sam, add_columns,
             output_parsed_alignments, output_stats, **kwargs)


def parse_py(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
             drop_readid, drop_seq, drop_sam, add_columns, 
             output_parsed_alignments, output_stats, **kwargs):
    instream = (_fileio.auto_open(sam_path, mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if sam_path else sys.stdin)
    outstream = (_fileio.auto_open(output, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)
    out_alignments_stream = (_fileio.auto_open(output_parsed_alignments, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output_parsed_alignments else None)
    out_stats_stream = (_fileio.auto_open(output_stats, mode='w', 
                               nproc=kwargs.get('nproc_out'),
                               command=kwargs.get('cmd_out', None)) 
                 if output_stats else None)

    if out_alignments_stream:
        out_alignments_stream.write('readID\tside\tchrom\tpos\tstrand\tmapq\tcigar\tdist_5_lo\tdist_5_hi\tmatched_bp\n')

    # generate empty PairCounter if stats output is requested:
    out_stat = PairCounter() if output_stats else None

    samheader, body_stream = _headerops.get_header(instream, comment_char='@')
    sam_chromsizes = _headerops.get_chromsizes_from_sam_header(samheader)
    chromosomes = _headerops.get_chrom_order(
        chroms_path, 
        list(sam_chromsizes.keys()))

    add_columns = [col for col in add_columns.split(',') if col]
    for col in add_columns:
        if not( (col in EXTRA_COLUMNS) or (len(col) == 2 and col.isupper())):
            raise Exception('{} is not a valid extra column'.format(col))

    columns =  (_pairsam_format.COLUMNS 
                + ([c+side for c in add_columns for side in ['1', '2']])
                )

    if drop_sam:
        columns.pop(columns.index('sam1'))
        columns.pop(columns.index('sam2'))

    header = _headerops.make_standard_pairsheader(
        assembly = assembly,
        chromsizes = [(chrom, sam_chromsizes[chrom]) for chrom in chromosomes],
        columns = columns, 
        shape = 'whole matrix' if kwargs['no_flip'] else 'upper triangle'

    )
    
    header = _headerops.insert_samheader(header, samheader) 
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l+'\n' for l in header))

    streaming_classify(body_stream, outstream, chromosomes, min_mapq, 
                       max_molecule_size, drop_readid, drop_seq, drop_sam, 
                       add_columns, out_alignments_stream, out_stat, **kwargs)

    # save statistics to a file if it was requested:
    if out_stat:
        out_stat.save(out_stats_stream)

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
    # close optional output streams if needed:
    if out_alignments_stream:
        out_alignments_stream.close()
    if out_stats_stream:
        out_stats_stream.close()


def parse_cigar(cigar):
    matched_bp = 0
    algn_ref_span = 0
    algn_read_span = 0
    read_len = 0
    clip5_ref = 0
    clip3_ref = 0

    if cigar != '*':
        cur_num = 0
        for char in cigar:
            charval = ord(char)
            if charval >= 48 and charval <= 57:
                cur_num = cur_num * 10 + (charval - 48)
            else:
                if char == 'M':
                    matched_bp += cur_num
                    algn_ref_span += cur_num
                    algn_read_span += cur_num
                    read_len += cur_num
                elif char == 'I':
                    algn_read_span += cur_num
                    read_len += cur_num
                elif char == 'D':
                    algn_ref_span += cur_num
                elif char == 'S' or char == 'H':
                    read_len += cur_num
                    if matched_bp == 0:
                        clip5_ref = cur_num
                    else:
                        clip3_ref = cur_num

                cur_num = 0

    return {
        'clip5_ref': clip5_ref,
        'clip3_ref': clip3_ref,
        'cigar': cigar,
        'algn_ref_span': algn_ref_span,
        'algn_read_span': algn_read_span,
        'read_len': read_len,
        'matched_bp': matched_bp,
    }


def empty_alignment():
    return {
        'chrom': _pairsam_format.UNMAPPED_CHROM,
        'pos5': _pairsam_format.UNMAPPED_POS,
        'pos3': _pairsam_format.UNMAPPED_POS, 
        'pos': _pairsam_format.UNMAPPED_POS, 
        'strand': _pairsam_format.UNMAPPED_STRAND, 
        'dist_to_5': 0, 
        'dist_to_3': 0, 
        'mapq': 0, 
        'is_unique': False, 
        'is_mapped': False, 
        'is_linear': True, 
        'cigar' : '*',
        'algn_ref_span': 0, 
        'algn_read_span': 0,
        'matched_bp': 0, 
        'clip3_ref': 0,
        'clip5_ref': 0, 
        'read_len': 0,
        'type':'N'
    }


def parse_algn(
        samcols, 
        min_mapq, 
        report_3_alignment_end=False, 
        sam_tags=None,
        store_seq=False):
    is_mapped = (int(samcols[1]) & 0x04) == 0
    mapq = int(samcols[4])
    is_unique = (mapq >= min_mapq)
    is_linear = not any([col.startswith('SA:Z:') for col in samcols[11:]])

    cigar = parse_cigar(samcols[5])

    if is_mapped: 
        if ((int(samcols[1]) & 0x10) == 0):
            strand = '+'
            dist_to_5 = cigar['clip5_ref']
            dist_to_3 = cigar['clip3_ref']
        else:
            strand = '-'
            dist_to_5 = cigar['clip3_ref']
            dist_to_3 = cigar['clip5_ref']

        if is_unique:
            chrom = samcols[2] 
            if strand == '+':
                pos5 = int(samcols[3])
                pos3 = int(samcols[3]) + cigar['algn_ref_span'] - 1
            else:
                pos5 = int(samcols[3]) + cigar['algn_ref_span'] - 1
                pos3 = int(samcols[3])
        else:
            chrom = _pairsam_format.UNMAPPED_CHROM
            strand = _pairsam_format.UNMAPPED_STRAND
            pos5 = _pairsam_format.UNMAPPED_POS
            pos3 = _pairsam_format.UNMAPPED_POS
    else:
        chrom = _pairsam_format.UNMAPPED_CHROM
        strand = _pairsam_format.UNMAPPED_STRAND
        pos5 = _pairsam_format.UNMAPPED_POS
        pos3 = _pairsam_format.UNMAPPED_POS

        dist_to_5 = 0
        dist_to_3 = 0

    algn = {
        'chrom': chrom,
        'pos5': pos5,
        'pos3': pos3,
        'strand': strand,
        'mapq': mapq,
        'is_mapped': is_mapped,
        'is_unique': is_unique,
        'is_linear': is_linear,
        'dist_to_5': dist_to_5,
        'dist_to_3': dist_to_3,
        'type': ('N' if not is_mapped else ('M' if not is_unique else 'U'))
    }

    algn.update(cigar)

    algn['pos'] = algn['pos3'] if report_3_alignment_end else algn['pos5']

    if sam_tags:
        for tag in sam_tags:
            algn[tag] = ''
            
        for col in samcols[11:]:
            for tag in sam_tags:
                if col.startswith(tag+':'):
                    algn[tag] = col[5:]
                    continue

    if store_seq:
        algn['seq'] = samcols[9]  

    return algn


def rescue_walk(algns1, algns2, max_molecule_size):
    """
    Rescue a single ligation that appears as a walk.

    Checks if a molecule with three alignments could be formed via a single 
    ligation between two fragments, where one fragment was so long that it 
    got sequenced on both sides.

    Uses three criteria:
    a) the 3'-end alignment on one side maps to the same chromosome as the 
    alignment fully covering the other side (i.e. the linear alignment)
    b) the two alignments point towards each other on the chromosome
    c) the distance between the outer ends of the two alignments is below
    the specified threshold.


    Alternatively, a single ligation get rescued when the 3' sub-alignment 
    maps to multiple locations or no locations at all.
    
    In the case of a successful rescue, tags the 3' sub-alignment with
    type='X' and the linear alignment on the other side with type='R'.

    Returns
    -------
    linear_side : int
        If the case of a successful rescue, returns the index of the side
        with a linear alignment.

    """
    
    # If both sides have one alignment or none, no need to rescue!
    n_algns1 = len(algns1)
    n_algns2 = len(algns2)

    if (n_algns1 <= 1) and (n_algns2 <= 1):
        return None

    # Can rescue only pairs with one chimeric alignment with two parts.
    if not (
        ((n_algns1 == 1) and (n_algns2 == 2))
        or ((n_algns1 == 2) and (n_algns2 == 1))
    ):
        return None

    first_read_is_chimeric = n_algns1 > 1
    chim5_algn  = algns1[0] if first_read_is_chimeric else algns2[0]
    chim3_algn  = algns1[1] if first_read_is_chimeric else algns2[1]
    linear_algn = algns2[0] if first_read_is_chimeric else algns1[0]
    
    # the linear alignment must be uniquely mapped
    if not(linear_algn['is_mapped'] and linear_algn['is_unique']):
        return None

    can_rescue = True
    # we automatically rescue chimeric alignments with null and non-unique 
    # alignments at the 3' side
    if (chim3_algn['is_mapped'] and chim5_algn['is_unique']):
        # 1) in rescued walks, the 3' alignment of the chimeric alignment must be on
        # the same chromosome as the linear alignment on the opposite side of the
        # molecule
        can_rescue &= (chim3_algn['chrom'] == linear_algn['chrom'])

        # 2) in rescued walks, the 3' supplemental alignment of the chimeric
        # alignment and the linear alignment on the opposite side must point
        # towards each other
        can_rescue &= (chim3_algn['strand'] != linear_algn['strand'])
        if linear_algn['strand'] == '+':
            can_rescue &= (linear_algn['pos5'] < chim3_algn['pos5'])
        else:
            can_rescue &= (linear_algn['pos5'] > chim3_algn['pos5'])

        # 3) in single ligations appearing as walks, we can infer the size of 
        # the molecule and this size must be smaller than the maximal size of 
        # Hi-C molecules after the size selection step of the Hi-C protocol
        if linear_algn['strand'] == '+':
            molecule_size = (
                chim3_algn['pos5']
                - linear_algn['pos5']
                + chim3_algn['dist_to_5']
                + linear_algn['dist_to_5']
            )
        else:
            molecule_size = (
                linear_algn['pos5']
                - chim3_algn['pos5']
                + chim3_algn['dist_to_5']
                + linear_algn['dist_to_5']
            )

        can_rescue &= (molecule_size <= max_molecule_size)
    
    if can_rescue:
        if first_read_is_chimeric:
            # changing the type of the 3' alignment on side 1, does not show up
            # in the output
            algns1[1]['type'] = 'X' 
            algns2[0]['type'] = 'R'
            return 2
        else:
            algns1[0]['type'] = 'R'
            # changing the type of the 3' alignment on side 2, does not show up
            # in the output
            algns2[1]['type'] = 'X'
            return 1
    else:
        return None


def _convert_gaps_into_alignments(sorted_algns, max_inter_align_gap):
    if (len(sorted_algns) == 1) and (not sorted_algns[0]['is_mapped']):
        return

    last_5_pos = 0
    for i in range(len(sorted_algns)):
        algn = sorted_algns[i]
        if (algn['dist_to_5'] - last_5_pos > max_inter_align_gap):
            new_algn = empty_alignment()
            new_algn['dist_to_5'] = last_5_pos
            new_algn['algn_read_span'] = algn['dist_to_5'] - last_5_pos
            new_algn['read_len'] = algn['read_len']
            new_algn['dist_to_3'] = new_algn['read_len'] - algn['dist_to_5']

            last_5_pos = algn['dist_to_5'] + algn['algn_read_span']

            sorted_algns.insert(i, new_algn)
            i += 2
        else:
            last_5_pos = max(last_5_pos, 
                             algn['dist_to_5'] + algn['algn_read_span'])
            i += 1


def _mask_alignment(algn):
    """
    Reset the coordinates of an alignment.
    """
    algn['chrom'] = _pairsam_format.UNMAPPED_CHROM
    algn['pos5'] = _pairsam_format.UNMAPPED_POS
    algn['pos3'] = _pairsam_format.UNMAPPED_POS
    algn['pos'] = _pairsam_format.UNMAPPED_POS
    algn['strand'] = _pairsam_format.UNMAPPED_STRAND

    return algn


def parse_sams_into_pair(sams1, 
                         sams2, 
                         min_mapq, 
                         max_molecule_size, 
                         max_inter_align_gap,
                         walks_policy,
                         report_3_alignment_end,
                         sam_tags,
                         store_seq):
    """
    Parse sam entries corresponding to a Hi-C molecule into alignments
    for a Hi-C pair. 

    Returns
    -------
    algn1, algn2: dict
        Two alignments selected for reporting as a Hi-C pair.

    algns1, algns2
        All alignments, sorted according to their order in on a read.

    """

    # Check if there is at least one SAM entry per side: 
    if (len(sams1) == 0) or (len(sams2) == 0):
        algns1 = [empty_alignment()]
        algns2 = [empty_alignment()]
        algns1[0]['type'] = 'X'
        algns2[0]['type'] = 'X'
        return algns1[0], algns2[0], algns1, algns2

    # Generate a sorted, gap-filled list of all alignments
    algns1 = [parse_algn(sam.rstrip().split('\t'), min_mapq, 
                         report_3_alignment_end, sam_tags, store_seq)
              for sam in sams1]
    algns2 = [parse_algn(sam.rstrip().split('\t'), min_mapq, 
                         report_3_alignment_end, sam_tags, store_seq)
              for sam in sams2]
    algns1 = sorted(algns1, key=lambda algn: algn['dist_to_5'])
    algns2 = sorted(algns2, key=lambda algn: algn['dist_to_5'])

    if max_inter_align_gap is not None:
        _convert_gaps_into_alignments(algns1, max_inter_align_gap)
        _convert_gaps_into_alignments(algns2, max_inter_align_gap)

    # Define the type of alignment on each side.
    # The most important split is between chimeric alignments and linear 
    # alignments.
    
    is_chimeric_1 = len(algns1) > 1
    is_chimeric_2 = len(algns2) > 1

    hic_algn1 = algns1[0]
    hic_algn2 = algns2[0]

    # Parse chimeras
    rescued_linear_side = None
    if is_chimeric_1 or is_chimeric_2:
        # Pick two alignments to report as a Hi-C pair.
        rescued_linear_side = rescue_walk(algns1, algns2, max_molecule_size)

        # if the walk is unrescueable:
        if rescued_linear_side is None:
            if walks_policy == 'mask':
                if is_chimeric_1 or is_chimeric_2:
                    hic_algn1 = _mask_alignment(dict(hic_algn1))
                    hic_algn2 = _mask_alignment(dict(hic_algn2))
                    hic_algn1['type'] = 'W'
                    hic_algn2['type'] = 'W'

            elif walks_policy == '5any':
                hic_algn1 = algns1[0]
                hic_algn2 = algns2[0]

            elif walks_policy == '5unique':
                hic_algn1 = algns1[0]
                for algn in algns1:
                    if algn['is_mapped'] and algn['is_unique']:
                        hic_algn1 = algn
                        break

                hic_algn2 = algns2[0]
                for algn in algns2:
                    if algn['is_mapped'] and algn['is_unique']:
                        hic_algn2 = algn
                        break

            elif walks_policy == '3any':
                hic_algn1 = algns1[-1]
                hic_algn2 = algns2[-1]

            elif walks_policy == '3unique':
                hic_algn1 = algns1[-1]
                for algn in algns1[::-1]:
                    if algn['is_mapped'] and algn['is_unique']:
                        hic_algn1 = algn
                        break

                hic_algn2 = algns2[-1]
                for algn in algns2[::-1]:
                    if algn['is_mapped'] and algn['is_unique']:
                        hic_algn2 = algn
                        break
            
            elif walks_policy == 'all':
                # TO BE IMPLEMENTED
                pass


            # lower-case reported walks on the chimeric side
            if walks_policy != 'mask':
                if is_chimeric_1:
                    hic_algn1 = dict(hic_algn1)
                    hic_algn1['type'] = hic_algn1['type'].lower()
                if is_chimeric_2:
                    hic_algn2 = dict(hic_algn2)
                    hic_algn2['type'] = hic_algn2['type'].lower()


    return hic_algn1, hic_algn2, algns1, algns2

def check_pair_order(algn1, algn2, chrom_enum):
    '''
    Check if a pair of alignments has the upper-triangular order or
    has to be flipped.
    '''

    # First, the pair is flipped according to the type of mapping on its sides.
    # Later, we will check it is mapped on both sides and, if so, flip the sides
    # according to these coordinates.

    has_correct_order = (
           (algn1['is_mapped'], algn1['is_unique']) 
        <= (algn2['is_mapped'], algn2['is_unique'])
    )

    # If a pair has coordinates on both sides, it must be flipped according to
    # its genomic coordinates.
    if ((algn1['chrom'] != _pairsam_format.UNMAPPED_CHROM)
        and (algn2['chrom'] != _pairsam_format.UNMAPPED_CHROM)):
        
        has_correct_order = (
                (chrom_enum[algn1['chrom']], algn1['pos']) 
             <= (chrom_enum[algn2['chrom']], algn2['pos']))

    return has_correct_order

def push_sam(line, drop_seq, sams1, sams2):
    """

    """

    sam = line.rstrip()
    if drop_seq:
        split_sam = sam.split('\t')
        split_sam[9] = '*'
        split_sam[10] = '*'
        sam = '\t'.join(split_sam)

        flag = split_sam[1]
        flag = int(flag)
    else:
        _, flag, _ = sam.split('\t', 2)
        flag = int(flag)


    if ((flag & 0x40) != 0):
        sams1.append(sam)
    else:
        sams2.append(sam)
    return


def write_all_algnments(read_id, all_algns1, all_algns2, out_file):
    for side_idx, all_algns in enumerate((all_algns1, all_algns2)):
        out_file.write(read_id)
        out_file.write('\t')
        out_file.write(str(side_idx+1))
        out_file.write('\t')
        for algn in sorted(all_algns, key=lambda x: x['dist_to_5']):
            out_file.write(algn['chrom'])
            out_file.write('\t')
            out_file.write(str(algn['pos5']))
            out_file.write('\t')
            out_file.write(algn['strand'])
            out_file.write('\t')
            out_file.write(str(algn['mapq']))
            out_file.write('\t')
            out_file.write(str(algn['cigar']))
            out_file.write('\t')
            out_file.write(str(algn['dist_to_5']))
            out_file.write('\t')
            out_file.write(str(algn['dist_to_5']+algn['algn_read_span']))
            out_file.write('\t')
            out_file.write(str(algn['matched_bp']))
            out_file.write('\t')

        out_file.write('\n')

def write_pairsam(
        algn1, algn2, read_id, sams1, sams2, out_file, 
        drop_readid, drop_sam, add_columns):
    """
    SAM is already tab-separated and
    any printable character between ! and ~ may appear in the PHRED field!
    (http://www.ascii-code.com/)
    Thus, use the vertical tab character to separate fields!

    """
    cols = [
        '.' if drop_readid else read_id,
        algn1['chrom'],
        str(algn1['pos']),
        algn2['chrom'],
        str(algn2['pos']),
        algn1['strand'],
        algn2['strand'],
        algn1['type'] + algn2['type']
    ]

    if not drop_sam:
        for sams in [sams1, sams2]:
            cols.append(
                _pairsam_format.INTER_SAM_SEP.join([
                    (sam.replace('\t', _pairsam_format.SAM_SEP)
                    + _pairsam_format.SAM_SEP 
                    + 'Yt:Z:' + algn1['type'] + algn2['type'])
                for sam in sams
                ])
            )

    
    for col in add_columns:
        # use get b/c empty alignments would not have sam tags (NM, AS, etc)
        cols.append(str(algn1.get(col, '')))
        cols.append(str(algn2.get(col, '')))

    out_file.write(_pairsam_format.PAIRSAM_SEP.join(cols) + '\n')


def streaming_classify(instream, outstream, chromosomes, min_mapq, max_molecule_size, 
                       drop_readid, drop_seq, drop_sam, add_columns, 
                       out_alignments_stream, out_stat, **kwargs):
    """

    """
    chrom_enum = dict(zip([_pairsam_format.UNMAPPED_CHROM] + list(chromosomes), 
                          range(len(chromosomes)+1)))
    sam_tags = [col for col in add_columns if len(col)==2 and col.isupper()]
    prev_read_id = ''
    sams1 = []
    sams2 = []
    line = ''
    store_seq = ('seq' in add_columns)
    
    instream = iter(instream)
    while line is not None:
        line = next(instream, None)

        read_id = line.split('\t', 1)[0] if line else None

        if not(line) or ((read_id != prev_read_id) and prev_read_id):
            algn1, algn2, all_algns1, all_algns2 = parse_sams_into_pair(
                sams1,
                sams2,
                min_mapq,
                max_molecule_size, 
                kwargs['max_inter_align_gap'],
                kwargs['walks_policy'],
                kwargs['report_alignment_end']=='3',
                sam_tags,
                store_seq
                )

            flip_pair = (not kwargs['no_flip']) and (
                    not check_pair_order(algn1, algn2, chrom_enum))

            if flip_pair:
                algn1, algn2 = algn2, algn1
                sams1, sams2 = sams2, sams1

            write_pairsam(
                algn1, algn2,
                prev_read_id, 
                sams1, sams2,
                outstream,
                drop_readid,
                drop_sam,
                add_columns)

            # add a pair to PairCounter if stats output is requested:
            if out_stat:
                out_stat.add_pair(algn1['chrom'],  int(algn1['pos']),  algn1['strand'],
                                  algn2['chrom'],  int(algn2['pos']),  algn2['strand'],
                                  algn1['type'] + algn2['type'])

            if out_alignments_stream:
                write_all_algnments(prev_read_id, all_algns1, all_algns2, out_alignments_stream)
            
            sams1.clear()
            sams2.clear()

        if line is not None:
            push_sam(line, drop_seq, sams1, sams2)
            prev_read_id = read_id

if __name__ == '__main__':
    parse()

def parse_alternative_algns(samcols):
    alt_algns = []
    for col in samcols[11:]:
        if not col.startswith('XA:Z:'):
            continue

        for SA in col[5:].split(';'):
            if not SA:
                continue
            SAcols = SA.split(',')

            chrom = SAcols[0] 
            strand = '-' if SAcols[1]<0 else '+'

            cigar = parse_cigar(SAcols[2])
            NM = SAcols[3] 

            pos = _pairsam_format.UNMAPPED_POS
            if strand == '+':
                pos = int(SAcols[1])
            else:
                pos = int(SAcols[1]) + cigar['algn_ref_span']

            alt_algns.append({
                'chrom': chrom,
                'pos': pos,
                'strand': strand,
                'mapq': mapq,
                'is_mapped': True,
                'is_unique': False,
                'is_linear': None,
                'cigar': cigar,
                'NM': NM,
                'dist_to_5': cigar['clip5_ref'] if strand == '+' else cigar['clip3_ref'],
            })

    return supp_algns

#def parse_supp(samcols, min_mapq):
#    supp_algns = []
#    for col in samcols[11:]:
#        if not col.startswith('SA:Z:'):
#            continue
#
#        for SA in col[5:].split(';'):
#            if not SA:
#                continue
#            SAcols = SA.split(',')
#            mapq = int(SAcols[4])
#            is_unique = (mapq >= min_mapq)
#
#            chrom = SAcols[0] if is_unique else _pairsam_format.UNMAPPED_CHROM
#            strand = SAcols[2] if is_unique else _pairsam_format.UNMAPPED_STRAND
#
#            cigar = parse_cigar(SAcols[3])
#
#            pos = _pairsam_format.UNMAPPED_POS
#            if is_unique:
#                if strand == '+':
#                    pos = int(SAcols[1])
#                else:
#                    pos = int(SAcols[1]) + cigar['algn_ref_span']
#
#            supp_algns.append({
#                'chrom': chrom,
#                'pos': pos,
#                'strand': strand,
#                'mapq': mapq,
#                'is_mapped': True,
#                'is_unique': is_unique,
#                'is_linear': None,
#                'cigar': cigar,
#                'dist_to_5': cigar['clip5_ref'] if strand == '+' else cigar['clip3_ref'],
#            })
#
#    return supp_algns
