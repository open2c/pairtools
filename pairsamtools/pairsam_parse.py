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
import copy

from . import _fileio, _pairsam_format, _headerops, cli, common_io_options
from .pairsam_stats import PairCounter


UTIL_NAME = 'pairsam_parse'


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
    help='Name of genome assembly (e.g. hg19, mm10)')
@click.option(
    "--min-mapq", 
    type=int, 
    default=1,
    show_default=True,
    help='The minimal MAPQ score of a mapped read')
@click.option(
    "--max-molecule-size", 
    type=int, 
    default=2000,
    show_default=True,
    help='The maximal size of a ligated Hi-C molecule; used in chimera rescue.')
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
    "--store-mapq", 
    is_flag=True,
    help='If specified, add two columns with MAPQ on both sides.')
@click.option(
    "--store-unrescuable-chimeras", 
    is_flag=True,
    help='If specified, do not mask unrescuable chimeras and report the 5\'-most alignments.')
@click.option(
    "--output-parsed-alignments", 
    type=str, 
    default="", 
    help='output file for all parsed alignments, including chimeras.'
        ' Useful for debugging and rnalysis of chimeras.'
        ' If file exists, it will be open in the append mode.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4-compressed.'
        ' By default, not used.'
        )
@click.option(
    "--output-stats", 
    type=str, 
    default="", 
    help='output file for various statistics of pairsam file. '
        ' By default, statistics is not generated.')
@click.option(
    '--report-alignment-end', 
    type=click.Choice(['5', '3']),
    default='5',
    help='specifies whether the 5\' or 3\' end of the alignment is reported as'
    ' the position of the Hi-C read.'

)

@common_io_options

def parse(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
          drop_readid, drop_seq, drop_sam, store_mapq, store_unrescuable_chimeras,
          output_parsed_alignments, output_stats, **kwargs):
    '''parse .sam and make .pairsam.

    SAM_PATH : input .sam file. If the path ends with .bam, the input is 
    decompressed from bam. By default, the input is read from stdin.
    '''
    parse_py(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
             drop_readid, drop_seq, drop_sam, store_mapq, store_unrescuable_chimeras,
             output_parsed_alignments, output_stats, **kwargs)


def parse_py(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
             drop_readid, drop_seq, drop_sam, store_mapq, 
             store_unrescuable_chimeras, output_parsed_alignments, output_stats, **kwargs):
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
        out_alignments_stream.write('side\tchrom\tpos\tstrand\tmapq\tcigar\tdist_5_lo\tdist_5_hi\tmatched_bp\n')

    # generate empty PairCounter if stats output is requested:
    out_stat = PairCounter() if output_stats else None


    samheader, body_stream = _headerops.get_header(instream, comment_char='@')
    sam_chromsizes = _headerops.get_chromsizes_from_sam_header(samheader)
    chromosomes = _headerops.get_chrom_order(
        chroms_path, 
        list(sam_chromsizes.keys()))

    columns =  (_pairsam_format.COLUMNS 
                + (['mapq1', 'mapq2'] if store_mapq else []))

    if drop_sam:
        columns.pop(columns.index('sam1'))
        columns.pop(columns.index('sam2'))

    header = _headerops.make_standard_pairsheader(
        assembly = assembly,
        chromsizes = [(chrom, sam_chromsizes[chrom]) for chrom in chromosomes],
        columns = columns
    )
    
    header = _headerops.insert_samheader(header, samheader) 
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l+'\n' for l in header))

    streaming_classify(body_stream, outstream, chromosomes, min_mapq, 
                       max_molecule_size, drop_readid, drop_seq, drop_sam, 
                       store_mapq, store_unrescuable_chimeras,
                       out_alignments_stream, out_stat, **kwargs)

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
    clip5 = 0
    clip3 = 0

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
                        clip5 = cur_num
                    else:
                        clip3 = cur_num

                cur_num = 0

    return {
        'clip5': clip5,
        'clip3': clip3,
        'cigar_str': cigar,
        'algn_ref_span': algn_ref_span,
        'algn_read_span': algn_read_span,
        'read_len': read_len,
        'matched_bp': matched_bp,
    }


def parse_algn(samcols, min_mapq):
    is_mapped = (int(samcols[1]) & 0x04) == 0

    mapq = int(samcols[4])
    is_unique = (mapq >= min_mapq)

    is_linear = not any([col.startswith('SA:Z:') for col in samcols[11:]])

    cigar = parse_cigar(samcols[5])

    chrom = samcols[2] if (is_mapped and is_unique) else _pairsam_format.UNMAPPED_CHROM

    strand = _pairsam_format.UNMAPPED_STRAND
    pos5 = _pairsam_format.UNMAPPED_POS
    pos3 = _pairsam_format.UNMAPPED_POS
    if is_mapped and is_unique: 
        if ((int(samcols[1]) & 0x10) == 0):
            strand = '+'
            pos5 = int(samcols[3])
            pos3 = int(samcols[3]) + cigar['algn_ref_span']
        else:
            strand = '-'
            pos5 = int(samcols[3]) + cigar['algn_ref_span']
            pos3 = int(samcols[3])

    return {
        'chrom': chrom,
        'pos5': pos5,
        'pos3': pos3,
        'strand': strand,
        'mapq': mapq,
        'is_mapped': is_mapped,
        'is_unique': is_unique,
        'is_linear': is_linear,
        'cigar': cigar,
        'dist_to_5': cigar['clip5'] if strand == '+' else cigar['clip3'],
    }


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
#                'dist_to_5': cigar['clip5'] if strand == '+' else cigar['clip3'],
#            })
#
#    return supp_algns


def rescue_chimeric_alignment(algns1, algns2, max_molecule_size):
    """

    Returns
    -------
    algn1_5, algn2_5, is_rescued

    """
    
    # If both alignments are non-chimeric, no need to rescue!
    n_unique_algns1 = sum([a['is_mapped'] and a['is_unique'] for a in algns1])
    n_unique_algns2 = sum([a['is_mapped'] and a['is_unique'] for a in algns2]) 

    assert (n_unique_algns1 >= 1) and (n_unique_algns2 >= 1)

    if (n_unique_algns1 == 1) and (n_unique_algns2 == 1):
        return True

    # Can rescue only pairs with one chimeric alignment with two parts.
    if not (
        ((n_unique_algns1 == 1) and (n_unique_algns2 == 2))
        or ((n_unique_algns1 == 2) and (n_unique_algns2 == 1))
    ):
        return False

    first_read_is_chimeric = n_unique_algns1 > 1
    chim5_algn  = algns1[0] if first_read_is_chimeric else algns2[0]
    chim3_algn  = algns1[1] if first_read_is_chimeric else algns2[1]
    linear_algn = algns2[0] if first_read_is_chimeric else algns1[0]

    can_rescue = True
    # in normal chimeras, the 3' alignment of the chimeric alignment must be on
    # the same chromosome as the linear alignment on the opposite side of the
    # molecule
    can_rescue &= (chim3_algn['chrom'] == linear_algn['chrom'])

    # in normal chimeras, the 3' supplemental alignment of the chimeric
    # alignment and the linear alignment on the opposite side must point
    # towards each other
    can_rescue &= (chim3_algn['strand'] != linear_algn['strand'])
    if linear_algn['strand'] == '+':
        can_rescue &= (linear_algn['pos5'] < chim3_algn['pos5'])
    else:
        can_rescue &= (linear_algn['pos5'] > chim3_algn['pos5'])

    # for normal chimeras, we can infer the size of the molecule and
    # this size must be smaller than the maximal size of Hi-C molecules after
    # the size selection step of the Hi-C protocol
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

    return can_rescue


def _alignment_sorting_order(algn):
    if (algn['is_mapped'] and algn['is_unique']):
        return algn['dist_to_5'] 
    else:
        return algn['cigar']['read_len'] + algn['dist_to_5']

def parse_sams_into_pair(chrom_enum, sams1, sams2, min_mapq, max_molecule_size, 
                         store_unrescuable_chimeras, report_alignment_end):
    """
    Possible pair types:
    ...

    Returns
    -------
    pair_type, algn1, algn2, flip_pair

    """
    sam1_repr_cols = sams1[0].rstrip().split('\t')
    sam2_repr_cols = sams2[0].rstrip().split('\t')

    algns1 = [parse_algn(sam.rstrip().split('\t'), min_mapq)
              for sam in sams1]
    algns2 = [parse_algn(sam.rstrip().split('\t'), min_mapq)
              for sam in sams2]
    algns1 = sorted(algns1, key=_alignment_sorting_order)
    algns2 = sorted(algns2, key=_alignment_sorting_order)

    algn1 = algns1[0]
    algn2 = algns2[0]

    is_null_1 = not algn1['is_mapped']
    is_null_2 = not algn2['is_mapped']
    is_multi_1 = not algn1['is_unique']
    is_multi_2 = not algn2['is_unique']
    is_chimeric_1 = sum([a['is_mapped'] and a['is_unique'] for a in algns1]) > 1
    is_chimeric_2 = sum([a['is_mapped'] and a['is_unique'] for a in algns2]) > 1

    side_type_1 = ('N' if is_null_1 
                       else ('M' if is_multi_1 
                                 else ('C' if is_chimeric_1
                                           else 'L')))
    side_type_2 = ('N' if is_null_2 
                       else ('M' if is_multi_2 
                                 else ('C' if is_chimeric_2
                                           else 'L')))

    # First, the pair is flipped according to the type of mapping on its sides.
    # Later, we will check it is mapped on both sides and, if so, flip the sides
    # according to these coordinates.

    flip_pair = (is_null_1, is_multi_1, is_chimeric_1) < (is_null_2, is_multi_2, is_chimeric_2)

    pair_type = (side_type_2+side_type_1) if flip_pair else (side_type_1+side_type_2)

    if is_multi_1:
        algn1['chrom'] = _pairsam_format.UNMAPPED_CHROM
        algn1['pos5'] = _pairsam_format.UNMAPPED_POS
        algn1['pos3'] = _pairsam_format.UNMAPPED_POS
        algn1['strand'] = _pairsam_format.UNMAPPED_STRAND
    if is_multi_2:
        algn2['chrom'] = _pairsam_format.UNMAPPED_CHROM
        algn2['pos5'] = _pairsam_format.UNMAPPED_POS
        algn2['pos3'] = _pairsam_format.UNMAPPED_POS
        algn2['strand'] = _pairsam_format.UNMAPPED_STRAND

    # Parse chimeras
    
    if ('C' in pair_type): 
       is_rescued = (pair_type == 'CL') and rescue_chimeric_alignment(
           algns1, algns2, max_molecule_size)

       if is_rescued:
           pair_type = 'CX'

       if not (is_rescued or store_unrescuable_chimeras):
           if is_chimeric_1:
               algn1 = copy.deepcopy(algn1)
               algn1['chrom'] = _pairsam_format.UNMAPPED_CHROM
               algn1['pos5'] = _pairsam_format.UNMAPPED_POS
               algn1['pos3'] = _pairsam_format.UNMAPPED_POS
               algn1['strand'] = _pairsam_format.UNMAPPED_STRAND
           if is_chimeric_2:
               algn2 = copy.deepcopy(algn2)
               algn2['chrom'] = _pairsam_format.UNMAPPED_CHROM
               algn2['pos5'] = _pairsam_format.UNMAPPED_POS
               algn2['pos3'] = _pairsam_format.UNMAPPED_POS
               algn2['strand'] = _pairsam_format.UNMAPPED_STRAND
        
    # If a pair has coordinates on both sides, it must be flipped according to
    # its genomic coordinates.

    if report_alignment_end == '5':
        algn1['pos'] = algn1['pos5']
        algn2['pos'] = algn2['pos5']
    else:
        algn1['pos'] = algn1['pos3']
        algn2['pos'] = algn2['pos3']

    if ((algn1['chrom'] != _pairsam_format.UNMAPPED_CHROM)
        and (algn1['chrom'] != _pairsam_format.UNMAPPED_CHROM)):
        
        flip_pair = (
            (chrom_enum[algn1['chrom']], algn1['pos']) 
                > (chrom_enum[algn2['chrom']], algn2['pos']))


    return pair_type, algn1, algn2, flip_pair, algns1, algns2


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


def write_all_algnments(all_algns1, all_algns2, out_file):
    for side_idx, all_algns in enumerate((all_algns1, all_algns2)):
        out_file.write(str(side_idx+1))
        out_file.write('\t')
        for algn in sorted(all_algns, key=lambda x: x['dist_to_5']):
            out_file.write(algn['chrom'])
            out_file.write('\t')
            out_file.write(str(algn['pos']))
            out_file.write('\t')
            out_file.write(algn['strand'])
            out_file.write('\t')
            out_file.write(str(algn['mapq']))
            out_file.write('\t')
            out_file.write(str(algn['cigar']['cigar_str']))
            out_file.write('\t')
            out_file.write(str(algn['dist_to_5']))
            out_file.write('\t')
            out_file.write(str(algn['dist_to_5']+algn['cigar']['algn_read_span']))
            out_file.write('\t')
            out_file.write(str(algn['cigar']['matched_bp']))
            out_file.write('\t')

        out_file.write('\n')

def write_pairsam(
        algn1, algn2, read_id, pair_type, sams1, sams2, out_file, 
        drop_readid, drop_sam, store_mapq):
    """
    SAM is already tab-separated and
    any printable character between ! and ~ may appear in the PHRED field!
    (http://www.ascii-code.com/)
    Thus, use the vertical tab character to separate fields!

    """
    if drop_readid:
        out_file.write('.')
    else:
        out_file.write(read_id)
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    out_file.write(algn1['chrom'])
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    out_file.write(str(algn1['pos']))
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    out_file.write(algn2['chrom'])
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    out_file.write(str(algn2['pos']))
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    out_file.write(algn1['strand'])
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    out_file.write(algn2['strand'])
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    out_file.write(pair_type)

    if not drop_sam:
        out_file.write(_pairsam_format.PAIRSAM_SEP)
        for i, sam in enumerate(sams1):
            out_file.write(sam.replace('\t', _pairsam_format.SAM_SEP))
            out_file.write(_pairsam_format.SAM_SEP + 'Yt:Z:')
            out_file.write(pair_type)
            if i < len(sams1) -1:
                out_file.write(_pairsam_format.INTER_SAM_SEP)

        out_file.write(_pairsam_format.PAIRSAM_SEP)
        for i, sam in enumerate(sams2):
            out_file.write(sam.replace('\t', _pairsam_format.SAM_SEP))
            out_file.write(_pairsam_format.SAM_SEP + 'Yt:Z:')
            out_file.write(pair_type)
            if i < len(sams2) -1:
                out_file.write(_pairsam_format.INTER_SAM_SEP)

    if store_mapq:
        out_file.write(_pairsam_format.PAIRSAM_SEP)
        out_file.write(str(algn1['mapq']))
        out_file.write(_pairsam_format.PAIRSAM_SEP)
        out_file.write(str(algn2['mapq']))

    out_file.write('\n')


def streaming_classify(instream, outstream, chromosomes, min_mapq, max_molecule_size, 
                       drop_readid, drop_seq, drop_sam, store_mapq, 
                       store_unrescuable_chimeras, out_alignments_stream, out_stat, **kwargs):
    """

    """
    chrom_enum = dict(zip([_pairsam_format.UNMAPPED_CHROM] + list(chromosomes), 
                          range(len(chromosomes)+1)))
    prev_read_id = ''
    sams1 = []
    sams2 = []
    line = ''
    
    instream = iter(instream)
    while line is not None:
        line = next(instream, None)

        read_id = line.split('\t', 1)[0] if line else None

        if not(line) or ((read_id != prev_read_id) and prev_read_id):
            pair_type, algn1, algn2, flip_pair, all_algns1, all_algns2 = parse_sams_into_pair(
                chrom_enum,
                sams1,
                sams2,
                min_mapq,
                max_molecule_size, 
                store_unrescuable_chimeras, 
                kwargs['report_alignment_end'])
            if flip_pair:
                write_pairsam(
                    algn2, algn1, 
                    prev_read_id, 
                    pair_type,
                    sams2, sams1,
                    outstream, 
                    drop_readid,
                    drop_sam,
                    store_mapq)
                # add a pair to PairCounter if stats output requested
                if out_stat:
                    out_stat.add_pair(algn2, algn1, pair_type)
            else:
                write_pairsam(
                    algn1, algn2,
                    prev_read_id, 
                    pair_type,
                    sams1, sams2,
                    outstream,
                    drop_readid,
                    drop_sam,
                    store_mapq)
                # add a pair to PairCounter if stats output requested
                if out_stat:
                    out_stat.add_pair(algn1, algn2, pair_type)
            if out_alignments_stream:
                write_all_algnments(all_algns1, all_algns2, out_alignments_stream)
            
            sams1.clear()
            sams2.clear()

        if line is not None:
            push_sam(line, drop_seq, sams1, sams2)
            prev_read_id = read_id

if __name__ == '__main__':
    parse()
