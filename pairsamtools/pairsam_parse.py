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

from . import _io, _pairsam_format, _headerops, cli

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
         'If the path ends with .gz, the output is bgzip-compressed. '
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
        ' If the path ends with .gz, the output is bgzip-compressed.'
        ' By default, not used.'
        )
def parse(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
          drop_readid, drop_seq, drop_sam, store_mapq, store_unrescuable_chimeras,
          output_parsed_alignments):
    '''parse .sam and make .pairsam.

    SAM_PATH : input .sam file. If the path ends with .bam, the input is 
    decompressed from bam. By default, the input is read from stdin.
    '''
    parse_py(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
             drop_readid, drop_seq, drop_sam, store_mapq, store_unrescuable_chimeras,
             output_parsed_alignments)


def parse_py(sam_path, chroms_path, output, assembly, min_mapq, max_molecule_size, 
             drop_readid, drop_seq, drop_sam, store_mapq, 
             store_unrescuable_chimeras, output_parsed_alignments):
    instream = (_io.open_sam_or_bam(sam_path, mode='r') 
                if sam_path else sys.stdin)
    outstream = (_io.open_bgzip(output, mode='w') 
                 if output else sys.stdout)
    out_alignments_stream = (_io.open_bgzip(output_parsed_alignments, mode='w') 
                 if output_parsed_alignments else None)
    if out_alignments_stream:
        out_alignments_stream.write('side\tchrom\tpos\tstrand\tmapq\tcigar\tdist_5_lo\tdist_5_hi\tmatched_bp\n')


    samheader, body_stream = _headerops.get_header(instream, comment_char='@')
    sam_chromsizes = _headerops.get_chromsizes_from_sam_header(samheader)
    chromosomes = _headerops.get_chrom_order(
        chroms_path, 
        list(sam_chromsizes.keys()))

    header = _headerops.make_standard_pairsheader(
        assembly=assembly,
        chromsizes=[(chrom, sam_chromsizes[chrom]) for chrom in chromosomes],
        columns=_pairsam_format.COLUMNS + (['mapq1', 'mapq2'] if store_mapq else [])
    )
    
    header = _headerops.insert_samheader(header, samheader) 
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    outstream.writelines((l+'\n' for l in header))

    streaming_classify(body_stream, outstream, chromosomes, min_mapq, 
                       max_molecule_size, drop_readid, drop_seq, drop_sam, 
                       store_mapq, store_unrescuable_chimeras, out_alignments_stream)

    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


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
    pos = _pairsam_format.UNMAPPED_POS
    if is_mapped and is_unique: 
        if ((int(samcols[1]) & 0x10) == 0):
            strand = '+'
            pos = int(samcols[3])
        else:
            strand = '-'
            pos = int(samcols[3]) + cigar['algn_ref_span']

    return {
        'chrom': chrom,
        'pos': pos,
        'strand': strand,
        'mapq': mapq,
        'is_mapped': is_mapped,
        'is_unique': is_unique,
        'is_linear': is_linear,
        'cigar': cigar,
        'dist_to_5': cigar['clip5'] if strand == '+' else cigar['clip3'],
    }


def parse_supp(samcols, min_mapq):
    supp_algns = []
    for col in samcols[11:]:
        if not col.startswith('SA:Z:'):
            continue

        for SA in col[5:].split(';'):
            if not SA:
                continue
            SAcols = SA.split(',')
            mapq = int(SAcols[4])
            is_unique = (mapq >= min_mapq)

            chrom = SAcols[0] if is_unique else _pairsam_format.UNMAPPED_CHROM
            strand = SAcols[2] if is_unique else _pairsam_format.UNMAPPED_STRAND

            cigar = parse_cigar(SAcols[3])

            pos = _pairsam_format.UNMAPPED_POS
            if is_unique:
                if strand == '+':
                    pos = int(SAcols[1])
                else:
                    pos = int(SAcols[1]) + cigar['algn_ref_span']

            supp_algns.append({
                'chrom': chrom,
                'pos': pos,
                'strand': strand,
                'mapq': mapq,
                'is_mapped': True,
                'is_unique': is_unique,
                'is_linear': None,
                'cigar': cigar,
                'dist_to_5': cigar['clip5'] if strand == '+' else cigar['clip3'],
            })

    return supp_algns


def rescue_chimeric_alignment(repr_algn1, repr_algn2, supp_algns1, supp_algns2,
                              max_molecule_size):
    """

    Returns
    -------
    algn1_5, algn2_5, is_rescued

    """
    # If both alignments are non-chimeric, no need to rescue!
    if (not supp_algns1) and (not supp_algns2):
        return True

    # If both alignments are chimeric, cannot rescue
    if (supp_algns1) and (supp_algns2):
        return False

    # Cannot rescue a chimeric alignment with multiple supplemental alignments.
    if (len(supp_algns1) > 1) or (len(supp_algns2) > 1):
        return False

    sup_algn = supp_algns1[0] if (supp_algns1) else supp_algns2[0]
    # if the supplemental alignment is non-unique, no need to rescue!
    if not sup_algn['is_unique']:
        return True

    first_read_is_chimeric = bool(supp_algns1)
    linear_algn = repr_algn2 if first_read_is_chimeric else repr_algn1
    repr_algn = repr_algn1 if first_read_is_chimeric else repr_algn2

    chim5_algn = (repr_algn
                  if (repr_algn['dist_to_5'] < sup_algn['dist_to_5'])
                  else sup_algn)
    chim3_algn = (sup_algn
                  if (repr_algn['dist_to_5'] < sup_algn['dist_to_5'])
                  else repr_algn)

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
        can_rescue &= (linear_algn['pos'] < chim3_algn['pos'])
    else:
        can_rescue &= (linear_algn['pos'] > chim3_algn['pos'])

    # for normal chimeras, we can infer the size of the molecule and
    # this size must be smaller than the maximal size of Hi-C molecules after
    # the size selection step of the Hi-C protocol
    if linear_algn['strand'] == '+':
        molecule_size = (
            chim3_algn['pos']
            - linear_algn['pos']
            + chim3_algn['dist_to_5']
            + linear_algn['dist_to_5']
        )
    else:
        molecule_size = (
            linear_algn['pos']
            - chim3_algn['pos']
            + chim3_algn['dist_to_5']
            + linear_algn['dist_to_5']
        )

    can_rescue &= (molecule_size <= max_molecule_size)

    return can_rescue


def parse_sams_into_pair(chromosomes, sams1, sams2, min_mapq, max_molecule_size, 
             store_unrescuable_chimeras):
    """
    Possible pair types:
    ...

    Returns
    -------
    pair_type, algn1, algn2, flip_pair

    """
    chrom_enum = dict(zip([_pairsam_format.UNMAPPED_CHROM] + list(chromosomes), range(len(chromosomes))))
    sam1_repr_cols = sams1[0].rstrip().split('\t')
    sam2_repr_cols = sams2[0].rstrip().split('\t')

    algn1 = parse_algn(sam1_repr_cols, min_mapq)
    algn2 = parse_algn(sam2_repr_cols, min_mapq)

    is_null_1 = not algn1['is_mapped']
    is_null_2 = not algn2['is_mapped']
    is_multi_1 = not algn1['is_unique']
    is_multi_2 = not algn2['is_unique']
    is_chimeric_1 = not algn1['is_linear']
    is_chimeric_2 = not algn2['is_linear']


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
        algn1['pos'] = _pairsam_format.UNMAPPED_POS
        algn1['strand'] = _pairsam_format.UNMAPPED_STRAND
    if is_multi_2:
        algn2['chrom'] = _pairsam_format.UNMAPPED_CHROM
        algn2['pos'] = _pairsam_format.UNMAPPED_POS
        algn2['strand'] = _pairsam_format.UNMAPPED_STRAND

    # Parse chimeras

    all_algns1 = [algn1]
    all_algns2 = [algn2]
    
    if ('C' in pair_type): 
        if (pair_type == 'CL') or store_unrescuable_chimeras:
            supp_algns1 = parse_supp(sam1_repr_cols, min_mapq)
            supp_algns2 = parse_supp(sam2_repr_cols, min_mapq)

            all_algns1 += supp_algns1
            all_algns2 += supp_algns2

            algn1_5 = min(all_algns1, key=lambda x:x['dist_to_5'])
            algn2_5 = min(all_algns2, key=lambda x:x['dist_to_5'])
    
        is_rescued = False
        if (pair_type == 'CL'):
           is_rescued = rescue_chimeric_alignment(
                algn1, algn2, supp_algns1, supp_algns2, max_molecule_size)

           if is_rescued:
               pair_type = 'CX'

        if (is_rescued or store_unrescuable_chimeras):
            algn1 = algn1_5
            algn2 = algn2_5
        else:
            if is_chimeric_1:
                algn1 = copy.deepcopy(algn1)
                algn1['chrom'] = _pairsam_format.UNMAPPED_CHROM
                algn1['pos'] = _pairsam_format.UNMAPPED_POS
                algn1['strand'] = _pairsam_format.UNMAPPED_STRAND
            if is_chimeric_2:
                algn2 = copy.deepcopy(algn2)
                algn2['chrom'] = _pairsam_format.UNMAPPED_CHROM
                algn2['pos'] = _pairsam_format.UNMAPPED_POS
                algn2['strand'] = _pairsam_format.UNMAPPED_STRAND
        
    # If a pair has coordinates on both sides, it must be flipped according to
    # its genomic coordinates.

    if ((algn1['chrom'] != _pairsam_format.UNMAPPED_CHROM)
        and (algn1['chrom'] != _pairsam_format.UNMAPPED_CHROM)):
        
        flip_pair = (
            (chrom_enum[algn1['chrom']], algn1['pos']) 
                > (chrom_enum[algn2['chrom']], algn2['pos']))


    return pair_type, algn1, algn2, flip_pair, all_algns1, all_algns2


def push_sam(line, drop_seq, sams1, sams2):
    """

    """

    if drop_seq:
        split_line = line.split('\t')
        split_line[9] = '*'
        split_line[10] = '*'
        line = '\t'.join(split_line)

        flag = split_line[1]
        flag = int(flag)
    else:
        _, flag, _ = line.split('\t', 2)
        flag = int(flag)


    if ((flag & 0x40) != 0):
        if ((flag & 0x800) == 0):
            sams1.insert(0, line)
        else:
            sams1.append(line)
    else:
        if ((flag & 0x800) == 0):
            sams2.insert(0, line)
        else:
            sams2.append(line)
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
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    if drop_sam:
        out_file.write('.')
    else:
        for i, sam in enumerate(sams1):
            out_file.write(sam[:-1].replace('\t', _pairsam_format.SAM_SEP))
            out_file.write(_pairsam_format.SAM_SEP + 'Yt:Z:')
            out_file.write(pair_type)
            if i < len(sams1) -1:
                out_file.write(_pairsam_format.INTER_SAM_SEP)
    out_file.write(_pairsam_format.PAIRSAM_SEP)
    if drop_sam:
        out_file.write('.')
    else:
        for i, sam in enumerate(sams2):
            out_file.write(sam[:-1].replace('\t', _pairsam_format.SAM_SEP))
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
                       store_unrescuable_chimeras, out_alignments_stream):
    """

    """

    prev_read_id = ''
    sams1 = []
    sams2 = []
    line = ''
    
    while line is not None:
        line = next(instream, None)

        read_id = line.split('\t', 1)[0] if line else None

        if not(line) or ((read_id != prev_read_id) and prev_read_id):
            pair_type, algn1, algn2, flip_pair, all_algns1, all_algns2 = parse_sams_into_pair(
                chromosomes,
                sams1,
                sams2,
                min_mapq,
                max_molecule_size, 
                store_unrescuable_chimeras)
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
            if out_alignments_stream:
                write_all_algnments(all_algns1, all_algns2, out_alignments_stream)
            
            sams1.clear()
            sams2.clear()

        if line is not None:
            push_sam(line, drop_seq, sams1, sams2)
            prev_read_id = read_id

if __name__ == '__main__':
    parse()
