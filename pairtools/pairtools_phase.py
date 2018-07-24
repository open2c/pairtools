import sys
import click
import re, fnmatch

from . import _fileio, _pairsam_format, cli, _headerops, common_io_options

UTIL_NAME = 'pairtools_phase'

@cli.command()
@click.argument(
    'pairs_path', 
    type=str,
    required=False)

@click.option(
    '-o', "--output", 
    type=str, 
    default="", 
    help='output file.'
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4c-compressed.'
        ' By default, the output is printed into stdout.')

@click.option(
    "--phase-suffixes", 
    nargs=2,
    #type=click.Tuple([str, str]),
    help='phase suffixes.'
    )

@click.option(
    "--clean-output", 
    is_flag=True,
    help='drop all columns besides the standard ones and phase1/2'
    )

@common_io_options

def phase(
    pairsam_path,
    output,
    phase_suffixes,
    clean_output,
    **kwargs
    ):
    '''Phase pairs mapped to a diploid genome.

    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz or .lz4, the
    input is decompressed by pbgzip/lz4c. By default, the input is read from stdin.

    '''
    phase_py(
        pairs_path, output, phase_suffixes, clean_output,
        **kwargs
    )

    
def phase_py(
    pairs_path, output, phase_suffixes, clean_output,
    **kwargs
    ):

    instream = (_fileio.auto_open(pairsam_path, mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if pairsam_path else sys.stdin)
    outstream = (_fileio.auto_open(output, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)
    old_column_names = _headerops.extract_column_names(header)

    if clean_output:
        new_column_names = [col for col in old_column_names
                            if col in _pairsam_format.COLUMNS]
        new_column_idxs = [i for i,col in enumerate(old_column_names)
                            if col in _pairsam_format.COLUMNS] + [
            len(old_column_names), len(old_column_names)+1]
    else:
        new_column_names = list(old_column_names)

    new_column_names.append('phase1')
    new_column_names.append('phase2')
    header = _headerops._update_header_entry(
        header, 'columns', ' '.join(new_column_names))

    if (   ('XB1' not in old_column_names) 
        or ('XB2' not in old_column_names)  
        or ('AS1' not in old_column_names)  
        or ('AS2' not in old_column_names)
        or ('XS1' not in old_column_names)  
        or ('XS2' not in old_column_names)
        ):
        raise ValueError(
            'The input pairs file must be parsed with the flag --add-columns XB,AS,XS --min-mapq 0')

    COL_XB1 = old_column_names.index('XB1')
    COL_XB2 = old_column_names.index('XB2')
    COL_AS1 = old_column_names.index('AS1')
    COL_AS2 = old_column_names.index('AS2')
    COL_XS1 = old_column_names.index('XS1')
    COL_XS2 = old_column_names.index('XS2')

    outstream.writelines((l+'\n' for l in header))


    def get_chrom_phase(chrom, phase_suffixes):
        if chrom.endswith(phase_suffixes[0]):
            return '0', chrom[:-len(phase_suffixes[0])]
        elif chrom.endswith(phase_suffixes[1]):
            return '1', chrom[:-len(phase_suffixes[1])]
        else:
            return '!', chrom


    def phase_side(chrom, XB, AS, XS, phase_suffixes):
        phase, chrom_base = get_chrom_phase(chrom, phase_suffixes)
        XBs = [i for i in XB.split(';') if len(i)>0]

        if AS > XS:
            return phase, chrom_base

        elif len(XBs) >= 1:
            if len(XBs) >= 2:
                alt2_chrom, alt2_pos, alt2_CIGAR, alt2_NM, alt2_AS = XBs[1].split(',')
                if alt2_AS == XS == AS:
                    return '!', '!'

            alt_chrom, alt_pos, alt_CIGAR, alt_NM, alt_AS = XBs[0].split(',')
            alt_phase, alt_chrom_base = get_chrom_phase(alt_chrom, phase_suffixes)

            alt_is_homologue = (
                (chrom_base == alt_chrom_base)
                and
                (
                    ((phase=='0') and (alt_phase=='1'))
                    or 
                    ((phase=='1') and (alt_phase=='0'))
                )
            )
            
            if alt_is_homologue:
                return '.', chrom_base

        return '!', '!'


    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        cols.append('!')
        cols.append('!')
        pair_type = cols[_pairsam_format.COL_PTYPE]

        if cols[_pairsam_format.COL_C1] != _pairsam_format.UNMAPPED_CHROM:

            phase1, chrom_base1 = phase_side(
                cols[_pairsam_format.COL_C1],
                cols[COL_XB1], 
                int(cols[COL_AS1]),
                int(cols[COL_XS1]),
                phase_suffixes
                )

            cols[-2] = phase1
            cols[_pairsam_format.COL_C1] = chrom_base1

            if chrom_base1 == '!':
                cols[_pairsam_format.COL_C1] = _pairsam_format.UNMAPPED_CHROM
                cols[_pairsam_format.COL_P1] = str(_pairsam_format.UNMAPPED_POS)
                cols[_pairsam_format.COL_S1] = _pairsam_format.UNMAPPED_STRAND
                pair_type = 'M' + pair_type[1]

        if cols[_pairsam_format.COL_C2] != _pairsam_format.UNMAPPED_CHROM:

            phase2, chrom_base2 = phase_side(
                cols[_pairsam_format.COL_C2],
                cols[COL_XB2], 
                int(cols[COL_AS2]),
                int(cols[COL_XS2]),
                phase_suffixes
                )

            cols[-1] = phase2
            cols[_pairsam_format.COL_C2] = chrom_base2

            if chrom_base2 == '!':
                cols[_pairsam_format.COL_C2] = _pairsam_format.UNMAPPED_CHROM
                cols[_pairsam_format.COL_P2] = str(_pairsam_format.UNMAPPED_POS)
                cols[_pairsam_format.COL_S2] = _pairsam_format.UNMAPPED_STRAND
                pair_type = pair_type[0] + 'M'

        cols[_pairsam_format.COL_PTYPE] = pair_type

        if clean_output:
            cols = [cols[i] for i in new_column_idxs]
           
        outstream.write(_pairsam_format.PAIRSAM_SEP.join(cols))
        outstream.write('\n')

    if instream != sys.stdin:
        instream.close()

    if outstream != sys.stdout:
        outstream.close()

if __name__ == '__main__':
    phase()
