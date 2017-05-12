from collections import OrderedDict
import sys
import copy
import itertools

from . import __version__, _pairsam_format


PAIRS_FORMAT_VERSION = '1.0.0'


def get_header(instream, comment_char='#'):
    '''Returns a header from the stream and an iterator for the remaining
    lines.

    Parameters
    ----------
    instream : a file object
        An input stream.

    comment_char : str
        The character prepended to header lines (use '@' when parsing sams, 
        '#' when parsing pairsams).

    Returns
    -------
    header : list of str
        The header lines, stripped of terminal spaces and newline characters.

    line_iterator : an iterator of str
        An iterator over the residual lines of the input.
    
    '''
    header = []
    if not comment_char:
        raise Exception('Please, provide a comment char!')
    line = None
    for line in instream:
        if line.startswith(comment_char):
            header.append(line.strip())
        else:
            break

    if line:
        return header, itertools.chain([line], instream)
    else:
        return header, instream


def extract_fields(header, field_name, save_rest=False):
    '''
    Extract the specified fields from the pairs header and returns
    a list of corresponding values, even if a single field was found.
    Additionally, can return the list of intact non-matching entries.
    '''

    fields = []
    rest = []
    for l in header:
        if l.lstrip('#').startswith(field_name+':'):
            fields.append(l.split(':',1)[1].strip())
        elif save_rest:
            rest.append(l)

    if save_rest:
        return fields, rest
    else:
        return fields

def extract_column_names(header):
    '''
    Extract column names from header lines.
    '''
    columns = extract_fields(header, 'columns')

    if len(columns) != 0:
        return columns[0].split(' ')
    else:
        return []


def get_chromsizes_from_sam_header(samheader):
    SQs = [l.split('\t') for l in samheader if l.startswith('@SQ')]
    chromsizes = [(sq[1][3:], int(sq[2][3:])) for sq in SQs]
    return OrderedDict(chromsizes)


def get_chrom_order(chroms_file, sam_chroms):
    """
    Produce an "enumeration" of chromosomes based on the list
    of chromosomes

    """
    chrom_enum = OrderedDict()
    i = 1
    with open(chroms_file, 'rt') as f:
        for line in f:
            chrom = line.split('\t')[0].strip()
            if chrom:
                chrom_enum[chrom] = i
                i += 1

    remaining = sorted(chrom for chrom in sam_chroms
                        if chrom not in chrom_enum.keys())
    for chrom in remaining:
        chrom_enum[chrom] = i
        i += 1

    return chrom_enum


def make_standard_pairsheader(
        assembly=None,
        chromsizes=None,
        columns=_pairsam_format.COLUMNS):
    header = []
    header.append(
        '## pairs format v{}'.format(PAIRS_FORMAT_VERSION))
    header.append('#shape: upper triangle')

    header.append('#genome_assembly: {}'.format(
        assembly if assembly is not None else 'unknown'))

    if chromsizes is not None:
        try:
            chromsizes = chromsizes.items()
        except AttributeError:
            pass
        for chrom, length in chromsizes:
            header.append('#chromsize: {} {}'.format(chrom, length))

    header.append('#columns: '+ ' '.join(columns))

    return header


def insert_samheader(header, samheader):
    new_header = [l for l in header if not l.startswith('#columns')]
    if samheader:
        new_header += ['#samheader: '+l for l in samheader]
    new_header += [l for l in header if l.startswith('#columns')]
    return new_header


def mark_header_as_sorted(header):
    header = copy.deepcopy(header)
    if not any([l.startswith('#sorted') for l in header]):
        if header[0].startswith('##'):
            header.insert(1, '#sorted: chr1-chr2-pos1-pos2')
        else:
            header.insert(0, '#sorted: chr1-chr2-pos1-pos2')
    for i in range(len(header)):
        if header[i].startswith('#chromosomes'):
            chroms = header[i][12:].strip().split(' ')
            header[i] = '#chromosomes: {}'.format(' '.join(sorted(chroms)))
    return header


def append_new_pg(header, ID='', PN='', VN=None, CL=None, force=False):
    header = copy.deepcopy(header)
    samheader, other_header = extract_fields(header, 'samheader', save_rest=True)
    new_samheader = _add_pg_to_samheader(samheader, ID, PN, VN, CL, force)
    new_header = insert_samheader(other_header, new_samheader)
    return new_header

def _update_header_entry(header, field, new_value):
    header = copy.deepcopy(header)
    found = False
    newline = '#{}: {}'.format(field, new_value)
    for i in range(len(header)):
        if header[i].startswith('#'+field):
            header[i] = newline
            found = True
    if not found:
        if header[-1].startswith('#columns'):
            header.insert(-1, newline)
        else:
            header.append(newline)
    return header


def _add_pg_to_samheader(samheader, ID='', PN='', VN=None, CL=None, force=False):
    '''Append a @PG record to an existing sam header. If the header comes
    from a merged file and thus has multiple chains of @PG, append the
    provided PG to all of the chains, adding the numerical suffix of the 
    branch to the ID.

    Parameters
    ----------

    header : list of str
    ID, PN, VN, CL : std
        The keys of a new @PG record. If absent, VN is the version of
        pairsamtools and CL is taken from sys.argv.
    force : bool
        If True, ignore the inconsistencies among @PG records of the existing 
        header.

    Returns
    -------
    new_header : list of str
        A list of new headers lines, stripped of newline characters.


    '''
    if VN is None:
        VN = __version__
    if CL is None:
        CL = ' '.join(sys.argv)

    pre_pg_header = [line.strip() for line in samheader
                     if line.startswith('@HD')
                         or line.startswith('@SQ')
                         or line.startswith('@RG')
                     ]

    post_pg_header = [line.strip() for line in samheader
                     if not line.startswith('@HD')
                         and (not line.startswith('@SQ'))
                         and (not line.startswith('@RG'))
                         and (not line.startswith('@PG'))
                     ]
        
    pg_chains = _parse_pg_chains(samheader, force=force)

    for i,br in enumerate(pg_chains):
        new_pg = {'ID':ID, 'PN':PN, 'VN':VN, 'CL':CL}
        new_pg['PP'] = br[-1]['ID']
        if len(pg_chains) > 1:
            new_pg['ID'] = new_pg['ID'] + '-' + str(i+1)
        new_pg['raw'] = _format_pg(**new_pg)
        br.append(new_pg)

    new_header = (
        pre_pg_header 
        + [pg['raw'] for br in pg_chains for pg in br]
        + post_pg_header)

    return new_header
        

def _format_pg(**kwargs):
    out = (
        ['@PG']
        + ['{}:{}'.format(field, kwargs[field])
           for field in ['ID', 'PN', 'CL', 'PP', 'DS', 'VN'] 
           if field in kwargs])
    return '\t'.join(out)


def _parse_pg_chains(header, force=False):
    pg_chains = []
    parsed_pgs = [
        dict(
            [field.split(':', maxsplit=1)    
                for field in l.strip().split('\t')[1:]]
            + [('raw', l.strip())]
        )
        for l in header
        if l.startswith('@PG')
        ]

    while True:
        if len(parsed_pgs) == 0:
            break
        
        for i in range(len(parsed_pgs)):
            pg = parsed_pgs[i]
            if 'PP' not in pg:
                pg_chains.append([pg])
                parsed_pgs.pop(i)
                break
            else:
                matching_chains = [
                    branch for branch in pg_chains
                    if branch[-1]['ID'] == pg['PP']
                ]
                if len(matching_chains) > 1:
                    if force:
                        matching_chains[0].append(pg)
                        parsed_pgs.pop(i)
                        break

                    else:
                        raise Exception(
                            'Multiple @PG records with the IDs identical to the PP field of another record:\n'
                            + '\n'.join([br[-1]['raw'] for br in matching_chains])
                            + '\nvs\n'
                            + pg['raw']
                            )


                if len(matching_chains) == 1:
                    matching_chains[0].append(pg)
                    parsed_pgs.pop(i)
                    break
                
            if force:
                pg_chains.append([pg])
                parsed_pgs.pop(i)
                break
            else:
                raise Exception(
                    'Cannot find the parental @PG record for the @PG records:\n'
                    + '\n'.join([pg['raw'] for pg in parsed_pgs])
                    )

    return pg_chains

def _merge_samheaders(samheaders, force=False):
    # first, append an HD line if it is present in any files
    # if different lines are present, raise an error
    HDs = set.union(*[set(line for line in samheader if line.startswith('@HD'))
                      for samheader in samheaders])
    if len(HDs) > 1 and not force:
        raise Exception('More than one unique @HD line is found in samheaders!')
    HDs = [list(HDs)[0]] if HDs else []

    # second, confirm that all files had the same SQ lines
    # add SQs from the first file, keeping its order
    SQs = [set(line for line in samheader if line.startswith('@SQ'))
           for samheader in samheaders]
    common_SQs = set.intersection(*SQs)

    SQs_same = all([len(samheader) == len(common_SQs) 
                    for samheader in SQs])

    if not SQs_same and not(force):
        raise Exception('The SQ (sequence) lines of the sam headers are not identical')
    SQs = [line for line in samheaders[0] if line.startswith('@SQ')]

    # third, append _all_ PG chains, adding a unique index according to the 
    # provided merging order
    PGs = []
    for i, samheader in enumerate(samheaders):
        for line in samheader:
            if line.startswith('@PG'):
                split_line = line.split('\t')
                for j in range(len(split_line)):
                    if (split_line[j].startswith('ID:') 
                    or split_line[j].startswith('PP:')):
                        split_line[j] = split_line[j] + '-' + str(i+1)
                PGs.append('\t'.join(split_line))

    # finally, add all residual unique lines
    rest = sum([
        list(set(line for line in samheader 
                if  (not line.startswith('@HD'))
                and (not line.startswith('@SQ'))
                and (not line.startswith('@PG'))
            ))
        for samheader in samheaders], 
        [])

    new_header = []
    new_header += HDs
    new_header += SQs 
    new_header += PGs
    new_header += rest

    return new_header


def _merge_pairheaders(pairheaders, force=False):
    new_header = []

    # first, add all keys that are expected to be the same among all headers
    keys_expected_identical = [
        '## pairs format',
        '#sorted:',
        '#shape:',
        '#genome_assembly:',
        '#columns:']

    for k in keys_expected_identical:
        lines = [[l for l in header if l.startswith(k)]
                 for header in pairheaders]
        same = all([l==lines[0] for l in lines])
        if not (same or force):
            raise Exception(
                'The following header entries must be the same '
                'the merged files: {}'.format(k))

        new_header += lines[0]

    # second, add the chromosome field.
    # the chromosomes are expected to be the same, but if merge is forced,
    # make a sorted list of unique chromosomes
    chromosome_headers = [[l for l in header if l.startswith('#chromosomes:')]
                          for header in pairheaders]
    if not (all(len(c)==0 for c in chromosome_headers)
            or all(len(c)==1 for c in chromosome_headers)
            or force):
        raise Exception(
                'All headers of the merged files must have the same number '
                '#chromosome entries (either 0 or 1)')

    if len(chromosome_headers[0]) == 1:
        chrom_lists = [set(ch[0].split(':',1)[1].strip().split(' ')) 
                       for ch in chromosome_headers]

        same = all([c==chrom_lists[0] for c in chrom_lists])
        if not(same or force):
            raise Exception('All headers must list the same chromosomes!')

        chroms = sorted(set.union(*chrom_lists))
        chrom_line = '#chromosomes: {}'.format(' '.join(sorted(chroms)))

        if new_header[-1].startswith('#columns'):
            new_header.insert(-1, chrom_line)
        else:
            new_header.append(chrom_line)

    # finally, add a sorted list of other unique fields
    other_lines = sorted(set(
        l for h in pairheaders for l in h
        if not any(l.startswith(k) 
                   for k in keys_expected_identical + ['#chromosomes'])))

    if new_header[-1].startswith('#columns'):
        new_header = new_header[:-1] + other_lines + [new_header[-1]]
    else:
        new_header = new_header + other_lines
    
    return new_header

                    
def merge_headers(headers, force=False):

    samheaders, pairheaders = zip(*[extract_fields(h, 'samheader', save_rest=True) for h in headers])
    # HD headers contain information that becomes invalid after processing
    # with distiller. Do not print into the output.
    new_pairheader = _merge_pairheaders(pairheaders, force=False)
    new_samheader = _merge_samheaders(samheaders, force=force)

    new_header = insert_samheader(new_pairheader, new_samheader)

    return new_header


#def _guess_genome_assembly(samheader):
#    PG = [l for l in samheader if l.startswith('@PG') and '\tID:bwa' in l][0]
#    CL = [field for field in PG.split('\t') if field.startswith('CL:')]
#
#    return ga

