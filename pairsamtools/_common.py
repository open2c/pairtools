import pipes
import copy
import collections
import itertools

COL_READID = 0
COL_C1 = 1
COL_C2 = 2
COL_P1 = 3
COL_P2 = 4
COL_S1 = 5
COL_S2 = 6
COL_PTYPE = 7
COL_SAM1 = 8
COL_SAM2 = 9

SAM_ENTRY_SEP = '\tNEXT_SAM\t'

def open_sam_or_bam(path, mode):
    '''Opens a file as a bam file is `path` ends with .bam, otherwise 
    opens it as a sam.
    '''
    if mode not in ['r','w']:
        raise Exception("mode can be either 'r' or 'w'")
    if path.endswith('.bam'):
        if mode =='w': 
            t = pipes.Template()
            t.append('samtools view -bS', '--')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append('samtools view -h', '--')
            f = t.open(path, 'r')
        else:
            raise Exception("Unknown mode : {}".format(mode))
        return f
    else:
        return open(path, mode)


def open_bgzip(path, mode):
    '''Opens a file as a bgzip file is `path` ends with .bam, otherwise 
    opens it as a text.
    '''
    if mode not in ['r','w']:
        raise Exception("mode can be either 'r' or 'w'")
    if path.endswith('.gz'):
        if mode =='w': 
            t = pipes.Template()
            t.append('bgzip -c', '--')
            f = t.open(path, 'w')
        elif mode =='r': 
            t = pipes.Template()
            t.append('zcat', '--')
            f = t.open(path, 'r')
        else:
            raise Exception("Unknown mode : {}".format(mode))
        return f
    else:
        return open(path, mode)


def get_header(instream, comment_char='#'):
    '''Returns a header from the stream and an iterator for the remaining
    lines.

    Parameters
    ----------
    instream : a file object
        An input stream.

    comment_char : str
        The character prepended to header lines (use '' when parsing sams, 
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
        comment_char = '@'
    for line in instream:
        if line.startswith(comment_char):
            header.append(line.strip())
        else:
            break
        
    return header, itertools.chain([line], instream)


def append_pg_to_sam_header(header, pg_dict, comment_char='#', force=False):
    '''Append a @PG record to an existing sam header. If the header comes
    from a merged file and thus has multiple branches of @PG, append the
    provided PG to all of the branches, adding the numerical suffix of the 
    branch to the ID.

    Parameters
    ----------

    header : list of str
    pg_dict : dict
        A dictionary describing the new @PG record. Must contain the keys "ID" 
        and "PN". The keys "CL", "DS", "VN" are processed if present.
    comment_char : str
        If provided, assume that all lines of the header must and do start
        with an extra `comment_char` character (pairsam-specific). 
    force : bool
        If True, ignore the inconsistencies among @PG records of the existing 
        header.

    Returns
    -------
    new_header : list of str
        A list of new headers lines, stripped of newline characters.


    '''
    pre_pg_header = [line.strip() for line in header
                     if line.startswith(comment_char + '@HD')
                         or line.startswith(comment_char + '@SQ')
                         or line.startswith(comment_char + '@RG')
                     ]

    post_pg_header = [line.strip() for line in header
                     if not line.startswith(comment_char + '@HD')
                         and (not line.startswith(comment_char + '@SQ'))
                         and (not line.startswith(comment_char + '@RG'))
                         and (not line.startswith(comment_char + '@PG'))
                     ]
        
    pg_branches = _parse_pg_branches(
            header, comment_char=comment_char, force=force)

    for i,br in enumerate(pg_branches):
        new_pg = copy.deepcopy(pg_dict)
        new_pg['PP'] = br[-1]['ID']
        if len(pg_branches) > 1:
            new_pg['ID'] = new_pg['ID'] + '-' + str(i+1)
        new_pg['raw'] = _format_pg(new_pg, comment_char=comment_char)
        br.append(new_pg)

    new_header = (
        pre_pg_header 
        + [pg['raw'] for br in pg_branches for pg in br]
        + post_pg_header)

    new_header = [l+'\n' for l in new_header]

    return new_header
        

def _format_pg(pg_dict, comment_char='#'):
    out = (
        [comment_char+'@PG'] 
        + ['{}:{}'.format(field,pg_dict[field])
           for field in ['ID', 'PN', 'CL', 'PP', 'DS', 'VN'] 
           if field in pg_dict])
    return '\t'.join(out)


def _parse_pg_branches(header, comment_char='#', force=False):
    pg_branches = []
    parsed_pgs = [
        dict(
            [field.split(':', maxsplit=1)    
                for field in l.strip().split('\t')[1:]]
            + [('raw', l.strip())]
        )
        for l in header
        if l.startswith(comment_char+'@PG')
        ]

    while True:
        if len(parsed_pgs) == 0:
            break
        
        for i in range(len(parsed_pgs)):
            pg = parsed_pgs[i]
            if 'PP' not in pg:
                pg_branches.append([pg])
                parsed_pgs.pop(i)
                break
            else:
                matching_branches = [
                    branch for branch in pg_branches
                    if branch[-1]['ID'] == pg['PP']
                ]
                if len(matching_branches) > 1:
                    if force:
                        matching_branches[0].append(pg)
                        parsed_pgs.pop(i)
                        break

                    else:
                        raise Exception(
                            'Multiple @PG records with the IDs identical to the PP field of another record:\n'
                            + '\n'.join([br[-1]['raw'] for br in matching_branches])
                            + '\nvs\n'
                            + pg['raw']
                            )


                if len(matching_branches) == 1:
                    matching_branches[0].append(pg)
                    parsed_pgs.pop(i)
                    break
                
            if force:
                pg_branches.append([pg])
                parsed_pgs.pop(i)
                break
            else:
                raise Exception(
                    'Cannot find the parental @PG record for the @PG records:\n'
                    + '\n'.join([pg['raw'] for pg in parsed_pgs])
                    )

    return pg_branches



    


