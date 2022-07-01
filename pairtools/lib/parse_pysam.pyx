from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libcalignedsegment cimport AlignedSegment, AlignmentHeader
from pysam.libchtslib cimport *
from pysam.libcutils cimport array_to_qualitystring

cdef class AlignmentFilePairtoolized(AlignmentFile):
    """ Modified class that loads each entry as pairtoolozed alignment. """

    def __next__(self):
        cdef int ret = self.cnext()
        if (ret >= 0):
            # Redefine the constructed object:
            return makeAlignedSegmentPairtoolized(self.b, self.header)
        elif ret == -2:
            raise IOError('truncated file')
        else:
            raise StopIteration

cdef AlignedSegmentPairtoolized makeAlignedSegmentPairtoolized(bam1_t *src,
                                       AlignmentHeader header):
    '''return an AlignedSegmentPairtoolized object constructed from `src`'''
    # note that the following does not call __init__
    # Redefine the constructed object:
    cdef AlignedSegmentPairtoolized dest = AlignedSegmentPairtoolized.__new__(AlignedSegmentPairtoolized)
    dest._delegate = bam_dup1(src)
    dest.header = header
    return dest

cdef class AlignedSegmentPairtoolized(AlignedSegment):
    """ In the pairtoolized class we inherit everything and
        add some useful properties and functions on top of that.
    """

    def is_unique(self, min_mapq):
        """true if read is unique mapping (by mapq)"""
        return self.mapq >= min_mapq

    property is_linear:
        """true if read is linear (SA is present in tages)"""
        def __get__(self):

            if self.has_tag('SA'):
                return False
            # for tag in self.tags:
            #     if 'SA'==tag[0]:
            #         return False
            return True

    property cigar_dict:
        """Parsed CIGAR as dictionary with interpretable fields"""

        def __get__(self):
            """Parse cigar tuples reported as cigartuples of pysam read entry.
            Reports alignment span, clipped nucleotides and more.
            See https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
            """
            matched_bp = 0
            algn_ref_span = 0
            algn_read_span = 0
            read_len = 0
            clip5_ref = 0
            clip3_ref = 0

            cigarstring = self.cigarstring
            cigartuples = self.cigartuples
            if cigartuples is not None:
                for operation, length in cigartuples:
                    if operation == 0:  # M, match
                        matched_bp += length
                        algn_ref_span += length
                        algn_read_span += length
                        read_len += length
                    elif operation == 1:  # I, insertion
                        algn_read_span += length
                        read_len += length
                    elif operation == 2:  # D, deletion
                        algn_ref_span += length
                    elif (
                            operation == 4 or operation == 5
                    ):  # S and H, soft clip and hard clip, respectively
                        read_len += length
                        if matched_bp == 0:
                            clip5_ref = length
                        else:
                            clip3_ref = length

            return {
                "clip5_ref": clip5_ref,
                "clip3_ref": clip3_ref,
                "cigar": cigarstring,
                "algn_ref_span": algn_ref_span,
                "algn_read_span": algn_read_span,
                "read_len": read_len,
                "matched_bp": matched_bp
            }


from cpython cimport array
import cython
cimport cython

cpdef list get_mismatches_c(str seq, array.array quals, list aligned_pairs):
    '''
    This function takes a SAM alignment and, for every mismatch between the read and reference sequences,
    returns a tuple (the reference bp, the read bp, PHRED quality of the bp, reference position, read position).
    
    Reference: https://github.com/gerlichlab/scshic_pipeline/blob/master/bin/seq_mismatches.pyx 
    '''

    cdef cython.int read_pos, ref_pos
    cdef str orig_bp, orig_bp_upper
    cdef list mismatches = []

    for read_pos, ref_pos, orig_bp in aligned_pairs:
        orig_bp_upper = orig_bp.upper()
        if (seq[read_pos] != orig_bp_upper):
            mismatches.append(
                (orig_bp_upper,
                 seq[read_pos],
                 quals[read_pos],
                 ref_pos,
                 read_pos)
            )

    return mismatches