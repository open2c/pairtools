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
