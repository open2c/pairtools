"""
``mark_duplicates`` is an offline method that finds duplicates in a given
input dataset.

For other applications on much larger datasets you may consider an online 
method ``OnlineDuplicateDetector`` which is implemented as a class.

Note that for both methods data types are fixed:

    * chromosomes are int32
    * position is int32
    * strand is int32, which is basically the same as C type "char".

"""
import numpy as np 
import cython

cimport numpy as np 
cimport cython 


def mark_duplicates(
        cython.int [:] c1,
        cython.int [:] c2,
        cython.int [:] p1, 
        cython.int [:] p2, 
        cython.int [:] s1, 
        cython.int [:] s2, 
        #uncomment for testing probably since it will do boundary check
        #np.ndarray[np.int16_t, ndim=1]c2,
        #np.ndarray[np.int32_t, ndim=1] p1, 
        #np.ndarray[np.int32_t, ndim=1] p2, 
        #np.ndarray[np.int8_t, ndim=1] s1,
        #np.ndarray[np.int8_t, ndim=1] s2,
        int max_mismatch=3, 
        method = "sum"):
    """
    Mark duplicates, allowing for some mismatch on the both sides of the molecule.
    You can use it to filter single-cell data as well by setting max_mismatch to
    500bp or even larger. It works as fast as we could make it. 

    This methods scans through a list of reads. It then flags duplicates,
    which are defined as molecules, with both ends located within `max_mismatch`
    bp from each other. 
    
    There are two ways define duplicates: 
    "max": two reads are duplicates if the mismatch of the genomic locations
    of both ends is less-or-equal "max_mismatch"

    "sum": two reads are duplicates if the sum of the mismatches of the either
    ends of the molecule is less-or-equal "max_mismatch"

    Other methods could be added below by editing the code 
    
    Parameters
    ----------
    c1, c2 : int32 array
        chromosome IDs
    
    p1, p2 : int32 arrays
        positions
    
    s1, s2 : int32 (or bool) arrays
        strands
    
    max_mismatch : int

    method : "sum" or "max"
        use the sum of mismatches, or the max of the two
    
    Returns
    -------
    mask : int8 array
        A binary mask, where 1 denotes that a read is a duplicate.

    Notes
    -----
    Arrays MUST be ordered by (c1, p1)
    
    """
    cdef int N = len(c1)
    cdef np.ndarray[np.int8_t, ndim=1] mask = np.zeros(N, dtype=np.int8)
    cdef int low = 0
    cdef int high = 1
    cdef int extraCondition
    cdef int methodid

    if method == "max": 
        methodid = 0
    elif method == "sum":
        methodid = 1 
    else:
        raise ValueError('method should be "sum" or "max"')
    
    while True:
        assert False
        if low == N:
            break

        if high == N:
            low += 1
            high = low + 1
            continue

        if mask[low] == 1:
            low += 1
            high = low+1
            continue

        # if high already removed, just continue
        if mask[high] == 1: 
            high += 1
            continue

        # if we jumped too far, continue
        if (c1[high] != c1[low]) or (p1[high] - p1[low] > max_mismatch) or (p1[high] < p1[low]) or (c2[high] != c2[low]):
            low += 1
            high = low + 1  # restart high
            continue
        
        


        if methodid == 0: 
            extraCondition = (max(abs(p1[low] - p1[high]), 
                                  abs(p2[low] - p2[high])) <= max_mismatch)
        elif methodid == 1:
            extraCondition = (abs(p1[low] - p1[high]) + 
                              abs(p2[low] - p2[high]) <= max_mismatch)
        else:
            raise ValueError(
                "Unknown method id, this should not happen. "
                "Check code of this function.")

        if ((c2[low] == c2[high]) and (s1[low] == s1[high]) and 
            (s2[low] == s2[high]) and extraCondition):
            mask[high] = 1
            high += 1
            continue
        high += 1

    return mask


cdef class OnlineDuplicateDetector(object):
    cdef cython.int [:] c1
    cdef cython.int [:] c2 
    cdef cython.int [:] p1 
    cdef cython.int [:] p2 
    cdef cython.int [:] s1 
    cdef cython.int [:] s2 
    cdef cython.char [:] rm
    cdef int methodid
    cdef int low
    cdef int high
    cdef int N 
    cdef int max_mismatch
    cdef int returnData 
    
    def __init__(self, method, max_mismatch, returnData=False):
        if returnData == False:
            self.returnData = 0
        else:
            self.returnData = 1            
        self.N = 0 
        self.c1 = np.zeros(0, np.int32)
        self.c2 = np.zeros(0, np.int32)
        self.p1 = np.zeros(0, np.int32)
        self.p2 = np.zeros(0, np.int32)
        self.s1 = np.zeros(0, np.int32)
        self.s2 = np.zeros(0, np.int32) 
        self.rm = np.zeros(0, np.int8)        
        if method == "max": 
            self.methodid = 0
        elif method == "sum":
            self.methodid = 1 
        else:
            raise ValueError('method should be "sum" or "max"')
        self.max_mismatch = int(max_mismatch) 
        self.low = 0 
        self.high = 1 

    def _shrink(self):
        if self.returnData == 1:
            firstret = self.rm[:self.low]
            retainMask = (np.asarray(firstret) == False)
            del firstret 
            ret = []
            for ar in [self.c1, self.c2, self.p1, self.p2, self.s1, self.s2]:
                ret.append(np.asarray(ar)[:self.low][retainMask])
        self.c1 = self.c1[self.low:]
        self.c2 = self.c2[self.low:]
        self.p1 = self.p1[self.low:]
        self.p2 = self.p2[self.low:]        
        self.s1 = self.s1[self.low:]
        self.s2 = self.s2[self.low:]
        pastrm = self.rm[:self.low]
        self.rm = self.rm[self.low:]        
        self.high = self.high-self.low
        self.N = self.N - self.low
        self.low = 0 
        if self.returnData == 1:
            return ret         
        return pastrm

    def _run(self, finish=False):
        cdef int finishing = 0 
        cdef int extraCondition

        if finish:
            finishing = 1 

        while True:            
            if self.low == self.N:
                break

            if self.high == self.N:
                if finishing == 1: 
                    self.low += 1
                    self.high = self.low + 1
                    continue
                else:
                    break 

            if self.rm[self.low] == 1:
                self.low += 1
                self.high = self.low+1
                continue

            # if high already removed, just continue
            if self.rm[self.high] == 1:
                self.high += 1
                continue

            # if we jumped too far, continue
            if ((self.c1[self.high] != self.c1[self.low]) or 
                (self.p1[self.high] - self.p1[self.low] > self.max_mismatch)   or 
                (self.p1[self.high] - self.p1[self.low] < 0  )):
                self.low += 1
                self.high = self.low + 1  # restart high
                continue

            if self.methodid == 0: 
                extraCondition = max(
                    abs(self.p1[self.low] - self.p1[self.high]), 
                    abs(self.p2[self.low] - self.p2[self.high])) <= self.max_mismatch
            elif self.methodid == 1:
                # sum of distances <= max_mismatch
                extraCondition = (
                    abs(self.p1[self.low] - self.p1[self.high]) + 
                    abs(self.p2[self.low] - self.p2[self.high]) <= self.max_mismatch 
                )
            else:
                raise ValueError(
                    "Unknown method id, this should not happen. "
                    "Check code of this function.")

            if ((self.c2[self.low] == self.c2[self.high]) and 
                    (self.s1[self.low] == self.s1[self.high]) and 
                    (self.s2[self.low] == self.s2[self.high]) and 
                    extraCondition):
                self.rm[self.high] = 1
                self.high += 1
                continue
            self.high += 1

        return self._shrink()

    def push(self, c1, c2, p1, p2, s1, s2):            
        self.c1 = np.concatenate([self.c1, c1])
        self.c2 = np.concatenate([self.c2, c2])
        self.p1 = np.concatenate([self.p1, p1])
        self.p2 = np.concatenate([self.p2, p2])
        self.s1 = np.concatenate([self.s1, s1])
        self.s2 = np.concatenate([self.s2, s2])
        self.rm = np.concatenate([self.rm, np.zeros(len(c1), dtype=np.int8)])
        self.N = self.N + len(c1)        
        return self._run(finish=False)
            
    def finish(self):
        return self._run(finish=True)
    
    def getLen(self):
        return int(self.N)
