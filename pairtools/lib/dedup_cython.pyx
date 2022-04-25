"""
Legacy code:
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


### Online deduplicator used in pairtools.dedup Cython:
cdef class OnlineDuplicateDetector(object):
    cdef cython.int [:] c1
    cdef cython.int [:] c2
    cdef cython.int [:] p1
    cdef cython.int [:] p2
    cdef cython.int [:] s1
    cdef cython.int [:] s2
    cdef cython.char [:] rm
    cdef cython.int [:] parent_idxs
    cdef int methodid
    cdef int low
    cdef int high
    cdef int N
    cdef int max_mismatch
    cdef int returnData
    cdef int keep_parent_id

    def __init__(self, method, max_mismatch, returnData=False, keep_parent_id=False):
        if returnData == False:
            self.returnData = 0
        else:
            self.returnData = 1
        if keep_parent_id == False:
            self.keep_parent_id = 0
        else:
            self.keep_parent_id = 1
            self.parent_idxs = np.zeros(0, np.int32)

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
        if self.returnData == 1:
            self.low = 0
            return ret
        if self.keep_parent_id == 1: # Return parent readIDs alongside with duplicates mask:
            pastidx = self.parent_idxs[:self.low]
            self.low = 0
            return pastrm, pastidx
        self.low = 0
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
                if self.keep_parent_id == 1:
                    self.parent_idxs[self.high] = self.low
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
        if self.keep_parent_id == 1:
            self.parent_idxs = np.concatenate([self.parent_idxs, np.zeros(len(c1), dtype=np.int32)])
        self.N = self.N + len(c1)
        return self._run(finish=False)

    def finish(self):
        return self._run(finish=True)

    def getLen(self):
        return int(self.N)