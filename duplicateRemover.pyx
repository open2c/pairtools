import numpy as np 
cimport numpy as np 
import cython
cimport cython 


"""
duplicateRemoveMask is a static method that performs duplicate removal. 

For other applications on much larger datasets you may consider a processive method processiveDuplicateRemoval
which is implemented as a class. 

Note that for both methods data types are fixed: chromosomes are int16, position is int32, and strand is bool / int8, which are basically the same C type "char". 
"""

def duplicateRemoveMask(cython.short [:] c1,
                        cython.short [:] c2,
                        cython.int [:] p1, 
                        cython.int [:] p2, 
                        cython.char [:] s1, 
                        cython.char [:] s2, 
               #uncomment for testing probably since it will do boundary check
               #np.ndarray[np.int16_t, ndim=1]c2,
               #np.ndarray[np.int32_t, ndim=1] p1, 
               #np.ndarray[np.int32_t, ndim=1] p2, 
               #np.ndarray[np.int8_t, ndim=1] s1,
               #np.ndarray[np.int8_t, ndim=1] s2,
               int offset=3, 
                method = "sum"):
    """
    
    This is a method that removes duplicates with allowing for some tolerance in distance. You can use it to filter single-cell data as wel by setting offset to 500bp or even larger. It works as fast as we could make it. 
    
    Parameters: 
    c1 c2: int16 array of chromosome IDs
    p1, p2: int32 arrays of position 
    s1, s2: int8 (or bool) arrays if strand 
    Arrays MUST be ordered by (c1, p1) 
    offset: int
    method: "sum" or "max" use sum of distances, or max of the two
    
    Returns: 
        mask of elements, where 1 denotes reads that needs to be removed 
    
    This methods scans through a list of reads. It then removes duplicates,
    which are defined as reads starting within offset from each other. 
    
    There are two ways of identifying duplicates: 
    "sum": two reads are duplicates if a sum of distance between their starts, and 
           distance between their ends, is more or equal than "offset"
           
    "max": two reads are duplicates if at least one distance is more-or-equal than "offset"
    
    Other methods could be added below by editing the code 
     
    
    """
    cdef int N = len(c1)
    cdef np.ndarray[np.int8_t, ndim=1] rm = np.zeros(N, dtype = np.int8)
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
        if low == N:
            break

        if high == N:
            low += 1
            high = low + 1
            continue

        if rm[low] == 1:
            low += 1
            high = low+1
            continue

        if rm[high] == 1:   # if high already removed, just continue
            high += 1
            continue
        if (c1[high] != c1[low]) or (p1[high] - p1[low] >= offset):   # if we jumped too far, continue
            low += 1
            high = low + 1  # restart high
            continue
        if methodid == 0: 
            extraCondition = max(abs(p1[low] - p1[high]), abs(p2[low] - p2[high])) < offset
        elif methodid == 1: 
            extraCondition = abs(p1[low] - p1[high]) + abs(p2[low] - p2[high]) < offset  #sum of distances <offset
        else:
            raise ValueError("Unknown method id, this should not happen. Check code of this function.")

        if (c2[low] == c2[high]) and (s1[low] == s1[high]) and (s2[low] == s2[high]) and extraCondition:
            rm[high] = 1
            high += 1
            continue
        high += 1

    return rm


cdef class processiveDuplicateRemoval(object):
    cdef cython.short [:] c1
    cdef cython.short [:] c2 
    cdef cython.int [:] p1 
    cdef cython.int [:] p2 
    cdef cython.char [:] s1 
    cdef cython.char [:] s2 
    cdef cython.char [:] rm
    cdef int methodid
    cdef int low
    cdef int high
    cdef int N 
    cdef int offset
    cdef int returnData 
    
    def __init__(self, method, offset, returnData=False):
        if returnData == False:
            self.returnData = 0
        else:
            self.returnData = 1            
        self.N = 0 
        self.c1 = np.zeros(0, np.int16)
        self.c2 = np.zeros(0, np.int16)
        self.p1 = np.zeros(0, np.int32)
        self.p2 = np.zeros(0, np.int32)
        self.s1 = np.zeros(0, np.int8)
        self.s2 = np.zeros(0, np.int8) 
        self.rm = np.zeros(0, np.int8)        
        if method == "max": 
            self.methodid = 0
        elif method == "sum":
            self.methodid = 1 
        else:
            raise ValueError('method should be "sum" or "max"')
        self.offset = int(offset) 
        self.low = 0 
        self.high = 1 

    def addChunk(self, c1, c2, p1, p2, s1, s2):            
        self.c1 = np.concatenate([self.c1, c1])
        self.c2 = np.concatenate([self.c2, c2])
        self.p1 = np.concatenate([self.p1, p1])
        self.p2 = np.concatenate([self.p2, p2])
        self.s1 = np.concatenate([self.s1, s1])
        self.s2 = np.concatenate([self.s2, s2])
        self.rm = np.concatenate([self.rm, np.zeros(len(c1), dtype=np.int8)])
        self.N = self.N + len(c1)        
        return self.run(finish=False)
    
    def shrink(self):
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
        
        
    def run(self, finish=False):
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

            if self.rm[self.high] == 1:   # if high already removed, just continue
                self.high += 1
                continue
            if (self.c1[self.high] != self.c1[self.low]) or (self.p1[self.high] - self.p1[self.low] >= self.offset):   # if we jumped too far, continue
                self.low += 1
                self.high = self.low + 1  # restart high
                continue
            if self.methodid == 0: 
                extraCondition = max(abs(self.p1[self.low] - self.p1[self.high]), abs(self.p2[self.low] - self.p2[self.high])) < self.offset
            elif self.methodid == 1: 
                extraCondition = abs(self.p1[self.low] - self.p1[self.high]) + abs(self.p2[self.low] - self.p2[self.high]) < self.offset  #sum of distances <offset
            else:
                raise ValueError("Unknown method id, this should not happen. Check code of this function.")

            if (self.c2[self.low] == self.c2[self.high]) and (self.s1[self.low] == self.s1[self.high]) and (self.s2[self.low] == self.s2[self.high]) and extraCondition:
                self.rm[self.high] = 1
                self.high += 1
                continue
            self.high += 1

        return self.shrink()
    
    def finish(self):
        return self.run(finish=True)
    
