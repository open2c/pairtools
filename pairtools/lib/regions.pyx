""" Moved from pairlib, library for fast regions assignment """
from cython.operator cimport dereference, postincrement, postdecrement

from cpython cimport array
import cython

from libcpp.map cimport map
from libcpp.algorithm cimport lower_bound, upper_bound
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cpdef np.ndarray assign_regs_c(np.ndarray chroms, np.ndarray pos, dict reg_dict):
    assert len(chroms) == len(pos)

    cdef int n = len(chroms)
    cdef np.ndarray[np.int64_t, ndim=2] result = -1 * np.ones((n, 3), dtype=np.int64)
    cdef map[string, vector[int]] reg_map = reg_dict

    cdef map[string, vector[int]].iterator reg_map_it = reg_map.begin()
    cdef map[string, vector[int]].iterator reg_map_end = reg_map.end()

    cdef vector[int].iterator lo_b, up_b
    cdef int position, reg_boundary_idx

    # this can be parallelized with prange
    for i in range(n):
        reg_map_it = reg_map.find(chroms[i])
        if reg_map_it != reg_map_end:
            position = pos[i]
            up_b = upper_bound(
                dereference(reg_map_it).second.begin(),
                dereference(reg_map_it).second.end(),
                position)
            reg_boundary_idx = up_b - dereference(reg_map_it).second.begin()

            if reg_boundary_idx % 2 == 1:
                lo_b = up_b
                postdecrement(lo_b)
                result[i, 0] = (reg_boundary_idx - 1) // 2
                result[i, 1] = dereference(lo_b)
                result[i, 2] = dereference(up_b)

    return result