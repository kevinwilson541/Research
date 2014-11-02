import pycuda.driver as driver
import pycuda.gpuarray as gpuarray
import numpy
import pycuda.autoinit
from pycuda.compiler import SourceModule
import time

swap_mod = SourceModule("""
    __device__ int find_max(int * v, int i, int col_size, int row_size) {
        int max = v[col_size*i + i], value = i;
        int j;
        for (j = i; j < row_size; ++j) {
            if (v[col_size*j + i] > max) {
                max = v[col_size*j + i];
                value = j;
            }
        }
        return value;
    }

    __global__ void swap(int * v, int i, int col_size, int row_size) {
        int j = find_max(v, i, col_size, row_size);
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        int tmp = v[col_size*i + col];
        v[col_size*i + col] = v[col_size*j + col];
        v[col_size*j + col] = tmp;
    }
    """)
swap_mod = swap_mod.get_function("swap")

normalize_mod = SourceModule("""
    __device__ int inverse(int a, int n) {
        int t = 0, r = n, newt = 1, newr = a;
        while (newr != 0) {
            int quotient = r / newr;
            int tmp = newt;
            newt = t - quotient*newt;
            t = tmp;
            int tmp2 = newr;
            newr = r - quotient*newr;
            r = tmp2;
        }
        if (t < 0) t += n;
        return t;
    }
    
    __global__ void normalize(int * v, int i, int col_size, int p) {
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        if (v[col_size*i + i] != 0) {
            int inverse_num = inverse(v[col_size*i + i], p);
            v[col_size*i + col] *= inverse_num;
            if (v[col_size*i + col] >= p) v[col_size*i + col] -= p;
        }
    }
    """)
normalize= normalize_mod.get_function("normalize")

rref_mod = SourceModule("""
    __global__ void reduce(int * v, int i, int j, int col_size, int p) {
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        int factor = v[j*col_size + i];
        v[j*col_size + col] = (v[j*col_size + col] - factor*v[i*col_size + col])%p;
        if (v[j*col_size + col] < 0) v[j*col_size + col] += p;
    }
    """)
rref = rref_mod.get_function("reduce")

def gaussian_mod_elimination(m, p):
    rows = m.shape[0]
    cols = m.shape[1]
    for i in range(rows):
        swap_mod(m, numpy.int32(i), numpy.int32(m.shape[1]), numpy.int32(m.shape[0]), block=(cols,1,1))
        normalize(m, numpy.int32(i), numpy.int32(m.shape[1]), numpy.int32(p), block=(cols,1,1))
        l = range(rows)
        l.remove(i)
        for j in l:
            rref(m, numpy.int32(i), numpy.int32(j), numpy.int32(m.shape[1]), numpy.int32(p), block=(cols,1,1))

