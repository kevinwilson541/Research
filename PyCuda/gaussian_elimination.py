import pycuda.driver as driver
import pycuda.gpuarray as gpuarray
import numpy
import pycuda.autoinit
from pycuda.compiler import SourceModule
import time

swap_mod = SourceModule("""
    __device__ int find_max(float * v, int i, int col_size, int row_size) {
        float max = v[col_size*i + i];
        int value = i;
        int j;
        for (j = i; j < row_size; ++j) {
            if (v[col_size*j + i] > max) {
                max = v[col_size*j + i];
                value = j;
            }
        }
        return value;
    }

    __global__ void swap(float * v, int i, int col_size, int row_size) {
        int j = find_max(v, i, col_size, row_size);
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        float tmp = v[col_size*i + col];
        v[col_size*i + col] = v[col_size*j + col];
        v[col_size*j + col] = tmp;
    }
    """)
swap_mod = swap_mod.get_function("swap")

normalize_mod = SourceModule("""
    __global__ void normalize(float * v, int i, int col_size) {
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        if (v[col_size*i + i] != 0) {
            v[col_size*i + col] /= v[col_size*i+i];
        }
    }
    """)
normalize = normalize_mod.get_function("normalize")

rref_mod = SourceModule("""
    __global__ void reduce(float * v, int i, int j, int col_size) {
        int col = blockIdx.x*blockDim.x + threadIdx.x;
        float factor = v[j*col_size + i];
        v[j*col_size + col] -= (float)factor*v[i*col_size + col];
    }
    """)
rref = rref_mod.get_function("reduce")

def gaussian_elimination(m):
    rows = m.shape[0]
    cols = m.shape[1]
    for i in range(rows):
        swap_mod(m, numpy.int32(i), numpy.int32(m.shape[1]), numpy.int32(m.shape[0]), block=(cols,1,1))
        normalize(m,numpy.int32(i),numpy.int32(m.shape[1]), block=(cols,1,1))
        l = range(rows)
        l.remove(i)
        for j in l:
            rref(m, numpy.int32(i), numpy.int32(j), numpy.int32(m.shape[1]), block=(cols,1,1))

#a = gpuarray.to_gpu(numpy.array([[i for j in range(5)] for i in range(5)]).astype(numpy.float32))
#gaussian_elimination(a)
#print a.get()
