import pycuda.driver as driver
import pycuda.gpuarray as gpuarray
import numpy
import pycuda.autoinit
from pycuda.compiler import SourceModule
import time

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
        normalize(m,numpy.int32(i),numpy.int32(m.shape[1]), block=(cols,1,1))
        l = range(rows)
        l.remove(i)
        for j in l:
            rref(m, numpy.int32(i), numpy.int32(j), numpy.int32(m.shape[1]), block=(cols,1,1))

def poly_fit(v):
    n = len(v)
    num_mat = numpy.array([[0 for j in range(n+1)] for i in range(n)]).astype(numpy.float32)
    for i in range(n):
        for j in range(n):
            num_mat[i][j] = numpy.float32(i**j)
    num_mat[0][0] = numpy.float32(1)
    for i in range(n):
        num_mat[i][n] = numpy.float32(v[i])

    gpu_mat = gpuarray.to_gpu(num_mat)
    gaussian_elimination(gpu_mat)
    num_mat = gpu_mat.get()
    coeffs = numpy.array([num_mat[i][n] for i in range(n)])
    return lambda x: numpy.dot(numpy.array([x**i for i in range(n)]), coeffs)


f = poly_fit([i for i in range(10)])
print [f(i) for i in range(20)]

