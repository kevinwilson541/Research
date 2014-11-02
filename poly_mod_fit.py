import numpy
import pycuda.gpuarray as gpuarray
import pycuda.autoinit
import gaussian_mod_elimination as gme

def poly_mod_fit(v, p):
    n = len(v)
    num_mat = numpy.array([[0 for j in range(n+1)] for i in range(n)]).astype(numpy.int32)
    for i in range(n):
        for j in range(n):
            num_mat[i][j] = numpy.int32(i**j % p)
    num_mat[0][0] = numpy.int32(1 % p)
    for i in range(n):
        num_mat[i][n] = numpy.int32(v[i] % p)

    gpu_mat = gpuarray.to_gpu(num_mat)
    gme.gaussian_mod_elimination(gpu_mat, p)
    num_mat = gpu_mat.get()
    coeffs = numpy.array([num_mat[i][n] for i in range(n)])
    return lambda x: numpy.dot(numpy.array([x**i % p for i in range(n)]), coeffs)

f = poly_mod_fit([i**2 % 11 for i in range(5)], 11)
print [f(i) % 11 for i in range(10)]
