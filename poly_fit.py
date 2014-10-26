import numpy
import pycuda.gpuarray as gpuarray
import pycuda.autoinit
import gaussian_elimination as ge

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
    ge.gaussian_elimination(gpu_mat)
    num_mat = gpu_mat.get()
    coeffs = numpy.array([num_mat[i][n] for i in range(n)])
    return lambda x: numpy.dot(numpy.array([x**i for i in range(n)]), coeffs)


f = poly_fit([i for i in range(10)])
print [f(i) for i in range(20)]

