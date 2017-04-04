from __init__ import *
import numba, time

def matrmult(a, b, c):
    matmul(a, b, c)


@numba.jit('void(double[:,:], double[:, :], double[:, :])', nopython=True, nogil=True)
def matrmultjitt(a, b, c):
    for i in range(c.shape[0]):
    	for j in range(c.shape[1]):
    		for k in range(a.shape[1]):
        		c[i, j] += a[i, k] * b[k, i]


@numba.cuda.jit
def matrmultgpuu(matrfrst, matrseco, matroutp):
    i, j = numba.cuda.grid(2)
    if i < matroutp.shape[0] and j < matroutp.shape[1]:
        for k in range(matrfrst.shape[1]):
            matroutp[i, j] += matrfrst[i, k] * matrseco[k, j]


matrfrst = zeros((625, 10))
matrseco = zeros((10, 500))
matroutp = zeros((625, 500))

timeinit = time.time()
matrmult(matrfrst, matrseco, matroutp)
print (time.time() - timeinit) * 1e3

timeinit = time.time()
matrmultjitt(matrfrst, matrseco, matroutp)
print (time.time() - timeinit) * 1e3

timeinit = time.time()
matrmultgpuu(matrfrst, matrseco, matroutp)
print (time.time() - timeinit) * 1e3


