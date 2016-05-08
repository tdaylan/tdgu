
# coding: utf-8

# In[ ]:

def prca(matr):
    M = (matr - mean(matr.T, 1)).T
    eigl, eigt = linalg.eig(cov(M))
    tranmatr = dot(eigt.T, M).T
    return eigl, tranmatr, eigt


get_ipython().magic(u'matplotlib inline')

npara = 2
nsamp = 1000
matr = zeros((nsamp, 2))
matr[:, 0] = randn(nsamp)
matr[:, 1] = randn(nsamp) * 10.
eigl, tranmatr, eigt = prca(matr)

plt.scatter(matr[:, 0], matr[:, 1])
plt.show()

plt.scatter(tranmatr[:, 0], tranmatr[:, 1])
plt.show()

