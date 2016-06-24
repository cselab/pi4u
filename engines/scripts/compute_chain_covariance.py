import time
import numpy as np
import sys
from sklearn.neighbors import NearestNeighbors


def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)


start = time.time()

# Read samples
# print '---------------> Reading samples...'
samples = np.loadtxt('dataset.txt')
(n,dim) = np.shape(samples)

testset = np.loadtxt('testset.txt')
(ntest,dim) = np.shape(testset)

# Find nearest neighbors
# print '---------------> Computing neighbors...'
nn = min(max(5*(dim+1),50), n)
nbrs = NearestNeighbors(n_neighbors=nn, algorithm='auto', metric='euclidean').fit(samples)
nn_distances, nn_indices = nbrs.kneighbors(testset)

# Compute covariances
res = np.zeros((ntest, dim*dim))
for k in range(ntest):
    s = samples[nn_indices[k,:],:].reshape((nn,dim))
    cc = np.cov(s.T)

#    OK = is_pos_def(cc)
#    if not(OK):
#		print 'Warning: Matrix #%d is not positive definite!' % k

    res[k,:] = cc.reshape((1,dim*dim))

# Dump result
my_fmt = ''
for i in range(dim*dim):
    my_fmt += '%le '
np.savetxt('nn_result.txt', res, fmt=my_fmt)

end = time.time()
print 'Elapsed time: %.6f secs' % (end-start)
