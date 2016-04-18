import numpy as np
import sys
from sklearn.neighbors import NearestNeighbors
import pickle

gen_id = int(sys.argv[1])
chain_id = int(sys.argv[2])

lib = 'pyKriging'
#lib = 'scikit'
if lib == 'pyKriging':
#	import pyKriging
	from pyKriging.krige import kriging
	from pyKriging.utilities import saveModel
elif lib == 'scikit':
	from sklearn.gaussian_process import GaussianProcess


# The Kriging support points
print '---------------> Getting support points...'
#filename = 'kriging_train_%03d_%03d.txt' % (gen_id, chain_id)
filename = 'full_db_%03d.txt' % (gen_id-1)
f = open(filename, 'r')
samples = [ map(float,line.split()) for line in f ]
samples = np.array(samples)
(n,dim) = np.shape(samples)
f.close()

surrogateId = samples[:,  dim-1].astype(int)
samplesY    = samples[:,  dim-2]
samplesX    = samples[:,0:dim-2]

# redefine shapes
dim = dim-2
n = n-sum(surrogateId)

# Take only real model evaluations
true_indices = []
for i in range(len(surrogateId)):
	if surrogateId[i] == 0:
		true_indices.append(i)
# Update samples
samplesX = samplesX[true_indices,:].reshape((n,dim))
samplesY = samplesY[true_indices  ].reshape((n    ))

## soften bad points
#mY = 2*np.min(samplesY)
#for i in range(0,n):
#	if samplesY[i] <= -1e12:
#		samplesY[i] = mY

print '---------------> Reading chain leader...'
filename = 'leader_%03d_%03d.txt' % (gen_id, chain_id)
f = open(filename, 'r')
leader = [ map(float,line.split()) for line in f ]
leader = np.array(leader)
f.close()

# Find nearest neighbors
print '---------------> Computing neighbors...'
nn = min(max(5*(dim+1),60), n)
nbrs = NearestNeighbors(n_neighbors=nn, algorithm='ball_tree', metric='euclidean').fit(samplesX)
nn_distances, nn_indices = nbrs.kneighbors(leader)

# Update samples
samplesX = samplesX[nn_indices,:].reshape((nn,dim))
samplesY = samplesY[nn_indices  ].reshape((nn    ))
print 'COM (hull):', np.sum(samplesX,0)/nn


def train(X,Y):
	print '---------------> Training the model...'

	if lib == 'pyKriging':
		model = kriging(X,Y)
		model.train()
#		print '---------------> MSE (mean, std):', model.calcuatemeanMSE()
	elif lib == 'scikit':
		model = GaussianProcess()
#		model = GaussianProcess(regr='linear', corr='generalized_exponential',
#			theta0=(1,2,3,4,5), thetaL=(0.1,0.1,0.1,0.1,0.1), thetaU=(2,3,4,5,6))
		model.fit(X,Y)

	return model

def save(model,filename, hull, hull_filename):
	print '---------------> Saving the model...'

	if lib == 'pyKriging':
		saveModel(model, filename)
	elif lib == 'scikit':
		pickle.dump(model, open(model_filename, 'wb'))

	pickle.dump(hull, open(hull_filename, 'wb'))


def build_hull(points):
	"""
	`hull` is either the `MxK` array of the 
	coordinates of `M` points in `K`dimensions for which Delaunay triangulation
	will be computed
	"""

	from scipy.spatial import Delaunay

	hull = Delaunay(points)

	return hull


# Now that we have our initial data, we can create an instance of a Kriging model
model = train(samplesX, samplesY)
hull = build_hull(samplesX)

model_filename = 'kriging_model_%03d_%03d.pkl' % (gen_id, chain_id)
hull_filename = 'kriging_hull_%03d_%03d.pkl' % (gen_id, chain_id)

save(model, model_filename, hull, hull_filename)
