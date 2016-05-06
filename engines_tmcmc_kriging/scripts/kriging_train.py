import numpy as np
import sys
from sklearn.neighbors import NearestNeighbors
import pickle
import time
import os


##################################### Functions: #######################################

def train(X,Y):
	print '---------------> Training the model...'

	if lib == 'pyKriging':
		model = kriging(X,Y)
		model.train()
#		print '---------------> MSE (mean, std):', model.calcuatemeanMSE()
	elif lib == 'scikit':
#		model = GaussianProcess()
#		model = GaussianProcess(regr='linear', corr='squared_exponential')
		model = GaussianProcess(
				regr='linear', corr='generalized_exponential',
			theta0=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
			thetaL=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
			thetaU=np.concatenate( ( 10.0 *np.ones((1,dim)), 1.0*np.ones((1,1)) ), axis=1)
			)
#		model = GaussianProcess(
#				regr='quadratic', corr='generalized_exponential',
#			theta0=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
#			thetaL=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
#			thetaU=np.concatenate( ( 10.0 *np.ones((1,dim)), 1.0*np.ones((1,1)) ), axis=1)
#			)
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

###################################################################################

start = time.time()

# Parse arguments
print '---------------> Reading arguments...'
my_arg = sys.argv[1].split()
gen_id = int(my_arg[0])
chain_id = int(my_arg[1])

leader = []
for i in range(2,len(my_arg)):
	leader.append(float(my_arg[i]))
print leader


#lib = 'pyKriging'
lib = 'scikit'
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

## take only the last generation
#nlast = n/gen_id
#samples = samples[n-nlast:n,:]
#n = nlast

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

# Select only the unique points (needed for the scikit)
u, u_indices = np.unique(samplesY, return_index=True)
n = len(u_indices)
samplesX = samplesX[u_indices,:].reshape((n,dim))
samplesY = samplesY[u_indices  ].reshape((n    ))

#print '---------------> Reading chain leader...'
#filename = 'leader_%03d_%03d.txt' % (gen_id, chain_id)
#f = open(filename, 'r')
#leader = [ map(float,line.split()) for line in f ]
#leader = np.array(leader)
#f.close()
#
## remove the unnecessary file
#try:
#	os.remove(filename)
#except OSError:
#	pass

# Find nearest neighbors
print '---------------> Computing neighbors...'
nn = min(max(5*(dim+1),150), n)
try:
	nbrs = NearestNeighbors(n_neighbors=nn, algorithm='auto', metric='euclidean').fit(samplesX)
except: # catch *all* exceptions
	sys.exit(-1)
nn_distances, nn_indices = nbrs.kneighbors(leader)

# Update samples
samplesX = samplesX[nn_indices,:].reshape((nn,dim))
samplesY = samplesY[nn_indices  ].reshape((nn    ))
print 'COM (hull):', np.sum(samplesX,0)/nn

## print neighbors
#filename = 'neighbors_%03d_%03d.txt' % (gen_id, chain_id)
#np.savetxt(filename, samplesX)

# Now that we have our initial data, we can create an instance of a Kriging model
try:
	model = train(samplesX, samplesY)
except: # catch *all* exceptions
	sys.exit(-2)

hull = build_hull(samplesX)

model_filename = 'kriging_model_%03d_%03d.pkl' % (gen_id, chain_id)
hull_filename = 'kriging_hull_%03d_%03d.pkl' % (gen_id, chain_id)

save(model, model_filename, hull, hull_filename)

end = time.time()
print 'Elapsed time: %.6f secs' % (end-start)

