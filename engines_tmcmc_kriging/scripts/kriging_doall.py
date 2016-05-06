import time
start = time.time()

import numpy as np
import sys
from scipy.spatial import Delaunay
from sklearn.neighbors import NearestNeighbors

#lib = 'pyKriging'
lib = 'scikit'
if lib == 'pyKriging':
	import pyKriging
	from pyKriging.utilities import loadModel
elif lib == 'scikit':
	from sklearn.gaussian_process import GaussianProcess


################################### Functions: ######################################

def train(X,Y):
	print '---------------> Training the model...'

	if lib == 'pyKriging':
		model = kriging(X,Y)
		model.train()
	elif lib == 'scikit':
		model = GaussianProcess()
#		model = GaussianProcess(
#				regr='linear', corr='generalized_exponential',
#			theta0=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
#			thetaL=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
#			thetaU=np.concatenate( ( 10.0 *np.ones((1,dim)), 1.0*np.ones((1,1)) ), axis=1)
#			)
#		model = GaussianProcess(
#				regr='quadratic', corr='generalized_exponential',
#			theta0=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
#			thetaL=np.concatenate( ( 0.01 *np.ones((1,dim)), 0.01*np.ones((1,1)) ), axis=1),
#			thetaU=np.concatenate( ( 10.0 *np.ones((1,dim)), 1.0*np.ones((1,1)) ), axis=1)
#			)
		model.fit(X,Y)

	return model

def predict(model, X):
	if lib == 'pyKriging':
		Y = model.predict(pred_x)
		mse = model.calcuatemeanMSE()
		err = mse[0]
	elif lib == 'scikit':
		Y, mse = model.predict(X, eval_MSE=True)
		Y = Y[0]
		err = mse[0]

	return Y, abs(err/Y)

def build_hull(points):
	"""
	`hull` is either the `MxK` array of the 
	coordinates of `M` points in `K`dimensions for which Delaunay triangulation
	will be computed
	"""

	hull = Delaunay(points)

	return hull


def in_hull(p, hull):
	# code from http://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
	"""
	Test if points in `p` are in `hull`
	`p` should be a `NxK` coordinates of `N` points in `K` dimensions
	`hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
	coordinates of `M` points in `K`dimensions for which Delaunay triangulation
	will be computed
	"""

	if not isinstance(hull,Delaunay):
		hull = Delaunay(hull)
		
	return hull.find_simplex(p)>=0

####################################################################################

#start = time.time()

# Parse arguments

print '---------------> Reading arguments...'
my_arg = sys.argv[1].split()
gen_id = int(my_arg[0])
chain_id = int(my_arg[1])
task_id = int(my_arg[2])
D = int(my_arg[3])

leader = []
for i in range(4,D+4):
	leader.append(float(my_arg[i]))

pred_x = []
for i in range(D+4,len(my_arg)):
	pred_x.append(float(my_arg[i]))

# The Kriging support points
print '---------------> Getting support points...'
filename = 'full_db_%03d.txt' % (gen_id-1)
f = open(filename, 'r')
samples = [ map(float,line.split()) for line in f ]
f.close()
samples = np.array(samples)
(n,dim) = np.shape(samples)

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

samplesX = samplesX[true_indices,:]#.reshape((n,dim))
samplesY = samplesY[true_indices  ]#.reshape((n    ))

# Select only the unique points (needed for the scikit)
u, u_indices = np.unique(samplesY, return_index=True)
n = len(u_indices)

samplesX = samplesX[u_indices,:]#.reshape((n,dim))
samplesY = samplesY[u_indices  ]#.reshape((n    ))

# Find nearest neighbors
print '---------------> Computing neighbors...'
nn = min(max(5*(dim+1),150), n)
nbrs = NearestNeighbors(n_neighbors=nn, algorithm='auto', metric='euclidean').fit(samplesX)
nn_distances, nn_indices = nbrs.kneighbors(leader)

# Update samples
samplesX = samplesX[nn_indices,:].reshape((nn,dim))
samplesY = samplesY[nn_indices  ].reshape((nn    ))
print 'COM (hull):', np.sum(samplesX,0)/nn


# Now that we have our initial data, we can create an instance of a Kriging model
model = train(samplesX, samplesY)


# Is the prediction point in the convex hull?
if in_hull(pred_x, samplesX):
	try:
		pred_y, err = predict(model,pred_x)
	except:
		pred_y = -1
		err = 1e12
		pass
	print '---------------> Predicted for', pred_x, ':', pred_y, 'with error:', err

else:
	pred_y = -1
	err = 1e12
	print 'Point', pred_x, 'not in the convex hull for (', gen_id, chain_id, ')! Returning a big-big error.'

end = time.time()
print 'Prediction: %.6f' % pred_y
print 'Error: %.6f' % err
print 'Elapsed time: %.6f secs' % (end-start)

