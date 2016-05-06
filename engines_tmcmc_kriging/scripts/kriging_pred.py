import numpy as np
import sys
import pickle
import time
import os


################################### Functions: ######################################

def read(model_filename, hull_filename):
	if lib == 'pyKriging':
		model = loadModel(model_filename)
	elif lib == 'scikit':
		model = pickle.load(open(model_filename,'rb'))

	hull = pickle.load(open(hull_filename,'rb'))

	return model, hull

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

def in_hull(p, hull):
	# code from http://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
	"""
	Test if points in `p` are in `hull`
	`p` should be a `NxK` coordinates of `N` points in `K` dimensions
	`hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
	coordinates of `M` points in `K`dimensions for which Delaunay triangulation
	will be computed
	"""

	from scipy.spatial import Delaunay

	if not isinstance(hull,Delaunay):
		hull = Delaunay(hull)
		
	return hull.find_simplex(p)>=0

####################################################################################

start = time.time()

# Parse arguments

print '---------------> Reading arguments...'
my_arg = sys.argv[1].split()
gen_id = int(my_arg[0])
chain_id = int(my_arg[1])
task_id = int(my_arg[2])

pred_x = []
for i in range(3,len(my_arg)):
	pred_x.append(float(my_arg[i]))


#lib = 'pyKriging'
lib = 'scikit'
if lib == 'pyKriging':
	import pyKriging
	from pyKriging.utilities import loadModel
elif lib == 'scikit':
	from sklearn.gaussian_process import GaussianProcess

## Read the prediction point
#filename = 'kriging_pred_x_%03d_%03d_%03d.txt' % (gen_id, chain_id, task_id)
#f = open(filename, 'r')
#pred_x = [ map(float,line.split()) for line in f ]
#pred_x = pred_x[0]
#f.close()
#
## remove the unnecessary file
#try:
#	os.remove(filename)
#except OSError:
#	pass


# Read the model
model_filename = 'kriging_model_%03d_%03d.pkl' % (gen_id, chain_id)
hull_filename = 'kriging_hull_%03d_%03d.pkl' % (gen_id, chain_id)
model, hull = read(model_filename, hull_filename)

# Is the prediction point in the convex hull?
if in_hull(pred_x, hull):
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

#filename = 'kriging_pred_y_%03d_%03d_%03d.txt' % (gen_id, chain_id, task_id)
#f = open(filename , 'w')
#f.writelines( "%f\n%f\n" % (pred_y, err) )
#f.close()

end = time.time()
print 'Prediction: %.6f' % pred_y
print 'Error: %.6f' % err
print 'Elapsed time: %.6f secs' % (end-start)

