def plot(filename,Nth=2,i=0,j=1):
	"""
	Simple demo of a scatter plot.
	"""
	import numpy as np
	import matplotlib.pyplot as plt

	#data = np.loadtxt('curgen_db_007.txt')
	data = np.loadtxt(filename)
	N=data.shape[0];

	print N 

	x = data[:,i]
	y = data[:,j]
	colors = data[:,Nth]
	plt.scatter(x, y, c=colors, alpha=0.5)
	#plt.scatter(x, y, c=colors, s=1.0, alpha=0.5)
	plt.show()
