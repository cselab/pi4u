def plot_gen(filename, dim=2, i=1, j=2, namei="", namej="", show=1, save=0):
	import numpy as np
	import matplotlib.pyplot as plt

	data = np.loadtxt(filename)
	N=data.shape[0];

	print(N)

	x = data[:,i-1]
	y = data[:,j-1]
	colors = data[:,dim]

	fig = plt.figure()

	plt.scatter(x, y, c=colors, marker='o', s=10.0, alpha=0.3)

	if (namei==""):
		vari_name='var-'+str(i)
	else:
		vari_name=namei

	if (namej==""):
		varj_name='var-'+str(j)
	else:
		varj_name=namej

	plt.xlabel(vari_name, fontsize=18)
	plt.ylabel(varj_name, fontsize=18)

#	plt.xlabel('var-'+str(i), fontsize=18)
#	plt.ylabel('var-'+str(j), fontsize=18)
	plt.colorbar()

	if (show == 1):
		plt.show()

	if (save == 1):
		fig.savefig(filename+'.png', dpi=fig.dpi)

	return
