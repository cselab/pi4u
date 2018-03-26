import numpy as np
import sys

p = np.loadtxt( sys.argv[1] )

data = np.loadtxt( sys.argv[2] )

A     = p[0]
omega = p[1]
phi   = p[2]

t = data[:,0]

y = A*np.sin(omega*t+phi)

res = np.column_stack((t,y))

np.savetxt('result.txt', res )
