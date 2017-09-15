from plot import plot 

#plot('curgen_db_007.txt',Nth=2,i=0,j=1)

import sys
filename=sys.argv[1]
Nth=int(sys.argv[2])
i=int(sys.argv[3])
j=int(sys.argv[4])

print filename
print Nth
print i
print j

plot(filename,Nth,i,j)
