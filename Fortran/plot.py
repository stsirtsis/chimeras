import csv
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
from pylab import savefig
# make these smaller to increase the resolution
dx, dy = 1, 1

# generate 2 2d grids for the x & y bounds
y, x = np.mgrid[slice(1, 101, dy),
                slice(1, 101, dx)]

z = np.sin(x)**10 + np.cos(10 + y*x) * np.cos(x)

with open('testing.dat', 'rb') as csvfile:
	spamreader=csv.reader(csvfile, delimiter=',', quotechar='|')
	i=0
	for row in spamreader:		
		j=0		
		for num in row:
			if (num != ""):
				z[i,j]=float(num)			
				j+=1
		i+=1
			


# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
z = z[:-1, :-1]
levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.get_cmap('inferno')

fig = plt.plot()

# contours are *point* based plots, so convert our bound into point
# centers
cf = plt.contourf(x[:-1, :-1] + dx/2.,
                  y[:-1, :-1] + dy/2., z, levels=levels,
                  cmap=cmap)
plt.colorbar(cf)
plt.title('Mean phase velocity')
savefig('mean_phase_velocity.png')
#plt.show()
