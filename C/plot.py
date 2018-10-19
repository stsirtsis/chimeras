import sys
from os import listdir
from os.path import isfile,join
import csv
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
from pylab import savefig
from pathlib import Path

resFiles=list(sys.argv)[1:]
#print resFiles
timeVar=0.0

for f in resFiles:
    # make these smaller to increase the resolution
    dx, dy = 1, 1

    # generate 2 2d grids for the x & y bounds
    y, x = np.mgrid[slice(1, 101, dy),
    slice(1, 101, dx)]

    z = np.zeros((100,100))
    timeVar=f.split("_")[10]
    if (f.split("_")[1]=="MPV"):
        myimg=Path("Plots/mean_phase_velocity_"+timeVar+".png")
    elif (f.split("_")[1]=="POT"):
        myimg=Path("Plots/potential_"+timeVar+".png")
    if not myimg.is_file():
        with open(f, 'rb') as csvfile:
            spamreader=csv.reader(csvfile, delimiter=',', quotechar='|')
            i=0
            for row in spamreader:
                j=0
                for num in row:
                    if (num != ""):
                        #print i,j
                        z[i,j]=float(num)
                        j=j+1
                i=i+1

        # x and y are bounds, so z should be the value *inside* those bounds.
        # Therefore, remove the last value from the z array.
        z = z[:-1, :-1]
	#MaxNLocator.MAX_TICKS=1500
        levels = MaxNLocator(nbins=256).tick_values(z.min(), z.max())


        # pick the desired colormap, sensible levels, and define a normalization
        # instance which takes data values and translates those into levels.
        cmap = plt.get_cmap('jet')

        fig = plt.plot()

        # contours are *point* based plots, so convert our bound into point
        # centers
        cf = plt.contourf(x[:-1, :-1] + dx/2.,
                          y[:-1, :-1] + dy/2., z, levels=levels,
                          cmap=cmap)

        plt.colorbar(cf, format='%.3f', ticks=[np.amax(z), (np.amax(z)+np.amin(z))/2, np.amin(z)])
        if (f.split("_")[1]=="MPV"):
            plt.title('Mean phase velocity at time '+timeVar+'.')
            savefig('Plots/mean_phase_velocity_'+timeVar+'.png')
        elif (f.split("_")[1]=="POT"):
            plt.title('Potential at time '+timeVar+'.')
            savefig('Plots/potential_'+timeVar+'.png')
        plt.clf()
