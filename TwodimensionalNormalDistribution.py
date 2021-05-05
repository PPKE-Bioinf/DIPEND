# coding: utf-8

import math
import numpy as np
import sys
# splot (1/(3*sqrt(6.28)))*exp(-0.5*((x-0)/3)**2)*(1/(3*sqrt(6.28)))*exp(-0.5*((y-0)/3)**2) w pm3d
res=5 # resolution in degrees
grid=int(1+360/res)
avg=0 # the average will always be the centre of the map first, will be adjusted later
# set the desired Ramachandran angles and their deviation here
#phi_avg=-60 
#psi_avg=-45
#std=15

phi_avg = int(sys.argv[1])
psi_avg = int(sys.argv[2])
std = int(sys.argv[3])
outfname = sys.argv[4]

# this will be the array with the final values
values=np.zeros((grid,grid))

# a function to take care of the circular nature of the Ramachandran map
def inrama (a):
    if a < -180:
        return a+360
    if a > 180:
        return a-360
    return a

# a sum of all values for scaling at the end (to add up to a probability of 1)
allvalues=0
offset=0

# Calculating the distribution as a product of two 1D distributions
for phi in range(-180,180,res):
    for psi in range (-180,180,res):
        z= (1/(std*math.sqrt(2*math.pi)))*math.exp(-0.5*((phi-avg)/std)**2)
        z*=(1/(std*math.sqrt(2*math.pi)))*math.exp(-0.5*((psi-avg)/std)**2)
        #print(x,y,z)
        # adjusting the angles to the desired region
        #phi_out=inrama(phi-phi_avg)
        #psi_out=inrama(psi-psi_avg)
        phi_out=inrama(phi+phi_avg)
        psi_out=inrama(psi+psi_avg)
        # calculating the coordinates in the grid
        x=int(phi_out/res)
        y=int(psi_out/res)

        # for very shallow distributions (large stdev) the values might not be 0 at the map boundaries.
        # we will subtract this from all values in order to use only the differences between map regions
        if phi == -180 and psi == -180:
            #print ("z at edge:",z)
            offset=z

        values[x][y]=z-offset
        allvalues += z-offset

cumulative_probability = 0.0

# writing the distribution to a file
# note that we write the map in the conventional format
# from -180 to 180 but use the adjusted values
with open(outfname, "w") as f:
    for phi in range(-180,180,res):
        for psi in range (-180,180,res):
            x=int(phi/res)
            y=int(psi/res)
            # using the scaled values ensures that we have a cumulative robability of 1
            cumulative_probability +=  values[x][y]/allvalues
            f.write("%d %d %.2e %.2e\n" % (phi,psi,values[x][y]/allvalues, cumulative_probability))

# you can plot the resulting file with gnuplot:
# splot 'proba.dat' w pm3d
