import math
import random
import time
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle


Area = np.zeros([1, 1], dtype=int, order='F')
# Total area of placed circles. Initialize to 0.
AreaTotal = 0
# Center coord. of placed circles: r1:x,r2:y.
CenterCoord = np.zeros(1, dtype=int, order='F')
Diameter = np.zeros(1, dtype=int, order='F')
CenterElev = np.zeros(1, dtype=int, order='F')

Color_iter = 1000
Step_1 = 1
# Define origin for plotting.
XCoordinates_Orig = 0


# Packing density of the bed; ranges from > 0 to <~0.70.
Pack = 0.52
# Length of the domain in the streamwise direction in millimeters.
x_max = 500


# Spacing of nodes in millimeters.
Grid_Spacing = 1
# Minimum grain diameter in millimeters.
min_diam = 4.0
# Size of largest grain relative to the cross-stream length.
Size_Frac = 0.06
max_diam = 6.0 #y_max / CircDiam_Fac


# Build x array for grid.
X_array = np.arange(0, (x_max) + Grid_Spacing, Grid_Spacing, dtype=float)
Length2 = int(len(X_array))
XCoordinates = np.repeat(X_array[None], Length1, axis=0)

#%% Following cell concerned with 1D bed packing:
def bed_complete(pack_idx):
    """ boolean check for whether bed is fully packed or not"""
    # similarly, if np.count_nonzero(bed_space) == 500
    if pack_idx >= 500:
        return 1
    else: return 0
  
    
def pack_bed(random_diam, bed_particles, particle_id, pack_idx):
    """ add a new particle to the particle set. Ensure parameters are
    maintained and tight packing requirements are met """
    # TODO: fix being over plot by 5 max -- can we make a new min? 
    bed_particles[particle_id] = [random_diam, pack_idx]
    # update build parameters
    pack_idx += random_diam
    particle_id += 1
    
    # TODO: related; use remaining_space to understand tight packing near edge
    remaining_space = x_max - pack_idx
    return particle_id, pack_idx, remaining_space


def build_streambed(bed_particles, min_diam, max_diam, current_id, pack_idx):
    """ Build streambed until packed. Store each particle diamter and 
    starting idx (for plotting purposes) """
    while True:
        random_diam = random.randint(min_diam, max_diam)
        current_id, pack_idx, rem_space = pack_bed(random_diam, bed_particles, current_id, pack_idx)

        if bed_complete(pack_idx):
            break
        else: continue
    # strip zero element particles tuples
    valid = ((bed_particles==0).all(axis=(1)))
    bed_particles = bed_particles[~valid]
    return bed_particles
 

pack_idx = 0
max_particles = int(math.ceil(x_max/min_diam))
bed_particles = np.zeros([max_particles, 2],dtype='int')
current_id = 0

bed_particles = build_streambed(bed_particles, min_diam, max_diam, current_id, pack_idx)   
# FOR TESTING:
print(bed_particles)

# TODO: encapsulate plotting in function
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlim((-2, x_max + 2))
ax.set_ylim((0, x_max/2 + 2))
resolution = 50
    

radius_array = np.asarray((bed_particles[:,0] / 2.0), dtype=float)
print(radius_array)


#Radius_Array = Radius_Array.reshape(-1)
x_center = (bed_particles[:,0] + bed_particles[:,1]) - radius_array
y_center_bed = np.zeros(np.size(x_center))
#XCenter = XCenter.reshape(-1)
plt.rcParams['image.cmap'] = 'gray'
## This method of plotting circles comes from Stack Overflow questions\32444037
## Note that the patches won't be added to the axes, instead a collection will.
patches = []
for x1, y1, r in zip(x_center, y_center_bed, radius_array):
    circle = Circle((x1, y1), r)
    patches.append(circle)
p = (PatchCollection(patches, color="#BDBDBD", alpha=0.9, linewidths=(0, )))
#colors = Color_iter*np.random.rand(len(CenterCoord[0]))
#p = (PatchCollection(patches, cmap=cm.coolwarm,
#                     alpha=1.0, linewidths=(0, )))
#p.set_array(colors)
#p.set_clim([5, 950])
ax.add_collection(p)
plt.show()
    
    
