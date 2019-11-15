import math
import random
import time
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

#%% initial parameter for particles and plotting
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
max_diam = 6.0 

#%% Initial Packing of Streambed 
def bed_complete(pack_idx):
    """ Boolean check for whether 1D bed is fully packed or not""" 
    # similarly, if np.count_nonzero(bed_space) == x_max
    if pack_idx >= x_max:
        return 1
    else: return 0
  
    
def pack_bed(random_diam, bed_particles, particle_id, pack_idx):
    """ Add a new particle to the particle set. Ensure parameters are
    maintained and packing requirements are met """
    # TODO: fix being over plot by 5 max -- can we make a new min? 
    bed_particles[particle_id] = [random_diam, pack_idx]
    # update build parameters
    pack_idx += random_diam
    particle_id += 1
    # TODO: related; use remaining_space to understand tight packing near edge
    remaining_space = x_max - pack_idx
    return particle_id, pack_idx, remaining_space


def build_streambed(bed_particles, min_diam, max_diam, current_id, pack_idx):
    """ Build streambed until packed. Return n-2 array of particle diameter 
    and starting idx, where array index = particle id """
    while True:
        random_diam = random.randint(min_diam, max_diam)
        current_id, pack_idx, rem_space = pack_bed(random_diam, bed_particles, current_id, pack_idx)
        # TODO: clarify if x_max should be altered or bed repacked when x_max not met exactly
        if bed_complete(pack_idx):
            break
        else: continue
    # strip zero element particles tuples
    valid = ((bed_particles==0).all(axis=(1)))
    bed_particles = bed_particles[~valid]
    return bed_particles
 
# TODO: Make this part not look gross?
### Calls/Script Section
pack_idx = 0
max_particles = int(math.ceil(x_max/min_diam))
bed_particles = np.zeros([max_particles, 2],dtype='int')
current_id = 0

bed_particles = build_streambed(bed_particles, min_diam, max_diam, current_id, pack_idx)   
radius_array = np.asarray((bed_particles[:,0] / 2.0), dtype=float)
### End Calls/Script Section

# TODO: encapsulate plotting in function
fig = plt.figure(1)
fig.set_size_inches(10.5, 6.5)
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlim((-2, x_max + 2))
ax.set_ylim((0, x_max/4 + 2))
resolution = 50

x_center = (bed_particles[:,0] + bed_particles[:,1]) - radius_array
y_center_bed = np.zeros(np.size(x_center))
plt.rcParams['image.cmap'] = 'gray'
## This method of plotting circles comes from Stack Overflow questions\32444037
## Note that the patches won't be added to the axes, instead a collection will.
patches = []
for x1, y1, r in zip(x_center, y_center_bed, radius_array):
    circle = Circle((x1, y1), r)
    patches.append(circle)
p = (PatchCollection(patches, color="#BDBDBD", alpha=0.9, linewidths=(0, )))
ax.add_collection(p)

    
#%% Bed is built. Place n particles in avaliable vertices

def determine_num_particles(pack_frac, num_vertices):
    """ determine the number of particles to introduce into model, based
    on the packing fraction and number of avaliable vertices """
    
    num_particles = num_vertices * pack_frac
    num_particles = int(math.ceil(num_particles))
#    print(num_vertices, num_particles)
    return num_particles


def place_model_particles(vertex_idx):
    """ Randomly choose vertices from vertex_idx to place n particles. 
    Returns n-3 array containing the center coordinate, diameter, elevation
    of each individual particle """
    
    num_vertices = np.size(vertex_idx)
    already_selected = [False] * num_vertices
    num_particles = determine_num_particles(Pack, num_vertices)
     
    for particle in range(num_particles):
        
        random_diam = random.randint(min_diam, max_diam)
        # select vertex and ensure it has not previously been selected
        random_idx = random.randint(1, np.size(vertex_idx)-1)
        while already_selected[random_idx]:
            random_idx = random.randint(1, np.size(vertex_idx)-1)
        random_vertex = vertex_idx[random_idx]
       
        # FOR TESTING: plots chosen vertexes
        plt.axvline(x=random_vertex, color='g', linestyle='-')
#        if particle_fits():
#            # update particle cell
#        else: 
#            fu_resize()
        # CHECK fit, -- resize -- then place
        # UPDATE particle center, diameter and elevation
        
# boolean array with vertex avalibility based on packed bed 
avaliable_vertices = np.zeros(x_max, dtype=bool)
avaliable_vertices[bed_particles[:,1]] = 1 #TODO: exclude start idx for first particle
# x-indexes of avaliable spots
vertex_idx = np.transpose(np.nonzero(avaliable_vertices))


### FOR TESTING: Plots the avaliable vertex lines 
for xc in vertex_idx:
    plt.axvline(x=xc, color='b', linestyle='-')
###

place_model_particles(vertex_idx)

plt.show()       
# show places particles - display with different colour

