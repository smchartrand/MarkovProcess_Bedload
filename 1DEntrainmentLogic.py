import math
import random
import time
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

############################################################################### 
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

############################################################################### 
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
 
    return particle_id, pack_idx


def build_streambed(bed_particles, min_diam, max_diam, current_id, pack_idx):
    """ Build streambed until packed. Return n-2 array of particle diameter 
    and starting idx, where array index = particle id """
    while True:
        random_diam = random.randint(min_diam, max_diam)
        current_id, pack_idx = pack_bed(random_diam, bed_particles, current_id, pack_idx)
        if bed_complete(pack_idx):
            break
        else: continue
    
    # update x_max once the bed is complete
    x_max = bed_particles[current_id-1][0] + bed_particles[current_id-1][1]
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


###############################################################################   
#%% Bed is built. Place/create n particles in avaliable vertices
##### NOTE: THIS FUNCTION SHOULD END UP IN ANOTHER SECTION. NOT APPRO HERE
def plot_stream(bed_particles, radius_array, chosen_vertex, x_lim, y_lim):
    
    fig = plt.figure(1)
    fig.set_size_inches(10.5, 6.5)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    # NOTE: xlim and ylim modified for aspec ratio -- WIP
    ax.set_xlim((-2, x_lim))
    ax.set_ylim((0, y_lim))
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
    ### FOR TESTING: Plots the avaliable vertex lines 
    for xc in vertex_idx:
        plt.axvline(x=xc, color='b', linestyle='-')
    for green in chosen_vertex:
        plt.axvline(x=green, color='g', linestyle='-')
    ### 
    plt.show()
    return
    
def determine_num_particles(pack_frac, num_vertices):
    """ determine the number of particles to introduce into model, based
    on the packing fraction and number of avaliable vertices """
    
    num_particles = num_vertices * pack_frac
    num_particles = int(math.ceil(num_particles))
#    print(num_vertices, num_particles)
    return num_particles


def fit_particle(particle_id, chosen_vertex, diam):
    """ Fit a particle of size diam at location of chosen_vertex. If diam is
    too large, resize until it will fit. Return the resulting center coord, 
    diameter and elevation of the fitted particle """
    # grab information about left and right bed particles
    left_bed_p = 0
    right_bed_p = 0
    
    # resize particle diam until it 'fits'
    while ~particls_fits():
            fu_resize()
            
    
    return p_center, p_diam, p_elev

def place_model_particles(vertex_idx):
    """ Randomly choose vertices from vertex_idx to place n particles. 
    Returns n-3 array containing the center coordinate, diameter and elevation
    of each individual particle """
    
    num_vertices = np.size(vertex_idx)
    already_selected = [False] * num_vertices
    num_particles = determine_num_particles(Pack, num_vertices)
    model_particles = np.zeros([max_particles, 3],dtype='int')
    # FOR TESTING:
    chosen_vertex = np.zeros(num_particles)
    
    
    for particle in range(num_particles):
        # select vertex and ensure it has not previously been selected
        random_idx = random.randint(1, np.size(vertex_idx)-1)
        while already_selected[random_idx]:
            random_idx = random.randint(1, np.size(vertex_idx)-1)
        vertex = vertex_idx[random_idx]
        # FOR TESTING: 
        chosen_vertex[particle] = vertex
        
        
#        random_diam = random.randint(min_diam, max_diam)
#        p_center, p_diam, p_elev = fit_particle(particle, chosen_vertex, \  
#                                                random_diam)
#        #  update cell in model_particles
#        model_particles[particle_id][0] = p_center
#        model_particles[particle_id][1] = p_diam
#        model_particles[particle_id][2] = p_elev
        
    return chosen_vertex
 
       
avaliable_vertices = np.zeros(x_max, dtype=bool)
# boolean array with vertex avalibility based on packed bed 
avaliable_vertices[bed_particles[1:,1]] = 1
# x-indexes of avaliable spots
vertex_idx = np.transpose(np.nonzero(avaliable_vertices))
chosen_vertex = place_model_particles(vertex_idx)
plot_stream(bed_particles, radius_array, chosen_vertex, 100, 100/4)

      
# show places particles - display with different colour
############################################################################### 
