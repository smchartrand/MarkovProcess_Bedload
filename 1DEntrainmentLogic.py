from __future__ import division
import math
import random
import time
import numpy as np
import sympy as sy

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

# FOR TESTING: 
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

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
    
    center = pack_idx + (random_diam/2)

    elevation = 0
    bed_particles[particle_id] = [random_diam, pack_idx, center, elevation]
  
    # update build parameters
    pack_idx += random_diam
    particle_id += 1
 
    return particle_id, pack_idx


def build_streambed(bed_particles, min_diam, max_diam, current_id, pack_idx):
    """ Build streambed until packed. Return n-4 array of particle diameter 
    starting idx, x coord and y coord. Array index = particle id """
    global x_max
    while True:
#        random_diam = random.randint(min_diam, max_diam)
        random_diam = 5.0
        current_id, pack_idx = pack_bed(random_diam, bed_particles, current_id, pack_idx)
        if bed_complete(pack_idx):
            break
        else: continue
    
    # bed can be packed +- 8mm from the default x_max of 500 depending on the 
    # packing pattern -- therefore update x_max once bed is complete to new +- 8mm size
    x_max = int(math.ceil(bed_particles[current_id-1][0] + bed_particles[current_id-1][1]))
    # strip zero element particles tuples
    valid = ((bed_particles==0).all(axis=(1)))
    bed_particles = bed_particles[~valid]

    return bed_particles
 
    
### Calls/Script Section
pack_idx = 0
max_particles = int(math.ceil(x_max/min_diam))
bed_particles = np.zeros([max_particles, 4],dtype='float')
current_id = 0

bed_particles = build_streambed(bed_particles, min_diam, max_diam, current_id, pack_idx)   
radius_array = np.asarray((bed_particles[:,0] / 2.0), dtype=float)
### End Calls/Script Section

###############################################################################   
#%% Bed is built. Place/create n particles in avaliable vertices
##### NOTE: plot_stream SHOULD END UP IN ANOTHER SECTION. NOT APPRO HERE


def plot_stream(bed_particles, model_particles, radius_array, chosen_vertex, x_lim, y_lim):
    """ Plot the complete stream from 0,0 to x_lim and y_lim. Bed particles 
    are plotted as light grey and model particles are dark blue. Allows
    for closer look at state of a subregion of the stream during simulation """
    plt.clf()
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
    
    x_center_m = model_particles[:,0]
    y_center_m = model_particles[:,2]
    patches1 = []
    for x2, y2, r in zip(x_center_m, y_center_m, model_particles[:,1]/2):
        circle = Circle((x2, y2), r)
        patches1.append(circle)
    p_m = (PatchCollection(patches1, facecolors='black', edgecolors='black'))
    ax.add_collection(p_m)
    ### FOR TESTING: Plots the avaliable and chosen vertex lines 
    # uncomment for loop sections to draw the corresponding lines on figure
#    for xc in vertex_idx:
#        plt.axvline(x=xc, color='b', linestyle='-')
#    for xxc in total_vertices:
#        plt.axvline(x=xxc, color='m', linestyle='-')
#    for green in chosen_vertex:
#        plt.axvline(x=green, color='g', linestyle='-')
    ### 
    plt.show()
    return
    

def determine_num_particles(pack_frac, num_vertices):
    """ Return the number of particles to introduce into model. Choice is based
    on the packing fraction and number of avaliable vertices """
    
    num_particles = num_vertices * pack_frac
    num_particles = int(math.ceil(num_particles))

    return num_particles


def place_particle(g1, g2, p_diam):
    """ put function definition here """
    
    # initialize neighbour (g1 and g2) and placed particle information using info from the storage arrays
    x1 = g1[2]
    x2 = g2[2]
    y1 = g1[3]
    y2 = g2[3]
    r1 = g1[0] / 2
    r2 = g2[0] / 2
    rp = p_diam / 2 
    
    # define symbols for symbolic system solution using SymPy
    x3, y3 = sy.symbols('x3 y3')
    
    # create the symbolic system equations
    eq1 = sy.Eq(sy.sqrt((x1-x3)**2 + (y1-y3)**2)-r1-rp, 0)
    eq2 = sy.Eq(sy.sqrt((x2-x3)**2 + (y2-y3)**2)-r2-rp, 0)
    
    # solve the system of equations
    sol_dict = sy.solve((eq1, eq2), (x3, y3))
        
    # iterate into the solution dictionary to recieve particle center (x,y)
    p_x = (sol_dict[1][0])
    p_y = (sol_dict[1][1])
 
    return p_x, p_diam, p_y

    

def find_neighbours_of(idx):
    """ For a given vertex, return the two neighbouring grains """ 
    # TODO: Need to update function to consider model particles _as well as_ 
    # bed particles as possible contact neighbours when being entrained
    grain_x = bed_particles[idx]
    grain_y = bed_particles[idx+1]
    
    g1 = grain_x
    g2 = grain_y
            
    return g1, g2


def place_model_particles(vertex_idx, bed_particles):
    """ Randomly choose vertices from vertex_idx to place n particles. 
    Returns n-3 array containing the center coordinate, diameter and elevation
    of each individual particle """

    # create a boolean area representing the avaliability of the vertices
    num_vertices = np.size(vertex_idx)
    already_selected = [False] * num_vertices
    
    # determine the number of model particles that should be introduced into the stream bed
    num_particles = determine_num_particles(Pack, num_vertices)
    # create an empty n-4 array to store model particle information
    model_particles = np.zeros([num_particles, 4], dtype='float')
    
    #### FOR TESTING:
    chosen_vertex = np.zeros(num_particles)
    ####
    
    new_vertices = np.zeros((num_particles,2))
    
    for particle in range(num_particles):  
#        random_diam = random.randint(min_diam, max_diam)
        random_diam = 5.0
        
        # the following lines select a vertex to place the current particle at, 
        # and ensure that it is not already occupied by another particle
        random_idx = random.randint(0, np.size(vertex_idx)-1)
        while already_selected[random_idx]:
            random_idx = random.randint(0, np.size(vertex_idx)-1)
        already_selected[random_idx] = True
        vertex = vertex_idx[random_idx]
        
        #### FOR TESTING: 
        chosen_vertex[particle] = vertex
        ####
        
        # once a vertex is chosen, this function identifies the neighbours
        g1, g2 = find_neighbours_of(random_idx)
        # get the particles initial x, y and diameter information in the bed
        p_x, p_diam, p_y = place_particle(g1, g2, random_diam)
        
        # intialize the particle information
        model_particles[particle][0] = p_x
        model_particles[particle][2] = p_y
        model_particles[particle][1] = p_diam
        model_particles[particle][3] = particle # id number for each particle
        
        # store each particles vertex information in new_vertices
        new_vertices[particle][0] = (p_x) - (p_diam/2)
        new_vertices[particle][1] = (p_x) + (p_diam/2)
    
    # flatten both vertex arrays so that they can be merged together by insort()
    new_vertices = new_vertices.flatten()
    old_vertices = vertex_idx.flatten()
    
    # retrieve all vertices present in the stream. Currently does not filter out unavaliable vertices
    total_vertices = insort(new_vertices, old_vertices)
    
    return model_particles, chosen_vertex, total_vertices

# from: https://stackoverflow.com/questions/12427146/combine-two-arrays-and-sort
def insort(a, b, kind='mergesort'):
    # took mergesort as it seemed a tiny bit faster for my sorted large array try.
    c = np.concatenate((a, b)) # we still need to do this unfortunatly.
    c.sort(kind=kind)
    flag = np.ones(len(c), dtype=bool)
    np.not_equal(c[1:], c[:-1], out=flag[1:])
    return c[flag]

 
    
### Calls/Script Section
avaliable_vertices = np.zeros(x_max, dtype=bool)

avaliable_vertices[bed_particles[1:,1].astype(int)] = 1
# x-indexes of avaliable vertices to place model particles at
bed_vertices = np.transpose(np.nonzero(avaliable_vertices))

model_particles, chosen_vertex, total_vertices = place_model_particles(bed_vertices, bed_particles)
plot_stream(bed_particles, model_particles, radius_array, chosen_vertex, 500, 100/4)
### End Calls/Script Section

############################################################################### 
#%% Model particles have been placed on the bed; stream build is complete.
# Divide stream into sampling regions and build sampling array to store data
bed_sampreg = 25
# ------------------ below is WIP, pulled directly from Quasi2D:
## Data storage arrays
#Nu_out_Store = np.zeros([Bed_sampreg, Loops], dtype=int, order='F')#[0]
## Bed sampling boundaries in the x-direction
#BS_boundaries = (XCoordinates_Orig + (x_max / Bed_sampreg) *
#                 np.arange(0, Bed_sampreg + 1, dtype=int))
#SSamp_len = len(BS_boundaries)
## Index to hold boundaries between subsampling regions to facilitate random
## sampling from within each region.
#SubSampl_idx = np.zeros([1, SSamp_len], dtype=int, order='F')
#xB_idx = np.zeros([1, 2], dtype=int, order='F')
#yB_idx = np.zeros([1, 2], dtype=int, order='F')

###############################################################################
# number of model iterations  
n_iterations = 10
lambda_1 = 1
# Particle travel time minimum (t).
T_pmin = 0
# Particle travel time maximum (t).
T_pmax = 1.0
e_events_store = np.zeros(n_iterations)

# TODO: Need to move n particles (n = e_events) per subregion, not from whole stream
def move_model_particles(e_events, rand_particles):
    # uniform distribution constrained by Fathel et al., 2015.
#    if n == 0 and Nu_in > 0:
#        continue
    global total_vertices
    T_p_init1 = (np.random.uniform(T_pmin, T_pmax, e_events).
                 reshape(1, e_events))
    # https:\\stackoverflow.com\questions\2106503\
    # Now distribute the random variables over a pseudo exponential
    # distribution based on inverse transform sampling.
    T_p_init2 = np.log(1 - T_p_init1) / (-lambda_1)
    T_p_alt2 = -10 * np.log(1-T_p_init1)
    # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
    # For now I am multuplying by 10 so units are consistent (). Rounding
    # the output to map the entrained particles to the 2D grid.
    L_x_init = np.round((T_p_init2 ** 2) * 10 * 2, 1)
    
    # update the particle information to be original_x + hop distance
    rand_particles[:,0] = rand_particles[:,0] + L_x_init

    # now, update model_particle array for each of the randomly selected particles
    for idx, particle_id in enumerate(rand_particles[:,3]):
        model_particles[int(particle_id)] = rand_particles[idx]
    
    # create empty n-2 array to hold the vertices of the model particles
    new_vertices = np.zeros((num_particles,2))
    # for each model particle, store what would be its vertex values in new_vertices
    for idx, particle in enumerate(model_particles):
        new_vertices[idx][0] = particle[0] - particle[1]/2
        new_vertices[idx][1] = particle[0] + particle[1]/2
    
    # flatten both vertex arrays so that we can merge them together 
    new_vertices = new_vertices.flatten()
    old_vertices = bed_vertices.flatten()
    
    # get array of all vertices in the stream. Does _not_ account for availability
    total_vertices = insort(new_vertices, old_vertices)
    
    ## FOR TESTING:
    # sleep to more easily see the bed migration in matplotlib
    time.sleep(1)  
    ###
    
    plot_stream(bed_particles, model_particles, radius_array, chosen_vertex, 500, 100/4)


for step in range(n_iterations):
    # calculate the number of entrainment events per unit area of bed
    e_events = np.random.poisson(lambda_1, None)

    if e_events == 0:
        e_events = 1 #???
        
    # total number of entrained particles
    e_grains = e_events * bed_sampreg
    # randomly select model particles to entrain then convert to numpy array
    rand_particles = random.sample(model_particles, e_events)
    rand_particles = np.array(rand_particles)
    
    move_model_particles(e_events, rand_particles)
    e_events_store[step] = e_events
    
    


