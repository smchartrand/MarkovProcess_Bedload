from __future__ import division
import math
import random
import time
import numpy as np
import sympy as sy
import copy
import collections
import parameters # Import the parameters defined in the parameters file

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

############################################################################### 
#%% Initial Packing of Streambed 
     
def pack_bed(random_diam, bed_particles, particle_id, pack_idx):
    """ This function will add a new particle to the bed_particle array
    that has random diameter, id=packing_idx, calculated center and elevation """
    
    center = pack_idx + (random_diam/2)

    elevation = 0
    bed_particles[particle_id] = [center, random_diam, elevation, pack_idx]
  
    # update build parameters
    pack_idx += random_diam
    particle_id += 1
 
    return particle_id, pack_idx


def build_streambed(bed_particles, diam):
    """ Build streambed until packed. Return n-4 array of particle diameter 
    starting idx, x coord and y coord. Array index = particle id """
    running_id = 0
    running_pack_idx = 0
    while True:
        running_id, running_pack_idx = pack_bed(diam, bed_particles, running_id, running_pack_idx)
        if bed_complete(running_pack_idx):
            break
        else: continue
    
    # bed can be packed +- 8mm from the default x_max of 500 depending on the 
    # packing pattern -- therefore update x_max once bed is complete to new +- 8mm size
    x_max = int(math.ceil(bed_particles[running_id-1][1]+ bed_particles[running_id-1][3]))
    # strip zero element particles tuples
    valid = ((bed_particles==0).all(axis=(1)))
    bed_particles = bed_particles[~valid]
    
    # create array of bed_vertices based on bed_particles array
    avaliable_vertices = np.zeros(x_max, dtype=bool)

    avaliable_vertices[bed_particles[1:,3].astype(int)] = 1
    # x-indexes of avaliable vertices to place model particles at
    bed_vertices = np.transpose(np.nonzero(avaliable_vertices))

    return bed_particles, bed_vertices


def bed_complete(pack_idx):
    """ Provided a packing index (current # of particles placed in bed), this 
    function will return a Boolean indicating whether the bed packing is complete """ 
    # similarly, if np.count_nonzero(bed_space) == x_max
    if pack_idx >= parameters.x_max:
        return 1
    else: return 0

#%% Bed is built. Place/create n particles in avaliable vertices
##### NOTE: plot_stream SHOULD END UP IN ANOTHER SECTION. NOT APPRO HERE
    

def determine_num_particles(pack_frac, num_vertices):
    """ Return the number of particles to introduce into model. Choice is based
    on the packing fraction and number of avaliable vertices """
    
    num_particles = num_vertices * pack_frac
    num_particles = int(math.ceil(num_particles))

    return num_particles


def place_particle(g1, g2, p_diam):
    """ put function definition here """
    
    # initialize neighbour (g1 and g2) and placed particle information using info from the storage arrays
    x1 = g1[0]
    x2 = g2[0]
    y1 = g1[2]
    y2 = g2[2]
    r1 = g1[1] / 2
    r2 = g2[1] / 2
    rp = p_diam / 2 
    
    # define symbols for symbolic system solution using SymPy
    x3, y3 = sy.symbols('x3 y3')
    
    # create the symbolic system equations
    eq1 = sy.Eq(sy.sqrt((x1-x3)**2 + (y1-y3)**2)-r1-rp, 0)
    eq2 = sy.Eq(sy.sqrt((x2-x3)**2 + (y2-y3)**2)-r2-rp, 0)
    
    # solve the system of equations
    sol_dict = sy.solve((eq1, eq2), (x3, y3))
        
    # iterate into the solution dictionary to recieve new particle center (x,y)
    p_x = (sol_dict[1][0])
    p_y = (sol_dict[1][1])
 
    return p_x, p_diam, p_y

    
def find_neighbours_of(particle_idx, model_particles, bed_particles):
    """ For a given vertex, return the two neighbouring grains of highest level """ 

    # look for particles with center R (radius) away from particle_idx 
    left_neighbour_center = particle_idx - (parameters.set_diam / 2)
    right_neighbour_center = particle_idx + (parameters.set_diam / 2)
    
    # Look for left neighbour:
    # first search model particles for matches
    left_neighbour_idx = np.where(model_particles[:,0] == left_neighbour_center)
    
    # if no such model particle found, search bed particle array
    if left_neighbour_idx[0].size == 0: 
        left_neighbour_idx = np.where(bed_particles[:, 0] == left_neighbour_center)
        if left_neighbour_idx[0].size == 0:
            raise ValueError('No left neighbours found for _(idx).')
        else:
            left_neighbour = bed_particles[left_neighbour_idx]
            print('bed l n')
    else:
        left_neighbour = model_particles[left_neighbour_idx]
        print('model l n')
        
    # Look for right neighbour:
    # first search model particles for matches
    right_neighbour_idx = np.where(model_particles[:,0] == right_neighbour_center)
    
    if right_neighbour_idx[0].size == 0:
        right_neighbour_idx = np.where(bed_particles[:,0] == right_neighbour_center)
        if right_neighbour_idx[0].size == 0:
            raise ValueError('No right neighbours found for _(idx).')
        else:
            right_neighbour = bed_particles[right_neighbour_idx]
    else:
        right_neighbour = model_particles[right_neighbour_idx]

    #TODO: assert/check that only one neighbour has been found
            
    return left_neighbour[0], right_neighbour[0]


def set_model_particles(bed_vertices, bed_particles):
    """ Randomly choose vertices from vertex_idx to place n particles. 
    Returns n-4 array representing the resulting model particles where
    [0] = center coordinate,
    [1] = diameter,
    [2] = elevation,
    [3] = uid """

    # create a boolean area representing the avaliability of the vertices
    num_vertices = np.size(bed_vertices)
    already_selected = [False] * num_vertices
    
    # determine the number of model particles that should be introduced into the stream bed
    num_particles = determine_num_particles(parameters.Pack, num_vertices)
    # create an empty n-4 array to store model particle information
    model_particles = np.zeros([num_particles, 4], dtype='float')
    
    #### FOR TESTING:
    chosen_vertex = np.zeros(num_particles)
    ####
    
    new_vertices = np.zeros((num_particles,2))
    nulled_vertices = np.zeros((num_particles))
    
    
    for particle in range(num_particles):  
        
        # the following lines select a vertex to place the current particle at, 
        # and ensure that it is not already occupied by another particle
        random_idx = random.randint(0, np.size(bed_vertices)-1)
        while already_selected[random_idx]:
            random_idx = random.randint(0, np.size(bed_vertices)-1)
        already_selected[random_idx] = True
        vertex = bed_vertices[random_idx]
        
        #### FOR TESTING: 
        chosen_vertex[particle] = vertex
        ####
        
        # once a vertex is chosen, this function identifies the neighbours
        g1, g2 = find_neighbours_of(vertex, model_particles, bed_particles)
        # get the particles initial x, y and diameter information in the bed
        p_x, p_diam, p_y = place_particle(g1, g2, parameters.set_diam)
        
        # intialize the particle information
        model_particles[particle][0] = p_x
        model_particles[particle][2] = p_y
        model_particles[particle][1] = p_diam
        model_particles[particle][3] = particle # id number for each particle
        
        nulled_vertices[particle] = vertex
        # store each particles vertex information in new_vertices
        new_vertices[particle][0] = (p_x) - (p_diam/2)
        new_vertices[particle][1] = (p_x) + (p_diam/2)
        
    # flatten both vertex arrays so that they can be merged together by insort()
    new_vertices = new_vertices.flatten()
    bed_vertices = bed_vertices.flatten()
        
    valid_vertices = compute_valid_vertices(bed_vertices, new_vertices, nulled_vertices)
    
    return model_particles, chosen_vertex, valid_vertices, nulled_vertices

# from: https://stackoverflow.com/questions/12427146/combine-two-arrays-and-sort
def insort(a, b, kind='mergesort'):
    # took mergesort as it seemed a tiny bit faster for my sorted large array try.
    c = np.concatenate((a, b)) # we still need to do this unfortunatly.
    c.sort(kind=kind)
    flag = np.ones(len(c), dtype=bool)
    np.not_equal(c[1:], c[:-1], out=flag[1:])
    return c[flag]


'''
compute_valid_vertices::method

This method computes the valid vertex set. Funciton is given NumPy arrays of the 
bed vertices, newly introduced vertices and nulled vertices.

Bed vertices is the static set of vertices introduced by the bed formation
established at the start of this model. The set of vertices are those verticle
lines tracing the touching points of the bed particles.

Newly introduced verticed (new_vertices) are those vertices that a newly placed
particle introduces into the vertex calculation. For example, a particle of 
radius 2.5 that is placed at axis 55 will introduce vertices 57.5 and 52.5 as 
two new possible vertices.

Nulled vertices are those vertices which are currently occupied by a particle.
'''
def compute_valid_vertices(bed_vertices, new_vertices, nulled_vertices):
    valid_vertices = copy.deepcopy(bed_vertices)
    
    # nulled set is the shared elements between nulled_vertes and bed_vertices
    nulled_bed_set = set(nulled_vertices)&set(bed_vertices)
    valid_vertices = list(set(bed_vertices)-nulled_bed_set)


    element_count = collections.Counter(new_vertices)
    valid_new_vertices = [item for item in element_count if element_count[item]>1]
    
    nulled_new_set = set(nulled_vertices)&set(valid_new_vertices)
    valid_new_vertices = list(set(valid_new_vertices)-nulled_new_set)

#    valid_vertices = valid_vertices.append(valid_new_vertices)
    valid_vertices = insort(valid_vertices, valid_new_vertices)
    return valid_vertices

# TODO: Need to move n particles (n = e_events) per subregion, not from whole stream
''' 
This needs a description.
'''
def move_model_particles(e_events, event_particles, valid_vertices, model_particles, bed_particles, bed_vertices, nulled_vertices):

    T_p_init1 = (np.random.uniform(parameters.T_pmin, parameters.T_pmax, e_events).
                 reshape(1, e_events))
    # https:\\stackoverflow.com\questions\2106503\
    # Now distribute the random variables over a pseudo exponential
    # distribution based on inverse transform sampling.
    T_p_init2 = np.log(1 - T_p_init1) / (-parameters.lambda_1)
    # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
    # For now I am multuplying by 10 so units are consistent (). Rounding
    # the output to map the entrained particles to the 2D grid.
    L_x_init = np.round((T_p_init2 ** 2) * 10 * 2, 1)
    
    # update the particle information to be original_x + hop distance
    unverified_hop = event_particles[:,0] + L_x_init
    particles_self_vertices = event_particles[:,0] + event_particles[:,1]/2
    reintroduced_vertices = event_particles[:,0]
    # find the closest valid vertex to unverified_hop location 
    verified_hop_placement = np.zeros(e_events)
    
    '''
    This loop will iterate over each proposed hop distance in unverified_hop[]
    and will identify the closest avaliable vertex to the proposed distance.
    Once the closest avaliable vertex has been identified, the nulled_vertices
    array is updated by removing the vertices that have just been vacated (they
    are no longer considered nulled) and then adding the vertices just recently occupied.
    
    If we want to take away a particles self vertices from the simulation then 
    simple uncomment the valid_without_self initialization and replace the 
    enumeration of valid_vertices with valid_without_self. This will force a 
    particle to be pushed to the next avaliable vartex instead of falling back 
    into it's original vartex (i.e as if it was never entrained at all). If this
    way is implemented, one can also comment out the verification loop which
    checks 'if verified_hop == particles_self_vertices[i]' as it is not necessary
    '''
    for i, _ in enumerate(unverified_hop[0]): 
#        valid_without_self = [x for x in valid_vertices if x not in particles_self_vertices]
        try:    
            # from: https://stackoverflow.com/questions/2236906/first-python-list-index-greater-than-x
            # iterate over the valid vertex set and identify the closest vertex equal to or greater than  unverified_hop
            verified_hop = next(x[1] for x in enumerate(valid_vertices) if np.any(x[1] >= unverified_hop[0][i]))
            
            # check that the verified hop distances is not equal to a particles self_vertex (being placed at it's own vertex)
            if verified_hop == particles_self_vertices[i]:
                print("Caught self placement at " + str(particles_self_vertices[i]) + " ... returning particle to original vertex")
                verified_hop = event_particles[i,0]
                
            verified_hop_placement[i] = verified_hop
            print("Particle " + str(int(event_particles[i,3])) + " being entrained from " + str(event_particles[i,0]) + " to " + str(verified_hop))
            # Need to check if particle is being placed between two particles here:
            # update the x-location in rand_particles with the verified_hop
            left_neighbour, right_neighbour = find_neighbours_of(verified_hop, model_particles, bed_particles)
            # get the particles initial x, y and diameter information in the bed
            p_x, p_diam, p_y = place_particle(left_neighbour, right_neighbour, parameters.set_diam)

            try:
                # take nulled vertices without the reintroduced vertex
                nulled_vertices = nulled_vertices[nulled_vertices != reintroduced_vertices[i]] 
                # append the newly occupied vertex to nulled vertices
                nulled_vertices = np.append(nulled_vertices, verified_hop)
            except ValueError:
                print("verified_hop value not in np array")
            except AttributeError:
                print(AttributeError)         
        
        except StopIteration:
            # particle has looped around. Need to add it to the waiting list
            verified_hop_placement[i] = -1
            print("Particle exceeded stream... sending to -1 axis") 


    event_particles[:,0] = verified_hop_placement
    

    
    # now, update model_particle array with new x_location for each of the randomly selected particles
    for idx, particle_id in enumerate(event_particles[:,3]):
        model_particles[int(particle_id)] = event_particles[idx]
    
    # create empty n-2 array to hold the vertices of the model particles
    new_vertices = np.zeros((np.size(model_particles[:,0]),2))
    # for each model particle, store what would be its new vertex values in new_vertices
    for idx, particle in enumerate(model_particles):
        new_vertices[idx][0] = particle[0] - particle[1]/2 # left vertex
        new_vertices[idx][1] = particle[0] + particle[1]/2 # right vertex
        
    new_vertices = new_vertices.flatten()
    bed_vertices = bed_vertices.flatten()
    
    # flatten both vertex arrays so that we can merge them together 
    valid_vertices = compute_valid_vertices(bed_vertices, new_vertices, nulled_vertices)   
   
    ## FOR TESTING:
    # sleep to more easily see the bed migration in matplotlib
    time.sleep(1)
    ###
    
    plot_stream(bed_particles, model_particles, 150, 100/4, valid_vertices)
    return valid_vertices, nulled_vertices





def plot_stream(bed_particles, model_particles, x_lim, y_lim, valid_vertices):
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
    
    radius_array = np.asarray((bed_particles[:,1] / 2.0), dtype=float)
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
    for xxc in valid_vertices:
        plt.axvline(x=xxc, color='m', linestyle='-')
##    for green in chosen_vertex:
#        plt.axvline(x=green, color='g', linestyle='-')
    ### 
    plt.show()
    return
