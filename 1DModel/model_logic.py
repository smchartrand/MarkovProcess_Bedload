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


# FOR TESTING: 
def pause():
    programPause = input("Press the <ENTER> key to continue...")

############################################################################### 
#%% Initial Packing of Streambed 
     
def add_bed_particle(diam, bed_particles, particle_id, pack_idx):
    """ Add 'particle' to the bed particle list.
    
    
    
    Calculates center and elevation of particle 
    from input. Maintains pack_idx and particle_id 
    for next particle iteration.
    
    Keyword arguments:
    diam -- diameter of the pa
    bed_particles -- current list of bed particles
    particle_id -- index into bed_particles list
    pack_idx -- left-most extent of particle
    """
    
    center = pack_idx + (diam/2)

    elevation = 0
    bed_particles[particle_id] = [center, diam, elevation, pack_idx]
  
    # update build parameters
    pack_idx += diam
    particle_id += 1
 
    return particle_id, pack_idx


def build_streambed(bed_particles, diam):
    """ Build the bed particle list until 'packed'.
    
    
    
    Handles calls to add_bed_particle, checks for 
    completness of bed and updates the x-extent
    of stream when the packing exceeds/under packs 
    within 8mm range.
    
    Note: the updates to x-extent are only required 
    when variable particle diameter is being used. 
    
    Returns complete bed_particles list and list of 
    vertices derived from bed_particles 
    """
    running_id = 0
    running_pack_idx = 0
    while True:
        running_id, running_pack_idx = add_bed_particle(diam, 
                                                        bed_particles, 
                                                        running_id, 
                                                        running_pack_idx)
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
    bed_vertices = np.zeros(x_max, dtype=bool)

    bed_vertices[bed_particles[1:,3].astype(int)] = 1
    # x-indexes of avaliable vertices to place model particles at
    bed_vertices = np.transpose(np.nonzero(bed_vertices))

    return bed_particles, bed_vertices


def bed_complete(pack_idx):
    """Check to see if bed is complete based on model params.""" 
    # similarly, if np.count_nonzero(bed_space) == x_max
    if pack_idx >= parameters.x_max:
        return 1
    else: return 0

#%% Bed is built. Place/create n particles in avaliable vertices
##### NOTE: plot_stream SHOULD END UP IN ANOTHER SECTION. NOT APPRO HERE

def determine_num_particles(pack_frac, num_vertices):
    """ Return the number of model particles to be used in model."""
    
#    num_particles = num_vertices * pack_frac
#    num_particles = int(math.ceil(num_particles))
    num_particles = num_vertices - 1

    return num_particles


def place_particle(placement_idx, particle_diam, model_particles, 
                   bed_particles):
    """ Calculate new X and Y of particle based on location in stream.
    
    
    
    Provided a particle's (pA) location (xA) in stream, 
    search for 2 neighbouring particles (n1, n2) that pA might
    come into contact with when placed at xA. 
    
    Calculate the Y and X position of pA
    with a system of equations that uses
    the position of n1 and n2.
    
    Keyword arguments:
    placement_idx -- considered particles locaiton (pA)
    particle_diam -- diameter of considered particle
    model_particles -- model particle list
    bed_particles -- bed particle list
    
    """
    
    #TODO: make sure 'find_neighbours' can never find itself
    # find the two neighbours of the particle
    left_neighbour, right_neighbour = find_neighbours_of(placement_idx, 
                                                         model_particles, 
                                                         bed_particles)
   
    # initialize neighbour (g1 and g2) and placed particle information using info from the storage arrays
    x1 = left_neighbour[0]
    y1 = left_neighbour[2]
    r1 = left_neighbour[1] / 2
    
    x2 = right_neighbour[0]
    y2 = right_neighbour[2]
    r2 = right_neighbour[1] / 2
    
    rp = particle_diam / 2 
    
    #print(f'Calculating placement using elevations {y1} and {y2}')
    
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
    
 
    return p_x, p_y


def update_particle_states(model_particles, bed_particles):
    """ Set each model particle state according 'active' rules.
    
    
    
    If any model particle (pX) has a particle 
    resting on it in the stream then pX must 
    be set to Inactive indicated by a boolean 0 
    in its model particle list slot.
    
    If pX does not have any particles resting
    on it then it is considered Active 
    indicated by a boolean 1.
    
    Note: bed particles are always considered
    Inactive in the model.
    
    
    Keyword arguments:ß
    model_particles -- model particle list
    bed_particles -- bed particle list
    
    """
    
    for particle in model_particles:
        # set particle state to active (1)
        set_state(particle, 1)
        
    in_stream_particles = model_particles[model_particles[:,0] != -1]
    for particle in in_stream_particles:
        
        left_neighbour, right_neighbour = find_neighbours_of(particle[0], 
                                                             model_particles, 
                                                             bed_particles)
        # note: this method below could be improved if find_neighbours_of 
        # would indicate if a neighbour belongs to the model or bed particles
        if left_neighbour[2] > 0 and left_neighbour[2] < particle[2]:
            lmodel_index = np.where(model_particles[:,3] == left_neighbour[3])
            model_particles[lmodel_index] = set_state(model_particles[lmodel_index][0], 0)
        else:
            continue
                 
        if right_neighbour[2] > 0 and right_neighbour[2] < particle[2]:
            rmodel_index = np.where(model_particles[:,3] == right_neighbour[3])
            model_particles[rmodel_index] = set_state(model_particles[rmodel_index][0], 0)
        else:
            continue
    
    return model_particles
        

def set_state(particle, status):
    """ Set particle state to desired status."""
    particle[4] = status
    #print(f'Particle {particle[3]} set to {status} (1=active, 0=inactive)')
    return particle

    
def find_neighbours_of(particle_idx, model_particles, bed_particles):
    """ Find the 2 neighbouring (left, right) particles of particle_idx.
    
    
    A particle's center must be within +/- radius of
    particle_idx to be considered a neighbour of 
    particle_idx. Of those particles that qualify, only
    the two with the greatest elevations (Y) will be taken.
    
    Keyword arguments:
        
    particle_idx -- the x-location considered for neighbours
    model_particles -- model particle list
    bed_particles -- bed particle list
    """ 

    # look for particles with center R (radius) away from particle_idx 
    left_neighbour_center = particle_idx - (parameters.set_radius)
    right_neighbour_center = particle_idx + (parameters.set_radius)
    
    #TODO: assert/check that only one neighbour has been found
    #     If multiple neighbours found this means a verticle stack has occured?
    
    # Look for left neighbour:
    # first search model particles for matches
    left_neighbour_idx_model = np.where(model_particles[:,0] == left_neighbour_center)
    
    # if no such model particle found, search bed particle array
    if left_neighbour_idx_model[0].size == 0: 
        left_neighbour_idx_bed = np.where(bed_particles[:, 0] == left_neighbour_center)
        if left_neighbour_idx_bed[0].size == 0:
            # TODO: need to handle these errors - exit simulation?
            raise ValueError(f'No left neighbours found for location {particle_idx}')
        else:
            left_neighbour = bed_particles[left_neighbour_idx_bed]
    else:
        left_neighbour = model_particles[left_neighbour_idx_model]

        
    # Look for right neighbour:
    # first search model particles for matches
    right_neighbour_idx_model = np.where(model_particles[:,0] == right_neighbour_center)
    
    if right_neighbour_idx_model[0].size == 0:
        right_neighbour_idx_bed = np.where(bed_particles[:,0] == right_neighbour_center)
        if right_neighbour_idx_bed[0].size == 0:
            raise ValueError(f'No right neighbours found for location {particle_idx}')
        else:
            right_neighbour = bed_particles[right_neighbour_idx_bed]
    else:
        right_neighbour = model_particles[right_neighbour_idx_model]  
        
    # print(f'Event particle landing at {particle_idx} has l-neighbour {left_neighbour[0][0]} and r-neighbour {right_neighbour[0][0]}')

    return left_neighbour[0], right_neighbour[0]

##TODO: extract into create and set functions
def set_model_particles(bed_vertices, bed_particles):
    """ Create model particle list and situate in model stream.
    
    
    
    Create list of n model particles based 
    the packing fraction.
    
    Randomly assign avaliable x-vertices 
    to each model particle. Avaliable
    vertices are derived from the list of
    bed particles. 
    
    The structure of a resulting particle:
        [0] = center coordinate,
        [1] = diameter,
        [2] = elevation,
        [3] = uid,
        [4] = active (boolean)
    
    
    Keyword arguments:
    bed_vertices -- list of vertices based on bed particles
    bed_particles -- bed particle list
    
     """

    # create a boolean area representing the avaliability of the vertices
    num_vertices = np.size(bed_vertices)
    already_selected = [False] * num_vertices
    
    # determine the number of model particles that should be introduced into the stream bed
    num_particles = determine_num_particles(parameters.Pack, num_vertices)
    
    ### FOR TESTING - REMOVE LINE BELOW AND UNCOMMENT ABOVE LINE BEFORE SUBMITTING
    # num_particles = np.size(bed_vertices) - 1
    
    # create an empty n-4 array to store model particle information
    model_particles = np.zeros([num_particles, 5], dtype='float')
    
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
  
        # place particle at the chosen vertex
        p_x, p_y = place_particle(vertex, parameters.set_diam, model_particles, bed_particles)
        
        # intialize the particle information
        model_particles[particle][0] = p_x
        model_particles[particle][1] = parameters.set_diam
        model_particles[particle][2] = p_y
        model_particles[particle][3] = particle # id number for each particle
        model_particles[particle][4] = 1 # each particle begins as active
        
        nulled_vertices[particle] = vertex
        # store each particles vertex information in new_vertices
        new_vertices[particle][0] = (p_x) - (parameters.set_diam/2)
        new_vertices[particle][1] = (p_x) + (parameters.set_diam/2)
        
    # flatten both vertex arrays so that they can be merged together by insort()
    new_vertices = new_vertices.flatten()
    bed_vertices = bed_vertices.flatten()
        
    available_vertices = compute_available_vertices(bed_vertices, model_particles, nulled_vertices)
    
    # update particle states so that supporting particles are inactive
    model_particles = update_particle_states(model_particles, bed_particles)
    
    return model_particles, chosen_vertex, available_vertices, nulled_vertices


# from: https://stackoverflow.com/questions/12427146/combine-two-arrays-and-sort
def insort(a, b, kind='mergesort'):
    # took mergesort as it seemed a tiny bit faster for my sorted large array try.
    c = np.concatenate((a, b)) # we still need to do this unfortunatly.
    c.sort(kind=kind)
    flag = np.ones(len(c), dtype=bool)
    np.not_equal(c[1:], c[:-1], out=flag[1:])
    return c[flag]


def compute_available_vertices(bed_vertices, model_particles, nulled_vertices):
    """ Compute the avaliable vertices in the model stream.

    
    Keyword arguments:
    
    bed_vertices -- list of vertices introduced by bed particles
    side_vertices -- list of vertices introduced by the L/R 
                        extent of model particles
    nulled_vertices -- list of vertices that have been occupied    

    """
    #TODO: .... just remove the nulled vertex
    
    cp_model_particles = model_particles.copy()
    no0_model_particles = cp_model_particles[cp_model_particles[:,0] != 0]
    fn_model_particles = no0_model_particles[no0_model_particles[:,0] != -1]
            
    # create empty n-2 array to hold the vertices of the model particles
    new_vertices = np.zeros((np.size(fn_model_particles[:,0]),2))
            
    for idx, particle in enumerate(fn_model_particles):
        new_vertices[idx][0] = particle[0] - parameters.set_radius # left vertex
        new_vertices[idx][1] = particle[0] + parameters.set_radius # right vertex
    
    new_vertices = new_vertices.flatten()
    bed_vertices = bed_vertices.flatten()

    

    # grab the shared elements between nulled_vertices and bed_vertices
    nulled_bed_set = set(nulled_vertices)&set(bed_vertices)
    available_bed_vertices = list(set(bed_vertices)-nulled_bed_set)

    element_count = collections.Counter(new_vertices)
    avaliable_by_touch = [item for item in element_count if element_count[item]>1]
    
#    total_nulled = nulled_bed_set.union(set(nulled_vertices))
#    count = collections.Counter(total_nulled)
#    for vertex in avaliable_by_touch:
#        occr = count[float(vertex)]
#        if occr != 0 and (occr % 2) == 0:
#            #print(f'Vertex {vertex} occured even number of times in nulled_vertices')
#            nulled_vertices = nulled_vertices[nulled_vertices != float(vertex)]
#            
#    for vertex in bed_vertices:
#        occr = count[float(vertex)]
#        if occr != 0 and (occr % 2) != 0:
#            #print(f'Vertex {vertex} occured odd number of times in nulled_vertices')
#            nulled_vertices = nulled_vertices[nulled_vertices != float(vertex)]    
    nulled_new_set = set(nulled_vertices)&set(avaliable_by_touch)

    valid_new_vertices = list((set(avaliable_by_touch)-nulled_new_set))
    
    

    available_vertices = insort(available_bed_vertices, valid_new_vertices)
    
    return available_vertices


def lift_event_particles(idx_of_entrainment, model_particles, bed_vertices, 
                         nulled_vertices):
    """ 
    Lift all particles selected for entrainment and recalculate vertices.
    """
    print(f'nulled before lift: {nulled_vertices}')
    # release all x locations of the entrainment particles
    for entrainment in idx_of_entrainment:
        particle = model_particles[entrainment]
        nulled_vertices = nulled_vertices[nulled_vertices != particle[0]]
    
    print(f'Nulled after lift: {nulled_vertices}')
    # calculate avaliable vertices as if event particles don't exist
    avaliable_vertices = compute_available_vertices(bed_vertices, 
                                                    model_particles, 
                                                    nulled_vertices)
        
    return avaliable_vertices, bed_vertices, nulled_vertices

# TODO: Need to move n particles (n = e_events) per subregion, not from whole stream
def move_model_particles(e_events, idx_of_entrainment, model_particles, 
                         bed_particles, bed_vertices, nulled_vertices):
    """ Move selected model particles for entrainment.
    
    
    
    Calculate jump distances for each particle and place 
    where stream vertices allow. 
    
    Follows the current logic:
        1. lift all event particles and recalculate vertices
        2. move all event particles to downstream vertex 
            nearest to it's unverified hop (x + L_x_init)
        3. Check if any event particles chose same vertex
            a. if yes, randomly select from vertices for
                n - (n-1) to be pushed forward
    
    Keyword arguments:
    e_events -- the number of entrainment events to occur
    event_particles -- the particles selected for entrainment
    avaliable_vertices -- list of avaliable vertices
    model_particles -- list of model particles
    bed_particles -- list of bed particles
    bed_vertices -- list of bed vertices
    nulled_vertices -- list of currently nulled vertices
    """

    T_p_init1 = (np.random.uniform(parameters.T_pmin, parameters.T_pmax, 
                                   e_events).reshape(1, e_events))
    # https:\\stackoverflow.com\questions\2106503\
    # Now distribute the random variables over a pseudo exponential
    # distribution based on inverse transform sampling.
    T_p_init2 = np.log(1 - T_p_init1) / (-parameters.lambda_1)
    # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
    # For now I am multuplying by 10 so units are consistent (). Rounding
    # the output to map the entrained particles to the 2D grid.
    L_x_init = np.round((T_p_init2 ** 2) * 10 * 2, 1)
    
    available_vertices, bed_vertices, nulled_vertices = lift_event_particles(
            idx_of_entrainment, model_particles, bed_vertices, nulled_vertices)
    
    
# TODO: Update comment for this section    
    for entrainment in idx_of_entrainment: 
        
        particle = model_particles[entrainment]
        unverified_hop = particle[0] + L_x_init 
        try:       
            # from: https://stackoverflow.com/questions/2236906/first-python-list-index-greater-than-x
            # iterate over the valid vertex set and identify the closest vertex equal to or greater than  unverified_hop
            verified_hop = next(x[1] for x in enumerate(available_vertices) if np.any(x[1] >= unverified_hop))
                
            print("Particle " + str(int(particle[3])) + " being entrained from " + str(particle[0]) + " to " + str(verified_hop))

            # get the particles initial x, y and diameter information in the bed
            new_x, new_y = place_particle(verified_hop, parameters.set_diam, model_particles, bed_particles)
            particle[0] = new_x
            particle[2] = new_y

        except StopIteration:
            print("Particle exceeded stream... sending to -1 axis") 
            new_x = -1
            particle[0] = new_x

    model_particles[entrainment] = particle

    # update particle states so that supporting particles are inactive
    model_particles = update_particle_states(model_particles, bed_particles)
    
    ## FOR TESTING:
    # sleep to more easily see the bed migration in matplotlib
    time.sleep(1)
    ###
    
    plot_stream(bed_particles, model_particles, 250, 100/4, available_vertices)
    return available_vertices, nulled_vertices


# Taken from original model
def plot_stream(bed_particles, model_particles, x_lim, y_lim, available_vertices):
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
    x_center = bed_particles[:,0]
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
    for xxc in available_vertices:
        plt.axvline(x=xxc, color='m', linestyle='-')
##    for green in chosen_vertex:
#        plt.axvline(x=green, color='g', linestyle='-')
    ### 
    plt.show()
    return
