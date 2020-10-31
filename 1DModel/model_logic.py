from __future__ import division
import math
import random
import time
import numpy as np
import sympy as sy
import copy
import parameters # Import the parameters defined in the parameters file
from collections import Counter

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle



# FOR TESTING: 
def pause():
    programPause = input("Press the <ENTER> key to continue...")
    
    
class Model():
    
    ITERATION_TEMPLATE = ('''\n
    --------------------- Iteration {step} ---------------------\n
    # of entrainment events: {e_events}\n
    Particles to be entrained: {particles}\n                          
    ''')
    
    def __init__(self, bed_particles, model_particles, e_events_store):
        self.bed_particles = bed_particles
        self.model_particles = model_particles
        self.e_events_store = e_events_store
    
    def run(self):        
        for step in range(parameters.n_iterations):
            # Calculate the number of entrainment events per-unit-area
            e_events = np.random.poisson(parameters.lambda_1, None)
            if e_events == 0:
                e_events = 1 #???
                  
            # e_events = e_events * bed_sampreg
            in_stream_particles = self.model_particles[
                                    self.model_particles[:,0] != -1]
            active_particles =  in_stream_particles[
                                    in_stream_particles[:,4] != 0]
            
            id_of_event_particles = []
            if e_events > len(active_particles):
                e_events = len(active_particles)
            random_sample = random.sample(range(len(active_particles)), 
                                          e_events)
            for index in random_sample:
                id_of_event_particles.append(int(active_particles[index][3])  )
            
            ii = np.where(self.model_particles[:,0] == -1)[0]
            for index in ii:
                # send particle to 0 (starting point)
                self.model_particles[index][0] = 0 
                id_of_event_particles.append(index)
                e_events = e_events + 1
            
            num_event_particles = len(id_of_event_particles) 
            if e_events != num_event_particles:
                # 1. log that the # entrainment events has been altered
                
                # 2. update e_events
                e_events = num_event_particles
            
            available_vertices = compute_available_vertices(self.model_particles, 
                                                            self.bed_particles)
            print(self.ITERATION_TEMPLATE.format(
                                        step=step, 
                                        e_events=e_events, 
                                        particles=id_of_event_particles))
            self.e_events_store[step] = e_events
            
            available_vertices = move_model_particles(e_events, 
                                                      id_of_event_particles, 
                                                      self.model_particles, 
                                                      self.bed_particles)
            

class NoSupportingParticlesFoundError(Exception):
    """ Exception raised when insufficient supporting particles found
        for a given particle
        
    
    Keyword arguments:
    particle - the particle that caused the error
    message - explanation of error
    
    """
    #TODO: have exception display the supporting particles (if one was found)
    def __init__(self, particle, message = """Zero or one supporting particles 
                                                      found for particle"""):
        self.particle = particle
        self.supporting = supporting_particles
        self.message = message
    
    def __str__(self):
        return f' {self.message}: {self.particle}'

############################################################################### 
#%% Initial Packing of Streambed 
     
def add_bed_particle(diam, bed_particles, particle_id, pack_idx):
    """ Add 'particle' to the bed particle list.
    
    
    
    Calculates center and elevation of particle 
    from input. Maintains pack_idx and particle_id 
    for next particle iteration.
    
    Builds particle of the following structure:
        [0] = center coordinate,
        [1] = diameter,
        [2] = elevation,
        [3] = uid,
        [4] = active (boolean)
    
    Keyword arguments:
    diam -- diameter of the pa
    bed_particles -- current list of bed particles
    particle_id -- index into bed_particles list
    pack_idx -- left-most extent of particle
    """
    
    center = pack_idx + (diam/2)
    state = 0

    elevation = 0
    bed_particles[particle_id] = [center, diam, elevation, pack_idx, state]
  
    # update build parameters
    pack_idx += diam
    particle_id += 1
 
    return particle_id, pack_idx


def build_streambed():
    """ Build the bed particle list.
    
    
    Handles calls to add_bed_particle, checks for 
    completness of bed and updates the x-extent
    of stream when the packing exceeds/under packs 
    within 8mm range.
    
    Note: the updates to x-extent are only required 
    when variable particle diameter is being used. 
    
    Return values:
    bed_particles -- list of bed particles
    bed_vertices -- list of available vertices 
                    based on bed list 
    """
    bed_particles = np.zeros([parameters.max_particles, 5],dtype=float)
    
    running_id = 0
    running_pack_idx = 0
    while True:
        running_id, running_pack_idx = add_bed_particle(parameters.set_diam, 
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


def place_particle(particle, particle_diam, model_particles, 
                   bed_particles):
    """ Calculate new X and Y of particle based on location in stream.
    
    
    
    Provided a particle's (pA) location (xA) in stream, 
    search for 2 neighbouring particles (n1, n2) that pA might
    come into contact with when placed at xA. 
    
    Calculate the Y and X 'resting' position of pA
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
    left_support, right_support = find_supporting_particles_of(particle, 
                                                        model_particles, 
                                                        bed_particles,
                                                        already_placed=False)
   
    # initialize neighbour (g1 and g2) and placed particle information using info from the storage arrays
    x1 = left_support[0]
    y1 = left_support[2]
    r1 = left_support[1] / 2
    
    x2 = right_support[0]
    y2 = right_support[2]
    r2 = right_support[1] / 2
    
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
    """ Set each model particle's current 'active' state.
    
    
    
    If any model particle (pX) has a particle 
    resting on it in the stream then pX must 
    be set to Inactive indicated by a boolean 0 
    in its model particle list slot.
    
    If pX does not have any particles resting
    on it then it is considered Active 
    indicated by a boolean 1.
    
    Note: bed particles are always considered
    Inactive.
    
    
    Keyword arguments:
    model_particles -- model particle list
    bed_particles -- bed particle list
    
    """
    
    for particle in model_particles:
        #TODO: enumerate this somehow
        set_state(particle, 1)
        
    in_stream_particles = model_particles[model_particles[:,0] != -1]
    for particle in in_stream_particles:
        
        left_neighbour, right_neighbour = find_supporting_particles_of(
                                                        particle, 
                                                        model_particles, 
                                                        bed_particles,
                                                        already_placed=True)
            
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

    
def find_supporting_particles_of(particle, model_particles, bed_particles,
                                 already_placed):
    """ Find the 2 supporting particles for a given particle 'particle'.
    
    #TODO: Update description
    A particle's center must be within +/- radius of
    particle_idx to be considered a neighbour of 
    particle_idx. Of those particles that qualify, only
    the two with the greatest elevations (Y) will be taken.
    
    Keyword arguments:   
    particle_idx -- the x-location considered for neighbours
    model_particles -- model particle list
    bed_particles -- bed particle list
    """ 
    
    if already_placed:
        # TODO: we know that a particle cannot be supported by a particle
        # not directly below it. Find a programmic way to identify layers and 
        # eleminate all particles not in the layer directly below a particle  
        considered_particles = model_particles[
                                    (model_particles[:,2] < particle[2])]
        all_particles = np.concatenate((considered_particles, 
                                               bed_particles), axis=0)
    else:
        all_particles = np.concatenate((model_particles, 
                                               bed_particles), axis=0)
        
    # look for particles with center R (radius) away from particle_idx 
    left_center = particle[0] - (parameters.set_radius)
    right_center = particle[0] + (parameters.set_radius)
   
    left_candidates = all_particles[all_particles[:,0]== left_center]
    left_support = left_candidates[left_candidates[:,2] 
                                   == np.max(left_candidates[:,2])]
    
    right_candidates = all_particles[all_particles[:,0] == right_center]
    right_support = right_candidates[right_candidates[:,2]
                                    == np.max(right_candidates[:,2])]
    # TODO: .size == 0 considered non-pythonic
    if right_support.size == 0 or left_support.size == 0:
        raise NoSupportingParticlesFoundError(particle)

    return left_support[0], right_support[0]

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
        
        # intialize the particle information
        model_particles[particle][0] = vertex
        model_particles[particle][1] = parameters.set_diam
        
        model_particles[particle][3] = particle # id number for each particle
        model_particles[particle][4] = 1 # each particle begins as active
        
        # place particle at the chosen vertex
        p_x, p_y = place_particle(model_particles[particle], parameters.set_diam, model_particles, bed_particles)
        
        model_particles[particle][0] = p_x
        model_particles[particle][2] = p_y
 
       
    available_vertices = compute_available_vertices(model_particles, bed_particles)
    
    # update particle states so that supporting particles are inactive
    model_particles = update_particle_states(model_particles, bed_particles)
    
    return model_particles, available_vertices


# from: https://stackoverflow.com/questions/12427146/combine-two-arrays-and-sort
def insort(a, b, kind='mergesort'):
    # took mergesort as it seemed a tiny bit faster for my sorted large array try.
    c = np.concatenate((a, b)) # we still need to do this unfortunatly.
    c.sort(kind=kind)
    flag = np.ones(len(c), dtype=bool)
    np.not_equal(c[1:], c[:-1], out=flag[1:])
    return c[flag]


def compute_available_vertices(model_particles, bed_particles, lifted=False,
                               lifted_particles=None):
    """ Compute the avaliable vertices in the model stream.

    #TODO: update desc and keyword to outline the logic of computation
    
    Keyword arguments: 
    model_particles    -- list of model particles
    bed_particles      -- list of bed particles
    lifted             -- boolean flag indicating if calculation should 
                            take into account particles being 'lifed' 
    lifted_particles   -- idx of the 'lifted' particles
    """
    nulled_vertices = []
    avail_vertices = []
    
    if lifted == True:
        model_particles_lifted = copy.deepcopy(model_particles)
        
        model_particles_lifted = np.delete(model_particles_lifted, 
                                           lifted_particles, 0)
        
        all_particles = np.concatenate((model_particles_lifted, 
                                        bed_particles), axis=0)
    else:    
        all_particles = np.concatenate((model_particles, 
                                        bed_particles), axis=0)
    
    elevations = np.unique(all_particles[:,2])
    # sort elevation in descending order
    # https://stackoverflow.com/questions/26984414/
    elevations[::-1].sort()
    
    for elevation in elevations:
        tmp_particles = all_particles[all_particles[:,2] == elevation]
        
        for particle in tmp_particles:    
            nulled_vertices.append(particle[0])
        
        right_vertices = tmp_particles[:,0] + parameters.set_radius
        left_vertices = tmp_particles[:,0] - parameters.set_radius
        tmp_shared_vertices = np.intersect1d(left_vertices, right_vertices)
        
        for vertex in tmp_shared_vertices:
            if vertex not in nulled_vertices:
                avail_vertices.append(vertex)
                
        del(tmp_shared_vertices)
        del(tmp_particles)
        
    available_vertices = np.array(avail_vertices)
    
    return available_vertices


# TODO: Need to move n particles (n = e_events) per subregion, not from whole stream
def move_model_particles(e_events, idx_of_event_particles, model_particles, 
                         bed_particles):
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
    e_events                -- the number of entrainment events to occur
    idx_of_event_particles  -- list of indexes into the model_particles array
                              of the entrainment particles
    model_particles         -- list of model particles
    bed_particles           -- list of bed particles
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
    L_x_init = list(L_x_init)
    
    # compute the avaliable vertices with event particles 'lifted'
    available_vertices = compute_available_vertices(model_particles,
                                                     bed_particles, 
                                                     True,
                                                     idx_of_event_particles)
    available_vertices = np.sort(available_vertices)
    av_copy = copy.deepcopy(available_vertices)
    for count, entrainment in enumerate(idx_of_event_particles): 
        particle = model_particles[entrainment]
        unverified_hop = particle[0] + L_x_init[0][count]    
        # from: https://stackoverflow.com/questions/2236906/first-python-list-index-greater-than-x
        # iterate over the valid vertex set and identify the closest vertex equal to or greater than  unverified_hop

        # available vertices _needs_ to be sorted for 575-593 to work
        # TODO: extract 575-593 into function
        forward_vertices = available_vertices[available_vertices >= unverified_hop]
        if forward_vertices.size < 1:
            print("Particle exceeded stream... sending to -1 axis") 
            new_x = -1
            particle[0] = new_x
            
        else:
            verified_hop = forward_vertices[0]
            
            print(f"""Particle {particle[3]} being entrained from {particle[0]} to {verified_hop} with desired hop {unverified_hop}""")
            av_copy = av_copy[av_copy != verified_hop]
            particle[0] = verified_hop
            # get the particles initial x, y and diameter information in the bed
            new_x, new_y = place_particle(particle, parameters.set_diam, model_particles, bed_particles)
            particle[0] = new_x
            particle[2] = new_y           

        model_particles[entrainment] = particle
        #print(f'model_particles after single entrainment: {model_particles}')
    
    #TODO: improve implementation
    unique_entrainment_choices = False
    while unique_entrainment_choices == False:
        event_particles = model_particles[idx_of_event_particles]
        # Note: np.unique() runs in O(nlogn) due to sorting implementation.
        # if model consistently uses _very_ large particle sets then 
        # Counter() duplicate method should be considered given RT O(n)
        unique, count = np.unique(event_particles[:,0], 
                                           return_counts=True)
        nonunique_choices = unique[count < 1]

        samev_particles = event_particles[np.in1d(event_particles[:,0], nonunique_choices)]
        
        if samev_particles.size == 0:
            unique_entrainment_choices = True
            break
        else:
            stay_particle = samev_particles[random.sample(range(len(samev_particles)), 1)]
            relocated_particles = samev_particles[samev_particles[:,3] != stay_particle[0][3]]
            for particle in relocated_particles:
                forward_vertices = av_copy[av_cpy > particle[0]]
                if forward_vertices.size < 1:
                    print("Particle exceeded stream... sending to -1 axis during Unique Choice correction") 
                    new_x = -1
                    particle[0] = new_x
                    
                else:
                    verified_hop = forward_vertices[0]
                    particle[0] = verified_hop
                    # get the particles initial x, y and diameter information in the bed
                    new_x, new_y = place_particle(particle, parameters.set_diam, model_particles, bed_particles)
                    particle[0] = new_x
                    particle[2] = new_y 
                    
            model_particles[entrainment] = particle
                
                
        
    # update particle states so that supporting particles are inactive
    model_particles = update_particle_states(model_particles, bed_particles)
    
    ## FOR TESTING:
    # sleep to more easily see the bed migration in matplotlib
    time.sleep(1)
    ###
    
    available_vertices = compute_available_vertices(model_particles, bed_particles)
    plot_stream(bed_particles, model_particles, 250, 100/4, available_vertices)
    return available_vertices


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
