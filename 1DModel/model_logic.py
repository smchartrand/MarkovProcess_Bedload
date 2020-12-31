from __future__ import division
import math
import random
import numpy as np
import sympy as sy
import copy
import parameters # Import the parameters defined in the parameters file

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from collections import defaultdict
from matplotlib.patches import Circle


# FOR TESTING: 
def pause():
    programPause = input("Press the <ENTER> key to continue...")
    
# TODO: refactor from class to simple struct (dict, maybe?)
class Subregion():
    """ Subregion class.
    
    Each instance of Subregion contains
    the name, left and right boundaries 
    of a subregion. 
    
    Name and boundaries are set during 
    instantiation and can be retrieved
    afterwards using helper methods.
    """
    def __init__(self, name, left_boundary, right_boundary):
        self.name = name
        self.left_boundary = left_boundary
        self.right_boundary = right_boundary
        
    def leftBoundary(self):
        return self.left_boundary
    
    def rightBoundary(self):
        return self.right_boundary
    
    def getName(self):
        return self.name
        
    
def get_event_particles(e_events, subregions, model_particles):
    """ Find and return list of particles to be entrained

    Keyword arguments:
    e_events -- Number of events requested per subregion 
    subregions -- List of Subregion objects
    model_particles -- The Model's model_particles array 

    Returns:
    total_e_events -- Number of events over entire stream
    event_particles -- List of particles to be entrained

    """
    if e_events == 0:
        e_events = 1 #???
    
    event_particles = []
    for subregion in subregions:
        # Filter array for only active, in-stream particles per subregion
        subregion_particles = model_particles[
                (model_particles[:,0] >= subregion.leftBoundary())
              & (model_particles[:,0] <= subregion.rightBoundary())]
        in_stream_particles = subregion_particles[
                                                subregion_particles[:,0] != -1]
        active_particles =  in_stream_particles[
                                                in_stream_particles[:,4] != 0]
        
        if e_events > len(active_particles):
            random_sample = random.sample(range(len(active_particles)), 
                                         len(active_particles))
        else:
            random_sample = random.sample(range(len(active_particles)), 
                                         e_events)
        subregion_event_ids = []  
        for index in random_sample:
            #NOTE: this only works because index=id in the model_particle array
            subregion_event_ids.append(int(active_particles[index][3])  )
        
        ghost_particles = np.where(model_particles[:,0] == -1)[0]
        for index in ghost_particles:
            model_particles[index][0] = 0 
            subregion_event_ids.append(index)
        
        if e_events != len(subregion_event_ids):
            msg = (
                     f'\nRequested {e_events} events in {subregion.getName()} ' 
                     f'but {len(subregion_event_ids)} occuring\n'
            )
            print(msg)

        event_particles = event_particles + subregion_event_ids
        
    return event_particles


def define_subregions(bed_length, num_subregions):
    """ Define subregions list.
    

    Keyword arguments:
    bed_length -- The length of the model bed.
    subregions -- The number of subregions to create.

    Returns:
    subregions_arr -- The np array of subregions

    """
    assert(math.remainder(bed_length, num_subregions) == 0)
    
    subregion_length = bed_length/num_subregions
    left_boundary = 0
    subregions_arr = []
    for region in range(num_subregions):       
        right_boundary = left_boundary + subregion_length
        subregion = Subregion(f'subregion_{region}', left_boundary, right_boundary)
        left_boundary = right_boundary
        
        subregions_arr.append(subregion)
    
    return subregions_arr

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
    
    # bed can be packed +/- from the default x_max depending on the 
    # packing pattern -- therefore update x_max once bed is complete to new +- 8mm size
    x_max = int(math.ceil(bed_particles[running_id-1][1]+ bed_particles[running_id-1][3]))
    # strip zero element particles tuples
    valid = ((bed_particles==0).all(axis=(1)))
    bed_particles = bed_particles[~valid]

    return bed_particles, x_max


def get_bed_vertices(bed_particles):
    bed_vertices = np.zeros(parameters.x_max, dtype=bool)
    bed_vertices[bed_particles[1:,3].astype(int)] = 1
    # x-indexes of avaliable vertices to place model particles at
    bed_vertices = np.transpose(np.nonzero(bed_vertices))
    
    return bed_vertices


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
    
    num_particles = num_vertices * pack_frac
    num_particles = int(math.ceil(num_particles))
    # num_particles = num_vertices - 1

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
   
    # TODO: make this prettier/more readable
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
    
    np.format_float_positional(p_x, precision=2)
 
    return p_x, p_y


def update_particle_states(model_particles, bed_particles):
    """ Set each model particle's current 'active' state.
    
    
    
    If any model particle (pX) has a particle 
    resting on it in the stream then pX must 
    be set to Inactive indicated by a boolean 0.
    
    If pX does not have any particles resting
    on top of it then it is considered Active 
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
    """ Find the 2 supporting particles for a given particle.
    
    Provided a particle (numpy array selection), this
    function will search the stream for particles
    that could be considered 'supporting' particles.
    
    A supporting particle is a particle that is placed
    directly below our particle of concern. Supporting
    particles are those particles which 'hold up' the
    particle of concern.
    
    This function can search for supporting particles
    in two scenarios:
        1. The particles is already placed/set in 
        the stream.
        2. The particles is looking to be placed 
        in the stream at it's defined x-location (i.e
        after an entrainment event)
    
    Searching for supporting particles at location x could 
    result in two different results depending on the 
    aforementioned scenario, hence the distinct methods.
    
    Keyword arguments:   
    particle -- array representing a particle 
    model_particles -- model particle list
    bed_particles -- bed particle list
    already_placed -- boolean flag indicating if particle has
                      already been placed in the stream, or
                      is looking to be placed
    
    Returns:
    left_support -- the left supporting particle
    right_support -- the right supporting particle
    """ 
    bad_search = False
    if already_placed: # only consider particle below current elevation
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
   
    left_candidates = all_particles[all_particles[:,0] == left_center]
    try:
        left_support = left_candidates[left_candidates[:,2] 
                                       == np.max(left_candidates[:,2])]
    except ValueError:
        bad_search = True
        print(f'\n\nERROR: no left supporting particles found at {left_center},'
              f'searching for an article at {particle[0]}\n\n')
    
    right_candidates = all_particles[all_particles[:,0] == right_center]
    try:
        right_support = right_candidates[right_candidates[:,2]
                                    == np.max(right_candidates[:,2])]
    except ValueError:
        bad_search = True
        print(f'\n\nERROR: No right supporting particles found at {right_center},'
              f'searching for an article at {particle[0]}\n\n')
    if bad_search:
        raise ValueError(f'Supporting particles for particle {particle[3]} not found in model_particles')

    return left_support[0], right_support[0]

##TODO: extract into create and set functions
def set_model_particles(bed_particles):
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

    bed_vertices = get_bed_vertices(bed_particles)
    # create a boolean area representing the avaliability of the vertices
    num_vertices = np.size(bed_vertices)
    already_selected = [False] * num_vertices
    
    # determine the number of model particles that should be introduced into the stream bed
    num_particles = determine_num_particles(parameters.Pack, num_vertices)
    # create an empty n-5 array to store model particle information
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
        p_x, p_y = place_particle(model_particles[particle], 
                                  parameters.set_diam, 
                                  model_particles, 
                                  bed_particles)
        
        model_particles[particle][0] = p_x
        model_particles[particle][2] = p_y
 

    
    # update particle states so that supporting particles are inactive
    model_particles = update_particle_states(model_particles, bed_particles)
    
    return model_particles


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
                          calculate using 'lifted' vertice. Default False
    lifted_particles   -- idx of the 'lifted' particles. Default None
    """
    nulled_vertices = []
    avail_vertices = []
    
    # If we are lifting particles, we need to consider the subset of particles
    # that includes every particles _except_ the particles being 
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

def run_entrainments(model_particles, bed_particles, event_particle_ids, lambda_1):
    # compute available vertices for this iteration, lifting event particles
    available_vertices = compute_available_vertices(model_particles, bed_particles, lifted=True, lifted_particles=event_particle_ids)
    unverified_entrainments = fathel_furbish_hops(event_particle_ids, model_particles, lambda_1)
    entrainment_dict, model_particles, available_vertices = move_model_particles(unverified_entrainments, model_particles, bed_particles, available_vertices)
    unique_entrainments, redo_particle_ids = check_unique_entrainments(entrainment_dict)
     
    while not unique_entrainments:
        redo_entrainments = model_particles[np.searchsorted(model_particles[:,3], redo_particle_ids)]
        entrainment_dict, model_particles, available_vertices = move_model_particles(redo_entrainments, model_particles, bed_particles, available_vertices)
        unique_entrainments, redo_particles = check_unique_entrainments(entrainment_dict)

    model_particles = update_particle_states(model_particles, bed_particles)
    
    return model_particles
  
        
def fathel_furbish_hops(event_particle_ids, model_particles, lambda_1):
    """ Given a list of particles
    """
    event_particles = model_particles[np.searchsorted(model_particles[:,3], event_particle_ids)]
    
    T_p_init1 = (np.random.uniform(parameters.T_pmin, parameters.T_pmax, 
                                   len(event_particle_ids)).reshape(1, len(event_particle_ids)))
    # https:\\stackoverflow.com\questions\2106503\
    # Now distribute the random variables over a pseudo exponential
    # distribution based on inverse transform sampling.
    T_p_init2 = np.log(1 - T_p_init1) / (-lambda_1)
    # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
    # For now I am multuplying by 10 so units are consistent (). Rounding
    # the output to map the entrained particles to the 2D grid.
    L_x_init = np.round((T_p_init2 ** 2) * 10 * 2, 1)
    L_x_init = list(L_x_init)
    
    event_particles[:,0] = event_particles[:,0] + L_x_init
    
    return event_particles
 
       
def move_model_particles(event_particles, model_particles, bed_particles, available_vertices):
    """ Given an array of particles, move them to next closest valid vertex
    within the model stream.

    Keyword arguments:
        event_particles -- list of particles to be moved in stream, index 0
                            should represent the desired entrainment location
                            i.e if particle[0] = 6 and the closest available 
                            vertex is 7, then the particle _desired_ to land
                            at 6 but is moved to 7.
        model_particles -- model array of model particles
        bed_particles -- model array of bed particles
        available_particles -- array of available vertices in the stream
    
    Returns:

    """
    entrainment_dict = {}
    for particle in event_particles: 
        orig_x = model_particles[model_particles[:,3] == particle[3]][0][0]
        verified_hop = find_closest_vertex(particle[0], available_vertices)
        
        if verified_hop == -1:
            exceed_msg = (
                f'Particle {int(particle[3])} exceeded stream...'
                f'sending to -1 axis'
            )
            print(exceed_msg) 
            particle[0] = verified_hop
        else:
            hop_msg = (
                f'Particle {int(particle[3])} entrained from {orig_x} '
                f'to {verified_hop}. Desired placement was: {particle[0]}'
            )
            print(hop_msg)
            particle[0] = verified_hop
            placed_x, placed_y = place_particle(particle, parameters.set_diam, 
                                          model_particles, bed_particles)
            particle[0] = placed_x
            particle[2] = placed_y
            
        entrainment_dict[particle[3]] = verified_hop
        model_particles[model_particles[:,3] == particle[3]] = particle
        
    available_vertices = np.setdiff1d(available_vertices, list(entrainment_dict.values()))
    
    return entrainment_dict, model_particles, available_vertices
    

def find_closest_vertex(desired_hop, available_vertices):
    """ Find the closest (greater than or equal) available
    vertex to the desired entrainment. 
    
    Keyword arguments:
    desired_hop -- float indicating the desired hop location
    available_location -- np array of available vertices in model
    
    Returns:
    vertex -- the closest available vertex that is >= desired_hop
    """    
    available_vertices = np.sort(available_vertices)
    forward_vertices = available_vertices[available_vertices >= desired_hop]
    
    if forward_vertices.size < 1:
        vertex = -1
    else:
        vertex = forward_vertices[0]
    return vertex    


# TODO: confirm the naming of function
def check_unique_entrainments(entrainment_dict):
    """ Check that all entrainments in the dictionary are unqiue. 
    
    This function will flag any input with non-unique
    entrainemnts. For a model with n entrainments, and k 
    particles entraining at the same vertex, a list of k-1 particles 
    will returned. The list represents those particles that 
    should be re-entrained (forced to a different vertex).
    
    Keyword arguments:
        entrainment_dict -- dictionary with key=id and value=vertex 
                            particle (id) is being entrained at
                            
    Returns:
        unique_flag -- boolean indicating if entrainment_dict had 
                        only unique entrainments (true) or at 
                        least one non-unique entrainemnt event (false)
        redo_list -- list of particles to be re-entrained in order to
                        achieve uniqueness. If unique_flag is True 
                        then list will be returned empty
    """
    redo_list = []
    unique_flag = True
    # create defaultdict struct to avoid missing key errors while grouping
    entrainment_groups = defaultdict(list)
    for p_id, vertex in entrainment_dict.items():
        entrainment_groups[vertex].append(p_id)
    
    entrainment_groups = dict(entrainment_groups)
    for vertex, p_id in entrainment_groups.items():
        if vertex == -1:
            pass
        elif len(p_id) > 1:
            unique_flag = False
            nonunique_msg = (
                f'More than one particle attempting to entrain at vertex: '
                f'{vertex}... Randomly selecting one particle to remain, all ' 
                f'others will be forced to the next available vertex.'
            )
            print(nonunique_msg)
            stay_particle = random.sample(p_id, 1)
            print(stay_particle)
            for particle in p_id:
                if particle != stay_particle[0]:
                    redo_list.append(int(particle))
        
    return unique_flag, redo_list     


# Taken from original model
def plot_stream(iteration, bed_particles, model_particles, x_lim, y_lim,
                available_vertices, to_file):
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
    plt.title(f'Iteration {iteration}')
    if to_file:
        filename = f'iter_{iteration}.png'
        path = 'plots/' + filename
        plt.savefig(path)
    else:
        plt.show()
    return
