import model_logic as logic
import math
import parameters
import numpy as np
import random

# FOR TESTING: 
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")
 
###############################################################################
    
max_particles = int(math.ceil( parameters.x_max / parameters.min_diam ))
bed_particles = np.zeros([max_particles, 4],dtype=float)

bed_particles, bed_vertices = logic.build_streambed(bed_particles, parameters.set_diam)   


###############################################################################  

model_particles, chosen_vertex, valid_vertices, nulled_vertices = logic.set_model_particles(bed_vertices, bed_particles)
logic.plot_stream(bed_particles, model_particles, 150, 100/4, valid_vertices)

###############################################################################
e_events_store = np.zeros(parameters.n_iterations)

    
''' 
The following for loop will run the model n (n_iterations) times through.

In each loop, the number of entrainment events is calculated (e_events). Then
a random sample of particles is taken from the avaliable model particle set, whose
size is equal to the number of calculated entrainemnt events. If there are no 
particles in 'queue' at the start of the stream, then the random particles are
given immediately to move_model_particles. If there are particles in 'queue' 
then these particles are added to the raondomly selected set, and all are
passed to move_model_particles.
'''
for step in range(parameters.n_iterations):
    # calculate the number of entrainment events per unit area of bed
    e_events = np.random.poisson(parameters.lambda_1, None)

    if e_events == 0:
        e_events = 1 #???
        
    # total number of entrained particles -- used later when bed regions implemented
    # e_grains = e_events * bed_sampreg
    # randomly select model particles to entrain per unit area of bed
    print(model_particles)
    event_particles = random.sample(list(model_particles), e_events)
    event_particles = np.array(event_particles)
    
    ''' for this loop, identify the particles at x=-1 and add them all to the
    event_particles array '''
    # find all particles at x=-1 (particles in queue)
    ii = np.where(model_particles == -1)[0]
    
    for index in ii:
        model_particles[index][0] = 0 # send particle to 0 (starting point)
        event_particles = np.vstack((event_particles, model_particles[index]))
        e_events = e_events + 1
    
    valid_vertices, nulled_vertices = logic.move_model_particles(e_events, event_particles, valid_vertices, model_particles, bed_particles, bed_vertices, nulled_vertices)
    e_events_store[step] = e_events
    
    
    


