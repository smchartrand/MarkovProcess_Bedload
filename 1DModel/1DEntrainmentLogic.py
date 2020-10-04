import model_logic as logic
import parameters
import numpy as np
import random



 
###############################################################################
    
bed_particles = np.zeros([parameters.max_particles, 5],dtype=float)

bed_particles, bed_vertices = logic.build_streambed(bed_particles, parameters.set_diam)   

###############################################################################  

model_particles, chosen_vertex, available_vertices, nulled_vertices = logic.set_model_particles(bed_vertices, bed_particles)
logic.plot_stream(bed_particles, model_particles, 250, 100/4, available_vertices)
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
    # print(model_particles)
    in_stream_particles = model_particles[model_particles[:,0] != -1]
    active_particles =  in_stream_particles[in_stream_particles[:,4] != 0]
    
    # randomly sample particles from model_particles to entrain
    idx_of_entrainment = random.sample(range(len(active_particles)), e_events)
    
    # identify particles at x=-1 and add them all to the event particles array
    # note: particles at x=-1 are those particles in 'queue'
    ii = np.where(model_particles[:,0] == -1)[0]
    for index in ii:
        model_particles[index][0] = 0 # send particle to 0 (starting point)
        idx_of_entrainment = np.append(index)
        e_events = e_events + 1
    
    num_event_particles = len(idx_of_entrainment)
    
    # 
    if e_events != num_event_particles:
        # 1. log that the # entrainment events has been altered
        
        # 2. update e_events
        e_events = num_event_particles
    
    
    available_vertices = logic.compute_available_vertices(model_particles, bed_particles)
    
    print(f'particles selected for entrainment: {idx_of_entrainment}')
    available_vertices = logic.move_model_particles(e_events, idx_of_entrainment, model_particles, bed_particles, available_vertices)
    e_events_store[step] = e_events

    
    
    


