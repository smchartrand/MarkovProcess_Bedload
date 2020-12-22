import model_logic as ml 
import parameters as pm
import numpy as np
import os 

ITERATION_HEADER = ("""\n
    --------------------- Iteration {iteration} ---------------------
    """)
   
ITERATION_TEMPLATE = ("""\n
    # of entrainment events: {e_events}\n
    Particles to be entrained: {particles}\n                          
                      """)
                      
if not os.path.exists('plots'):
    os.makedirs('plots')
   
bed_particles, bed_length = ml.build_streambed()   
model_particles = ml.set_model_particles(bed_particles)

subregions = ml.define_subregions(bed_length, pm.num_subregions)
   
for iteration in range(pm.n_iterations):
    print(ITERATION_HEADER.format(iteration=iteration))  
    # Calculate the number of entrainment events per-unit-area
    e_events = np.random.poisson(pm.lambda_1, None)
    total_e_events, event_particles = ml.get_event_particles(
                                                            e_events, 
                                                            subregions, 
                                                            model_particles)
    print(ITERATION_TEMPLATE.format(
                                e_events=total_e_events, 
                                particles=event_particles))     
    
    model_particles = ml.move_model_particles(total_e_events, 
                                              event_particles, 
                                              model_particles, 
                                              bed_particles)
    
    ### FOR TESTING re-calculate avail vertice for plotting
    available_vertices = ml.compute_available_vertices(model_particles, 
                                                       bed_particles)
    ml.plot_stream(iteration, bed_particles, model_particles, 500, 30, 
                   available_vertices, to_file=True)
        
    