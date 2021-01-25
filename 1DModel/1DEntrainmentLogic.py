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

particle_flux_dict = {}
particle_flux_list = []
for iteration in range(pm.n_iterations):
    print(ITERATION_HEADER.format(iteration=iteration))  
    # Calculate the number of entrainment events per-unit-area
    e_events = np.random.poisson(pm.lambda_1, None)
    event_particles = ml.get_event_particles(e_events, subregions, model_particles)
    print(ITERATION_TEMPLATE.format(
                                e_events=len(event_particles), 
                                particles=event_particles))     
    
    model_particles, particle_flux = ml.run_entrainments(model_particles, bed_particles, event_particles, pm.normal_dist)
    particle_flux_dict[iteration] = particle_flux
    particle_flux_list.append(particle_flux)
    
    # Re-calculate avail vertice for plotting. Can be 
    available_vertices = ml.compute_available_vertices(model_particles, 
                                                       bed_particles)
    ml.plot_stream(iteration, bed_particles, model_particles, pm.x_max, 10, 
                   available_vertices, to_file=False)
    
    
        