import model_logic as ml
import parameters
import numpy as np



 
    
bed_particles, bed_vertices = ml.build_streambed()   
model_particles, available_vertices = ml.set_model_particles(bed_vertices, 
                                                             bed_particles)
e_events_store = np.zeros(parameters.n_iterations)
ml.plot_stream(bed_particles, model_particles, 250, 100/4, available_vertices)


model = ml.Model(bed_particles, model_particles, e_events_store)   
model.run()
    
    


