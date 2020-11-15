import model_logic as ml 
import parameters as param
    
bed_particles, x_max = ml.build_streambed()   
model_particles = ml.set_model_particles(bed_particles)


model = ml.Model(bed_particles, model_particles, x_max, param.num_subregions,
                 param.lambda_1, param.n_iterations)   
model.run()

### Data processing below
    
    