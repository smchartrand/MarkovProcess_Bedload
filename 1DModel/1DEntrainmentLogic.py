import model_logic as ml 
    
bed_particles, x_max = ml.build_streambed()   
model_particles = ml.set_model_particles(bed_particles)


model = ml.Model(bed_particles, model_particles, x_max)   
model.run()
    
    