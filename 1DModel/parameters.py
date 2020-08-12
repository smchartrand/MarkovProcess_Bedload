import math

#%% initial parameter for particles and plotting

''' Particle and Streambed Parameters '''
# Packing density of the bed; ranges from > 0 to <~0.70.
Pack = 0.52

# Length of the domain in the streamwise direction in millimeters.
x_max = 250

# Spacing of nodes in millimeters.
Grid_Spacing = 1

# Minimum grain diameter in millimeters.
# IF GRAIN DIAM VARIABLE USE THESE:
#   min_diam = 4.0
#   max_diam = 6.0 
# IF CONSTANT DIAM USE THIS:
set_diam = 5.0
set_radius = (set_diam / 2.0)

# number of bed particles based on stream length
max_particles = int(math.ceil( x_max / set_diam ))


''' Model Run Parameters '''
# number of model iterations  
n_iterations = 50
lambda_1 = 1

# Particle travel time minimum (t).
T_pmin = 0
# Particle travel time maximum (t).
T_pmax = 1.0


