"""
Parameters file. 

Here you can set the parameter values
of the run. Above each paramter is a
comment describing the parameter, 
it's required type and
the values that have been tested 
for it. 

Most parameters have suggested ranges
for input to stay within. If thes user
selects beyond the strict ranges or
sets to an incorrect type
the model will not run.

Note: this file will hopefully be 
converted to a YML in future 
iterations of the model.
"""

#%% Mutable parameters.

# The packing density of the model particles
# RANGE: > 0 to <~0.70.
Pack = 0.52

# Length of the domain in the streamwise direction (mm).
# RANGE: > 0
# TYPE: INT
# TESTED ON: 50-300 mm
x_max = 100
    
# Grain diameter (mm) 
# RANGE: > 0
# TESTED: 0.5, 1.0, 1.5, 2.0
set_diam = 1

# The number of bed subregions
# RANGE: > 0
# TYPE: INT
# TESTED: 1, 2, 3, 4, 5
num_subregions = 2

# The maximum number of levels permitted in stream
# TESTED: 2, 3, 4, 5,
level_limit = 4

# The number of iterations to run
# RANGE: > 0
# TYPE: INT
# TESTED: 5-400
n_iterations = 100

# 
lambda_1 = 1

# Minimum particle travel time (s)
# RANGE: > 0, < T_pmax
# TYPE: Float
# TESTED: 0.0
T_pmin = 0

# Maximum particle travel time (s).
# RANGE: > 0, > T_pmin
# TYPE: Float
# TESTED: 1.0
T_pmax = 1.0

# Distribution used during entrainment hop calculation
# If False, distribution = logNormal
# If True, distribution = Normal
# TYPE: Boolean
normal_dist=False

#%% Parameters that should not be touched.
set_radius = (set_diam / 2.0)
