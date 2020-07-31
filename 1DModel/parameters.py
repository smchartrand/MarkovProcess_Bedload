#%% initial parameter for particles and plotting
# Area = np.zeros([1, 1], dtype=int, order='F')

# # Total area of placed circles. Initialize to 0.
# AreaTotal = 0

# # Center coord. of placed circles: r1:x,r2:y.
# CenterCoord = np.zeros(1, dtype=int, order='F')

# Diameter = np.zeros(1, dtype=int, order='F')

# CenterElev = np.zeros(1, dtype=int, order='F')

Color_iter = 1000

Step_1 = 1

# Define origin for plotting.
XCoordinates_Orig = 0

# Packing density of the bed; ranges from > 0 to <~0.70.
Pack = 0.52

# Length of the domain in the streamwise direction in millimeters.
x_max = 150

# Spacing of nodes in millimeters.
Grid_Spacing = 1

# Minimum grain diameter in millimeters.
# IF VARIABLE USE THESE:
min_diam = 4.0
max_diam = 6.0 
# IF SET DIAM USE THIS:
set_diam = 5.0
set_radius = (set_diam / 2.0)

# Size of largest grain relative to the cross-stream length.
Size_Frac = 0.06



# number of model iterations  
n_iterations = 50
lambda_1 = 1
# Particle travel time minimum (t).
T_pmin = 0
# Particle travel time maximum (t).
T_pmax = 1.0

