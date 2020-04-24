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
x_max = 500

# Spacing of nodes in millimeters.
Grid_Spacing = 1

# Minimum grain diameter in millimeters.
min_diam = 4.0

set_diam = 5.0

# Size of largest grain relative to the cross-stream length.
Size_Frac = 0.06

max_diam = 6.0 