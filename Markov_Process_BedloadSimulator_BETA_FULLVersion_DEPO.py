# -*- coding: utf-8 -*-
# """
# Created on Thu May 24 12:25:15 2018
# @author: schartrand
# """
import math
import random
import time
import numpy as np
# import matplotlib.cm as cm
#
from numba import jit
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import axes3d
from scipy.spatial import Voronoi, voronoi_plot_2d
# Timer
Start_time = time.time()
#
# BEGIN CODE
# This code was developed by Shawn Chartrand and David Furbish, March 2018.
# The code is the basis of Markov birth-death simulations of bedload transport
# in rivers under relatively low transport rates. The code is split into two
# parts. Part 1 builds the random streambed surface. Part 2 entrains, deposits
# and moves particles along the streambed surface.
# The code randomly populates a rectangular area of specified dimension with
# circles of random and nonuniform size. Placement of circles proceeds using
# a grid mask to identify open coordinates. Placement of circles stops when a
# specified packing fraction (circle volume\total volume) is reached. After
# bed surface creation, grains are stochastically entrained and deposited.
###############################################################################
# PART 1 ######################################################################
###############################################################################
# STEP ONE: DEFINE DATA STRUCTURES, INITIALIZE SOME VARIABLES.
# DATA LISTS
# Circle (i.e. a particle) area.
Area = np.zeros([1, 1], dtype=int, order='F')
# Total area of placed circles. Initialize to 0.
AreaTotal = 0
# Center coord. of placed circles: r1:x,r2:y.
CenterCoord = np.zeros([2, 1], dtype=int, order='F')
# Final placed circle diameters.
Diameter = np.zeros([1, 1], dtype=int, order='F')
# Final placed circle center elevation (considered as a sphere).
CenterElev = np.zeros([1, 1], dtype=int, order='F')
# COUNTER AND VARAIBLES
# Defines size of colormap.
Color_iter = 1000
# Milestones used to display calculation progress.
Milestones = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 100]
# Initialize the step counter for part 1.
Step_1 = 1
# Define origin for plotting.
XCoordinates_Orig = 0
YCoordinates_Orig = 0
# END SECTION
###############################################################################
# STEP TWO: BUILD THE VIRTUAL DOMAIN WITH ALL NECESSARY DETAILS
# KEY INPUT PARAMETER VALUES
# Packing density of the bed; ranges from > 0 to <~0.70
Pack = 0.52
# Length of the domain in the streamwise direction in millimeters
x_max = 50
# Length of the domain in the cross-stream direction in millimeters
y_max = 50
# Spacing of nodes
Grid_Spacing = 0.5
# Minimum grain diameter
MinDiam = 3.0
# Size of largest grain relative to the cross-stream length
Size_Frac = 0.06
# Time threshold to restart script (seconds)
Time_threshold = 10
# DEPENDENT CALCULATIONS
# Total virtual streambed area.
Total_VArea = x_max * y_max
# Max grain size fraction of height.
CircDiam_Fac = 1 / Size_Frac
# Specify the max circle diameter.
Max_Diam = 6.0 #y_max / CircDiam_Fac
# END SECTION
###############################################################################
# STEP THREE: BUILD MATRIX OF COORDINATE ADDRESSES FOR CENTER OF CIRCLES
# Build x array for grid.
X_array = np.arange(0, (x_max) + Grid_Spacing, Grid_Spacing, dtype=float)
# Build y array for grid.
Y_array = np.arange(0, (y_max) + Grid_Spacing, Grid_Spacing, dtype=float)
# Length of array dim 1 and dim 2
Length1 = int(len(Y_array))
Length2 = int(len(X_array))
# For logical mask.
XYCoordinateslogical = np.zeros((Length1, Length2), dtype=int)
# Eliminate periphery cells to limit grain size to minimum specified diameter
MatIdx = np.int(MinDiam / Grid_Spacing)
XYCoordinateslogical[0:MatIdx, :] = 1
XYCoordinateslogical[:, 0:MatIdx] = 1
XYCoordinateslogical[Length1 - MatIdx:Length1, :] = 1
XYCoordinateslogical[:, Length2 - MatIdx:Length2] = 1
# Length of matrix dim 1 and dim 2
Dim1 = len(XYCoordinateslogical)
Dim2 = len(XYCoordinateslogical[0])
# This defines the logical array length for the mask.
SIZE = Dim1 * Dim2
# This variable is used for the logical mask.
idx = np.linspace(0, SIZE-1, SIZE, dtype=int)
# BUILD A COORDINATE MATRICES.
YCoordinates = np.repeat(Y_array[None], Length2, axis=0).transpose()
XCoordinates = np.repeat(X_array[None], Length1, axis=0)
# END SECTION
###############################################################################
# STEP FOUR: DEFINE TWO FUNCTIONS USED TO PLACE and REVISE CIRCLE DIAMETERS


def fu_radcheck(x_Center, y_Center, radius):
    # These operations are used to check if the randomly placed
    # circle falls outside the rectangular area
    x_circ_maxck = x_max - (x_Center + radius)
    x_circ_minck = (x_Center - radius) - XCoordinates_Orig
    y_circ_maxck = y_max - (y_Center + radius)
    y_circ_minck = (y_Center - radius) - YCoordinates_Orig
    circ_array = np.arange(0, 4)
    # Build array to check the cumulative status of overlap
    circ_check = (np.reshape(np.asarray([x_circ_maxck, x_circ_minck,
                                         y_circ_maxck, y_circ_minck]), (-1)))
    circ_check_logical = (circ_check <= 0)
    # Check to see if the circle radius needs to be changed.
    if np.sum(circ_check_logical) > 0:
        # Use logical indexing to speed things up.
        circ_logidx = np.where(circ_check_logical == 1)
        circ_check_idx = circ_array[circ_logidx]
        circ_check_radrev = min(circ_check[circ_check_idx])
        # Revise the circle radius
        radrev = radius + circ_check_radrev
        radius = radrev
    if Step_1 == 1:
        newCircle_Found = 1
    else:
        newCircle_Found = 0
    return radius, newCircle_Found


def fu_diamrevise(x_Center, y_Center, radius):
    # The following operations determines if the test circle overlaps
    # with existing circles, but only to the resolution of the
    # underlying grid. The outcome is a logical matrix. Values of 1
    # mean that grid point overalps between the test circle and an
    # existing circle.
    Ctr_dist_build = (((CenterCoord[0, :] - x_Center) ** 2 +
                       (CenterCoord[1, :] - y_Center) ** 2) ** 0.5)
    radius_check = (Diameter / 2) + radius
    radius_overlap = (radius_check >= Ctr_dist_build)
    # Now actually check the overalp condition.
    if np.sum([radius_overlap]) == 0:
        if radius >= (MinDiam / 2):
            # The new circle does not overlap so proceed.
            newCircle_Found = 1
        else:
            newCircle_Found = 0
    elif np.sum([radius_overlap]) == 1:
        # The new circle overlaps with one other circle
        overlap = (np.arange(0, len(radius_overlap[0]), dtype=int).
                   reshape(1, len(radius_overlap[0])))
        overlap_logix = (radius_overlap == 1)
        idx_true = overlap[overlap_logix]
        radius_temp = Ctr_dist_build[idx_true] - (Diameter[0, idx_true] / 2)
        # New circle radius can be negative when inside another circle.
        if radius_temp >= (MinDiam / 2):
            newCircle_Found = 1
            radius = radius_temp
        else:
            newCircle_Found = 0
    elif np.sum([radius_overlap]) >= 2:
        # Overlap with two or more circles.
        # Set newCircle flag to zero so the step is repeated.
        newCircle_Found = 0
    return radius, newCircle_Found
# END SECTION
###############################################################################
# STEP FIVE: RANDOMLY PLACE CIRCLES IN RECTANGULAR AREA


while (AreaTotal / Total_VArea) < Pack:
    # Initial calculations for all steps.
    # Flatten the logical matrix to an array
    coord_mat_transform = XYCoordinateslogical.flatten('F')
    # Random circle diameter
    diameter_Test = random.randint(MinDiam, Max_Diam)
    # Index used to track open coordinate
    test = (coord_mat_transform == 0)
    # Index used to track open x,y coordinate pairs
    idx2 = idx[test]
    # Random number used to select and open x,y coordinate pair
    randCoord = random.randint(1, (len(idx2)-1))
    # Choose open coord. with random number.
    idx3 = idx2[randCoord]
    # Specify random x and y center coordinate.
    x_Center = (XCoordinates[np.unravel_index([idx3], [Length1, Length2],
                                              order='F')])
    y_Center = (YCoordinates[np.unravel_index([idx3], [Length1, Length2],
                                              order='F')])
    radius = diameter_Test * 0.5
    if radius < (MinDiam * 0.5):
        print("--- YIKES ---")
        break
    # Check to see if circle plots outside of rectangular domain
    ###############
    # CALL FUNCTION
    ###############
    radius, newCircle_Found = fu_radcheck(x_Center, y_Center, radius)
    # Check to see if radius needs revision for steps greater than the first
    if Step_1 != 1:
        ###############
        # CALL FUNCTION
        ###############
        radius, newCircle_Found = fu_diamrevise(x_Center, y_Center, radius)
    # Only write data if a circle was placed.
    if newCircle_Found == 1:
        # This operation eliminates available grid points for placement.
        mask = ((((XCoordinates - x_Center) ** 2 +
                  (YCoordinates - y_Center)** 2) ** 0.5 <= radius))
        # Track which grid points are available for circles and
        # which are not and pass to the next step in the loop.
        XYCoordinateslogical[np.where(mask == 1)] = 1
        # Write data to variables.
        Circle_Area = (np.reshape(np.asfortranarray(math.pi * (radius ** 2),
                                                    dtype=float), (1, +1)))
        Area = np.hstack((Area, Circle_Area))
        AreaTotal = np.sum(Area)
        CenterCoord = np.hstack((CenterCoord, [x_Center, y_Center]))
        Diameter = np.hstack((Diameter, np.reshape(radius * 2, (1, +1))))
        CenterElev = np.hstack((CenterElev, np.reshape(radius, (1, +1))))
        if Step_1 == 1:
            Area = np.delete(Area, 0, 1)
            CenterCoord = np.delete(CenterCoord, 0, 1)
            Diameter = np.delete(Diameter, 0, 1)
            CenterElev = np.delete(CenterElev, 0, 1)
        # Calculate the percent completion of simulation and print to console.
        percentage_complete = (100.0 * (AreaTotal / Total_VArea) / Pack)
        while len(Milestones) > 0 and percentage_complete >= Milestones[0]:
            print "{}% complete".format(Milestones[0])
            # remove that milestone from the list
            Milestones = Milestones[1:]
        if (AreaTotal / Total_VArea) < Pack:
            # Advance counter.
            Step_1 = Step_1 + 1
    # This operation restarts the script if the execution time is too long.
    # This occurs if the underlying grid is too coarse.
    Current_time = time.time()
    if (Current_time - Start_time) > Time_threshold:
        Start_time = time.time()
        Step_1 = 1
        percentage_complete = 0
        # Circle area.
        Area = np.zeros([1, 1], dtype=int, order='F')
        # Total area of placed circles. Initialize to 0.
        AreaTotal = 0
        # Center coord. of placed circles: r1:x,r2:y.
        CenterCoord = np.zeros([2, 1], dtype=int, order='F')
        # Final placed circle diameters.
        Diameter = np.zeros([1, 1], dtype=int, order='F')
        # Final placed circle center elevation (considered as a sphere).
        CenterElev = np.zeros([1, 1], dtype=int, order='F')
        # For logical mask.
        XYCoordinateslogical = np.zeros((Length1, Length2), dtype=int)
        # Milestones used to display calculation progress.
        milestones = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 100]
        print("--- RESTARTING LOOP ---")
# END SECTION
###############################################################################
# STEP SIX: PRINT RESULTS TO SCREEN AND WRITE PERTINENT RESULTS TO FILE
if np.any(Diameter == 0):
    print("--- WARNING - Some particles have a zero diameter ---")
else:
    print("--- GOOD NEWS - All particles have a diameter ---")
print("--- %s seconds ---" % (time.time() - Start_time))
print("--- %s total particles in domain ---" % (Step_1))
# STEP SIX: PLOT THE RESULTS
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlim((-2, x_max + 2))
#ax.set_xticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
ax.set_ylim((-2, y_max + 2))
ax.minorticks_on()
ax.tick_params(axis='x', which='minor', direction='out')
ax.tick_params(axis='y', which='minor', direction='out')
plt.xlabel('downstream')
plt.ylabel('cross-stream')
plt.title('N = %i particles' %Step_1, fontsize=12, style='italic')
resolution = 50  # the number of vertices
Radius_Array = np.asarray((Diameter * 0.5), dtype=float)
Radius_Array = Radius_Array.reshape(-1)
XCenter = np.reshape(CenterCoord[0, :], (1, len(CenterCoord[0])))
XCenter = XCenter.reshape(-1)
YCenter = np.reshape(CenterCoord[1, :], (1, len(CenterCoord[0])))
YCenter = YCenter.reshape(-1)
plt.rcParams['image.cmap'] = 'gray'
# This method of plotting circles comes from Stack Overflow questions\32444037
# Note that the patches won't be added to the axes, instead a collection will.
patches = []
for x1, y1, r in zip(XCenter, YCenter, Radius_Array):
    circle = Circle((x1, y1), r)
    patches.append(circle)
p = (PatchCollection(patches, color="#BDBDBD", alpha=0.9, linewidths=(0, )))
# colors = Color_iter*np.random.rand(len(CenterCoord[0]))
# p = (PatchCollection(patches, cmap=cm.coolwarm,
#                     alpha=1.0, linewidths=(0, )))
# p.set_array(colors)
# p.set_clim([5, 950])
ax.add_collection(p)
plt.show()
# Need to insert some commands to write files to subdirectory.
# os.makedirs('\\ScriptTest')
# fig.savefig('.\ScriptTest\Starting_Bed.pdf', format='pdf', dpi=2400)
fig.savefig('./ScriptTest/TotalBedMotions00000.png', format='png', dpi=600)
# Save initial results for archiving and plotting.
np.save('.\ScriptTest\XCenter_Initial', XCenter)
np.save('.\ScriptTest\YCenter_Initial', YCenter)
np.save('.\ScriptTest\Radius_Array_Initial', Radius_Array)
# END OF PART 1 VIRTUAL BED CREATION
###############################################################################
# PART 2 ######################################################################
###############################################################################
# STEP ONE: DEFINE DATA STRUCTURES, INITIALIZE SOME VARIABLES.
# DATA LISTS
# Parameter to store coordinate and diameter data from step 1
Grains = np.vstack((CenterCoord, Diameter, CenterElev))
np.save('.\ScriptTest\ParticleDetails_Initial', Grains)
# Upstream particle supply rate (particles\t).
Nu_in = np.zeros([1, 1], dtype=int, order='F')
# Particle emigration rate (particles\t) at boundaries.
Nu_out = [0]
# Particle emigration rate (particles\t) at downstream boundary.
Nu_out_ds = [0]
# Particle birth or entrainment rate within the control volume (particle\t).
# Set birth rate as a constant value for trial runs
Lambda_1 = 3
# Entrainment events per unit area of the bed and unit time (events\t)
E_events = np.zeros([1, 1], dtype=int, order='F')
# Particle travel time minimum (t).
T_pmin = 0
# Particle travel time maximum (t).
T_pmax = 1.0
# Particle hop distance (L).
L_x = np.zeros([1, 1], dtype=int, order='F')
E_events_Store = [0]
E_grains_Store = [0]
L_x_Store = [0]
Hop_locx_Store = [0]
# Initialize the step counter for part 2.
Step_2 = 1
plot_step = 1
# Temporary loop completion criteria
Loops = 2
# END SECTION
###############################################################################
# STEP TWO: DIVIDE THE BED INTO SAMPLING REGIONS AND BUILD SAMPLING ARRAY
# Number of bed sampling regions within the control volume V_c
Bed_sampreg = 2
# Data storage arrays
Nu_out_Store = np.zeros([Bed_sampreg, Loops], dtype=int, order='F')#[0]
# Bed sampling boundaries in the x-direction
BS_boundaries = (XCoordinates_Orig + (x_max / Bed_sampreg) *
                 np.arange(0, Bed_sampreg + 1, dtype=int))
SSamp_len = len(BS_boundaries)
# Index to hold boundaries between subsampling regions to facilitate random
# sampling from within each region.
SubSampl_idx = np.zeros([1, SSamp_len], dtype=int, order='F')
xB_idx = np.zeros([1, 2], dtype=int, order='F')
yB_idx = np.zeros([1, 2], dtype=int, order='F')
# END SECTION
###############################################################################
# STEP TWO: DEFINE A SEARCH FUNCTION AND ONE TO SPECIFY PARTICLE KINEMATICS


# This search code using numba is from Stack Overflow
# https://stackoverflow.com/questions/7632963/
# I have modified the search criteria to suit the present purpose
@jit(nopython=True)
def fu_search(item, vec):
    """return the index of the first occurence of item in vec"""
    length = len(vec)
    for nn in xrange(length):
        if vec[nn] >= item:
            return nn
        if np.max(vec) < item:
            return length - 1
    return -1
    # Clear the loop variable from memory
    del nn


def fu_transport(E_events, BS_boundaries, Grains_sort):
    """"Randomly move particles across the virtual streambed"""
    for n in range(0, SSamp_len):
        if n == 0:
            # Initialize variables to store entrainment and hop information
            Entrain_idx = np.zeros([1, 1], dtype=int, order='F')
            Entrain_locsx = np.zeros([1, 1], dtype=int, order='F')
            Entrain_locsy = np.zeros([1, 1], dtype=int, order='F')
            Rand_ELocs_Store = np.zeros([1, 1], dtype=int, order='F')
            # Particle travel time (t).
            T_p = np.zeros([1, 1], dtype=int, order='F')
            T_p_2 = np.zeros([1, 1], dtype=int, order='F')
            # Particle travel distance (mm).
            L_x = np.zeros([1, 1], dtype=int, order='F')
            Entrain_diam = np.zeros([1, 1], dtype=int, order='F')
            # Keep track of the particle details array length in dimension 0
            Length_idx = len(Grains_sort[0])
        ###############
        # CALL FUNCTION
        ###############
        # Find indices within sort coordinates array for subsampling
        SubSampl_idx[0, n] = fu_search(BS_boundaries[n], Grains_sort[0, :])
        # Specify random entrainment sampling indices
        if n == 0:
            continue
        elif n == Bed_sampreg:
            Rand_ELocs = (random.sample(xrange(SubSampl_idx[0, n-1] + 1,
                                               Length_idx), E_events))
            Rand_ELocs = np.array(Rand_ELocs).reshape(1, E_events)
        else:
            Rand_ELocs = (random.sample(xrange(SubSampl_idx[0, n-1] + 1,
                                               SubSampl_idx[0, n]), E_events))
            Rand_ELocs = np.array(Rand_ELocs).reshape(1, E_events)
        # Save the random locations for replacement operation below
        Rand_ELocs_Store = np.hstack((Rand_ELocs_Store, Rand_ELocs))
        # Actual entrainment x-coordinates
        Loc_x = Grains_sort[0, Rand_ELocs]
        # Actual entrainment y-coordinates
        Loc_y = Grains_sort[1, Rand_ELocs]
        # Associated entrained particle diameters
        E_diam = Grains_sort[2, Rand_ELocs]
        # Calculate the travel time as a randomly sampled variable from a
        # uniform distribution constrained by Fathel et al., 2015.
        if n == 0 and Nu_in > 0:
            continue
        else:
            T_p_init1 = (np.random.uniform(T_pmin, T_pmax, E_events).
                         reshape(1, E_events))
        # https:\\stackoverflow.com\questions\2106503\
        # Now distribute the random variables over a pseudo exponential
        # distribution based on inverse transform sampling.
        T_p_init2 = np.log(1 - T_p_init1) / (-Lambda_1)
        T_p_alt2 = -10 * np.log(1-T_p_init1)
        # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
        # For now I am multuplying by 10 so units are consistent (). Rounding
        # the output to map the entrained particles to the 2D grid.
        L_x_init = np.round((T_p_init2 ** 2) * 10 * 2, 1)
        # Write the entrainment location data to three variables
        Entrain_idx = np.hstack((Entrain_idx, Rand_ELocs))
        Entrain_locsx = np.hstack((Entrain_locsx, Loc_x))
        Entrain_locsy = np.hstack((Entrain_locsy, Loc_y))
        T_p = np.hstack((T_p, T_p_init2))
        T_p_2 = np.hstack((T_p_2, T_p_alt2))
        L_x = np.hstack((L_x, L_x_init))
        Entrain_diam = np.hstack((Entrain_diam, E_diam))
        # Delete the initialization zero entry
        if n == (Bed_sampreg-1):
            Entrain_idx = np.delete(Entrain_idx, 0, 1)
            Entrain_locsx = np.delete(Entrain_locsx, 0, 1)
            Entrain_locsy = np.delete(Entrain_locsy, 0, 1)
            Rand_ELocs_Store = np.delete(Rand_ELocs_Store, 0, 1)
            T_p = np.delete(T_p, 0, 1)
            L_x = np.delete(L_x, 0, 1)
            Entrain_diam = np.delete(Entrain_diam, 0, 1)
    return n, Entrain_locsx, Entrain_locsy, L_x, T_p, T_p_2, Rand_ELocs_Store, Entrain_diam, Entrain_idx

###############################################################################
# STEP TWO: ENTRAIN, MOVE AND DEPOSIT PARTICLES IN CONTROL VOLUME
# Particles are entrained, moved and deposited for a specified number of Loops.


while Step_2 <= Loops:
    # E_events is the entrainment events per unit area of the bed
    E_events = np.random.poisson(Lambda_1, None)
    if E_events == 0:
        continue
        E_events = 1
    # Save the E_events to a list
    E_events_Store.extend([E_events])
    # Total number of entrained particles not counting immigrants across V_c
    E_grains = E_events * Bed_sampreg
    # Save the E_events to a list
    E_grains_Store.extend([E_grains])
    if Step_2 == 1:
        # Make it easy to figure out where to randomly entrain particles.
        Grains_sort = Grains[:, Grains[0].argsort()]
    else:
        # Make it easy to figure out where to randomly entrain particles
        # within the CV.
        Grains_sort = Grains_sort[:, Grains_sort[0].argsort()]
    ###############
    # CALL FUNCTION
    ###############
    n, Entrain_locsx, Entrain_locsy, L_x, T_p, T_p_2, Rand_ELocs_Store, Entrain_diam, Entrain_idx = fu_transport(E_events, BS_boundaries, Grains_sort)
    # This is the destination of partciles after the hop distance is applied.
    # The result returns the destination in the streamwise coordinate only.
    Hop_locx = Entrain_locsx + L_x
    Hop_locy = Entrain_locsy
    Hop_diam = Entrain_diam
    L_x_Store.extend([L_x])
    Hop_locx_Store.extend([Hop_locx])
    if (Step_2 == Loops) & (n == SSamp_len - 1):
        del L_x_Store[0]
        del Hop_locx_Store[0]
        del E_events_Store[0]
        del E_grains_Store[0]
    # Compile the grains which pass the subregions and the downstream end
    # Begin the range at 1 to skip the boundary at x = 0.
    for nn in range(1, SSamp_len):
        Nu_out = (np.count_nonzero((Entrain_locsx <= BS_boundaries[nn]) &
                                   (Hop_locx > BS_boundaries[nn])))
        Nu_out_Store[(nn - 1), (Step_2 - 1)] = Nu_out
    # Clear the loop variables from memory
    del n, nn
    # Deal with particles that cross the downstream boundary
    Nu_in = np.count_nonzero(Hop_locx > x_max)
    Hop_idx = np.argwhere(Hop_locx > x_max)
    if Nu_in == 1:
        Hop_locx[Hop_idx[0,0], Hop_idx[0,1]] = (Hop_locx[Hop_idx[0,0],
                 Hop_idx[0,1]] - x_max)
    elif Nu_in > 1:
        for n in range(0, Nu_in) :
            Hop_locx[Hop_idx[n,0], Hop_idx[n,1]] = (Hop_locx[Hop_idx[n,0],
                     Hop_idx[n,1]] - x_max)
    Hop_len = len(Hop_locx[0])
    ###############################################################
    # Move particles to physically reasonable positions in 3D space
    ###############################################################
    for n in range(0, Hop_len):
        # Specify the search neighborhood bounds for particle placement.
        # X-coordinate first.
        if n == 0:
            LxB = 0
        else:
            LxB = Hop_locx[0, n] - Max_Diam
            if LxB < 0:
                LxB = 0
        HxB = Hop_locx[0, n] + Max_Diam
        if HxB > x_max:
            HxB = x_max
        Overlap_xbound = np.round((LxB, HxB), 1)
        # Y-coordinate second
        LyB = Hop_locy[0, n] - Max_Diam
        if LyB < 0:
            LyB = 0
        HyB = Hop_locy[0, n] + Max_Diam
        if HyB > y_max:
            HyB = y_max
        Overlap_ybound = np.round((LyB, HyB), 1)
        for nn in range(0, len(Overlap_xbound)):
            ###############
            # CALL FUNCTION
            ###############
            # Find indices within sort coordinates array for subsampling.
            # X-coordinate first.
            xB_idx[0, nn] = fu_search(Overlap_xbound[nn], Grains_sort[0, :])
            xstart = xB_idx[0, 0]
            xstop = xB_idx[0, 1]
            #xB_idx = np.arange(xstart, xstop, 1)
            # Y-coordinate second.
            ystart = Overlap_ybound[0]
            ystop = Overlap_ybound[1]
            yB_idx = (np.where((Grains_sort[1, :] >= ystart) &
                               (Grains_sort[1, :] <= ystop))[0])
            # Find where the x and y indices match to focus the search
            # neighborhood this operation returns the indices for yB_idx where
            # the values of yB_idx are in between the xB_idx start and stop
            # values (indices of the xB_idx array). Need to grab the index
            # values after the first operation.
            Sear_yidx = np.where((yB_idx >= xstart) & (yB_idx <= xstop))[0]
            #Rmv_idx = np.where
            # This is the key search index.
            Sear_nidx = yB_idx[Sear_yidx]
        # Delete the entry for the entrained grain
        Rem_idx = np.where(Sear_nidx == Entrain_idx[0, n])[0]
        Sear_nidx = np.delete(Sear_nidx, Rem_idx)
        # Need to narrow down the search vector for the Hop_loc arrays
        Ctr_dist_trans = (((Grains_sort[0, Sear_nidx] - Hop_locx[0, n]) ** 2 +
                     (Grains_sort[1, Sear_nidx] - Hop_locy[0, n]) ** 2) ** 0.5)
        # I think this should be within the if statement above
        # Build the second part of the particle overlp test check
        Nei_grains = (Grains_sort[2, Sear_nidx] / 2) + (Entrain_diam[0, n] / 2)
        Grains_overlap = (Nei_grains >= Ctr_dist_trans)
        # Now actually check the overalp condition.
        if np.sum([Grains_overlap]) == 0:
            continue
        else:
            # Sub-step 1: Build a grid within the search neighborhood
            # This is the grid spacing
            Nei_len = len(Grains_overlap)
            Nei_GS = np.array(0.10)
            if Nei_GS > Grid_Spacing:
                print("--- YIKES. Neighborhood grid is too coarse ---")
            # Build x array for grid. This may throw an error at some point.
            XNei_array = (np.round(np.arange(Grains_sort[0, xstart], Grains_sort[0, xstop] + 0.1,Nei_GS, dtype = float),2))
            # Build y array for grid.
            YNei_array = np.linspace(HyB, LyB, (HyB-LyB)/Nei_GS + 1.0 )
            # The next 3 lines are used for debugging. Remove after working.
            X, Y = np.meshgrid(XNei_array, YNei_array)
            XNei = X.reshape((np.prod(X.shape),))
            YNei = Y.reshape((np.prod(Y.shape),))
            # Length of array dim 1 and dim 2
            LengthNei1 = int(len(YNei_array))
            LengthNei2 = int(len(XNei_array))
            # For logical mask.
            XYNei_Coordlogical = np.zeros((LengthNei1, LengthNei2), dtype=int)
            XYNei_Coordlogical_alt = np.zeros((LengthNei1, LengthNei2), dtype=int)
            # Fill up the grid with particles based on the search above.
            # Somewhere I have a significant digit issue. Need to find it.
            Nei_XY = np.round(zip(Grains_sort[0, Sear_nidx],Grains_sort[1, Sear_nidx]), 1)
            HElev = np.zeros((LengthNei1, LengthNei2), dtype=float)
            Surf = np.zeros((LengthNei1, LengthNei2), dtype=float)
            # Straight array based index of X and Y dimensions
            X_Idx = np.arange(0, LengthNei2, 1)
            Y_Idx = np.arange(0, LengthNei1, 1)
            for nnn in range(0, Nei_len):
                # Sometimes this operation throws an error because the NEI_XY
                # value does not occur in the Nei_array - this is due to the
                # search algorithm. Also there seem to be rounding errors with
                # the hop_locx sometimes leading to > 1 decimal place.
                XNei_Idx = np.where(XNei_array==Nei_XY[nnn][0])
                YNei_Idx = np.where(YNei_array==Nei_XY[nnn][1])
                XYNei_Coordlogical[YNei_Idx, XNei_Idx] = 1
                # Fill in the mask with the circle coordinates <= radii
                mask = ((((X - Nei_XY[nnn][0]) ** 2 +
                  (Y - Nei_XY[nnn][1])** 2) ** 0.5 <= (Nei_grains[nnn] * 0.5)))
                # Track which grid points are available for circles and
                # which are not and pass to the next step in the loop.
                XYNei_Coordlogical[np.where(mask == 1)] = 1
                # Use the X, Y mesh above to build the surface but only grab
                # coordinates for one particle at a time
                #Elev_XY = zip(X[mask], Y[mask])
                #
                xx = np.empty_like (X)
                xx[:] = X
                xx[np.where(mask==0)] = 0
                xxx=np.abs(xx-X[YNei_Idx,XNei_Idx])
                xxx[xxx>(2*Max_Diam)]=0
                yy = np.empty_like (Y)
                yy[:] = Y
                yy[np.where(mask==0)] = 0
                yyy=np.abs(yy-Y[YNei_Idx,XNei_Idx])
                yyy[yyy>(2*Max_Diam)]=0
                Nei_HYP_alt = np.sqrt(xxx**2 + yyy**2)
                # I have to invert the angles from 0 - I have the operation backwards somewhere....
                # Several lines of code that follow were inspired by the following post at Stack Overflow
                # https://stackoverflow.com/questions/33976911/
                phi_alt = np.arccos(0) - np.arccos(Nei_HYP_alt / (Nei_grains[nnn] * 0.5))
                phi_alt[np.where(mask==0)] = 0
                HElev = np.cos(phi_alt) * (Nei_grains[nnn] * 0.5)
                HElev[np.where(mask==0)] = 0
                Surf = np.sum(np.stack((Surf, HElev), axis=2), axis=2)
                #vor = Voronoi(xx,yy)
#                if nnn == 1:#(Nei_len - 1):
#                    cmap = plt.cm.binary
#                    #fig, ax = plt.subplots()
#                    fig = plt.figure(10)
#                    ax = fig.add_subplot(1, 1, 1, aspect='equal')
#                    pos = ax.imshow(Surf, extent=(XNei_array[0], XNei_array[-1], YNei_array[-1], YNei_array[0]), cmap=cmap)
#                    fig.colorbar(pos, ax=ax)
                # Surf = np.sum(np.place(Surf, mask, Hemi_elev_alt)
#                mask_flat = np.nonzero(mask)

                # ax.set_axis_off()
                #
#                X_dist = np.abs(X[mask] - X[YNei_Idx,XNei_Idx])
#                Y_dist = np.abs(Y[mask] - Y[YNei_Idx,XNei_Idx])
#                Nei_HYP = np.sqrt(X_dist**2 + Y_dist**2)
#                # Angle of inclination
#                phi = np.arccos(Nei_HYP / (Nei_grains[nnn] * 0.5))
#                # Elevation on the hemisphere
#                Hemi_elev = np.cos(phi) * (Nei_grains[nnn] * 0.5)
#                Elev_Surf[Elev_XY] = Hemi_elev

    del n
    # Store the new x-location for the entrained particles
    Grains_sort[0, Rand_ELocs_Store] = Hop_locx
    Grains_sort[1, Rand_ELocs_Store] = Hop_locy
# END SECTION
###############################################################################
# STEP SIX: PRINT RESULTS TO SCREEN AND WRITE PERTINENT RESULTS TO FILE
    print("--- %s seconds ---" % (time.time() - Start_time))
#   fig = plt.figure(2)
#   ax = fig.add_subplot(1, 1, 1, aspect='equal')
#   ax.set_xlim((-5, x_max + 5))
#   ax.set_xticks([0, 40, 80, 120, 160, 200])
#   ax.set_ylim((-5, y_max+5))
#   resolution = 50  # the number of vertices
#   # Radius_Array = np.asarray((Diameter / 2), dtype=float)
#   # Radius_Array = Radius_Array.reshape(-1)
#   # XCenter = np.reshape(CenterCoord[0, :], (1, len(CenterCoord[0])))
#   # XCenter = XCenter.reshape(-1)
#   # YCenter = np.reshape(CenterCoord[1, :], (1, len(CenterCoord[0])))
#   # YCenter = YCenter.reshape(-1)
#   Radius_Array = np.asarray((Particles_sort[2,:] / 2), dtype=float)
#   Radius_Array = Radius_Array.reshape(-1)
#   XCenter = np.reshape(Particles_sort[0, :], (1, len(Particles_sort[0])))
#   XCenter = XCenter.reshape(-1)
#   YCenter = np.reshape(Particles_sort[1, :], (1, len(Particles_sort[0])))
#   YCenter = YCenter.reshape(-1)
#   plt.rcParams['image.cmap'] = 'gray'
#   # This method of plotting circles comes from Stack Overflow questions\32444037
#   # Note that the patches won't be added to the axes, instead a collection will.
#   patches = []
#   for x1, y1, r in zip(XCenter, YCenter, Radius_Array):
#       circle = Circle((x1, y1), r)
#       patches.append(circle)
##   colors = np.random.uniform(490,510,len(CenterCoord[0]))
#   p = (PatchCollection(patches, color = "#BDBDBD",
#                        alpha=0.9, linewidths=(0, )))
##   p.set_array(colors)
##   p.set_clim([450, 500])
#   ax.add_collection(p)
#   # This method of plotting circles comes from Stack Overflow questions\32444037
#   # Note that the patches won't be added to the axes, instead a collection will.
#   ax = fig.add_subplot(1, 1, 1, aspect='equal')
#   ax.set_xlim((-5, x_max + 5))
#   ax.set_ylim((-5, y_max + 5))
#   resolution = 50  # the number of vertices
    XCenter_entrain = np.asarray(Entrain_locsx).ravel()
    YCenter_entrain = np.asarray(Entrain_locsy).ravel()
    E_radius_entrain = np.asarray(Entrain_diam).ravel() / 2
    Hop_XCenter = np.asarray(Hop_locx).ravel()
    Hop_YCenter = np.asarray(Hop_locy).ravel()
#   patches = []
#   for x1, y1, r in zip(XCenter_entrain, YCenter_entrain, E_radius_entrain):
#       circle = Circle((x1, y1), r)
#       patches.append(circle)
#   #colors = 1000*np.random.rand(len(CenterCoord[0]))
#   p1 = (PatchCollection(patches, color = "#63B8FF",
#                        alpha=0.9, linewidths=(0, )))#cm.Reds
#   ax.add_collection(p1)
#   plt.show()
#   fig.savefig('.\ScriptTest\TotalBedMotions0000%d.png' % plot_step, format='png', dpi=600)
#
#   plot_step = plot_step + 1

    fig = plt.figure(3)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_xlim((-1+LxB, HxB + 1))
    #ax.set_xlim((-5, x_max + 5))
    #ax.set_xticks([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])
    ax.set_ylim((-1+LyB, HyB+1))
    #ax.set_ylim((-5, y_max+5))
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', direction='out')
    ax.tick_params(axis='y', which='minor', direction='out')
    resolution = 50  # the number of vertices
    # xLH = (XNei_array[0],XNei_array[1],XNei_array[1],XNei_array[0],XNei_array[0])
    # yLH = (YNei_array[0],YNei_array[0],YNei_array[1],YNei_array[1],YNei_array[0])
    # ax.plot(xLH,yLH, '-o', c='r')
    Radius_Array = np.asarray((Grains_sort[2,:] / 2), dtype=float)
    Radius_Array = Radius_Array.reshape(-1)
    XCenter = np.reshape(Grains_sort[0, :], (1, len(Grains_sort[0])))
    XCenter = XCenter.reshape(-1)
    YCenter = np.reshape(Grains_sort[1, :], (1, len(Grains_sort[0])))
    YCenter = YCenter.reshape(-1)
    plt.rcParams['image.cmap'] = 'gray'
    # This method of plotting circles comes from Stack Overflow questions\32444037
    # Note that the patches won't be added to the axes, instead a collection will.
    patches = []
    for x1, y1, r in zip(XCenter, YCenter, Radius_Array):
        circle = Circle((x1, y1), r)
        patches.append(circle)
    #colors = 1000*np.random.rand(len(CenterCoord[0]))
    p2 = (PatchCollection(patches, color = "#C1C1C1",
                         alpha=0.75, linewidths=(0, )))#cm.binary
    # p.set_array(colors)
    # p.set_clim([5, 950])
    ax.add_collection(p2)
    # p.set_array(colors)
    # p.set_clim([5, 950])
    patches = []
    for x2, y2, r2 in zip(Hop_XCenter, Hop_YCenter, E_radius_entrain):
        circle = Circle((x2, y2), r2)
        patches.append(circle)
    # colors = 1000*np.random.rand(len(CenterCoord[0]))
    p3 = (PatchCollection(patches, color = "#5B5B5B",
                         alpha=0.75, linewidths=(0, )))#cm.Blues
    # p1.set_array(colors)
    # p1.set_clim([5, 950])
    ax.add_collection(p3)
    black = (0,0,0,1)
    Ex = np.reshape(Entrain_locsx, (1, len(Entrain_locsx[0])))
    Ex = Ex.reshape(-1)
    Ey = np.reshape(Entrain_locsy, (1, len(Entrain_locsy[0])))
    Ey = Ey.reshape(-1)
    patches = []
    for x3, y3, r3 in zip(Ex, Ey, E_radius_entrain):
        circle = Circle((x3, y3), r3)
        patches.append(circle)
    p4 = (PatchCollection(patches, edgecolors= ('deepskyblue',),
                          facecolors= 'None', linewidths=(1.25, )))#cm.Blues
    ax.add_collection(p4)
    if np.sum([Grains_overlap]) > 0:
        ax.scatter(XNei, YNei,s=0.01, c='b',zorder=1)
    plt.show()
    fig.savefig('./ScriptTest/TotalBedMotionsSqPts0000%d.png' % plot_step, format='png', dpi=600)
    if np.sum([Grains_overlap]) > 0:
        cmap = plt.cm.binary
        fig = plt.figure(4)
        ax = fig.add_subplot(1, 1, 1, aspect='equal')
        pos = ax.imshow(Surf, extent=(XNei_array[0], XNei_array[-1], YNei_array[-1], YNei_array[0]), cmap=cmap)
        fig.colorbar(pos, ax=ax)
    plot_step = plot_step + 1

    ##
    # This method of plotting circles comes from Stack Overflow questions\32444037
   # Note that the patches won't be added to the axes, instead a collection will.
#   fig = plt.figure(3)
#   ax = fig.add_subplot(1, 1, 1, aspect='equal')
#   ax.set_xlim((-5, x_max + 5))
#   ax.set_xticks([0, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400, 440, 480, 520, 560, 600, 640, 680, 720, 760, 800])
#   ax.set_ylim((-1, y_max+1))
#   resolution = 50  # the number of vertices
#   XCenter_entrain = np.asarray(Entrain_locsx).ravel()
#   YCenter_entrain = np.asarray(Entrain_locsy).ravel()
#   E_radius_entrain = np.asarray(E_radius).ravel()
#   Hop_XCenter = np.asarray(Hop_locx).ravel()
#   Hop_YCenter = np.asarray(Hop_locy).ravel()
#   patches = []
#   for x1, y1, r in zip(XCenter_entrain, YCenter_entrain, E_radius_entrain):
#       circle = Circle((x1, y1), r)
#       patches.append(circle)
#   colors = 1000*np.random.rand(len(CenterCoord[0]))
#   p = (PatchCollection(patches, cmap=cm.Reds,
#                        alpha=0.9, linewidths=(0, )))
#   p.set_array(colors)
#   p.set_clim([5, 950])
#   patches = []
#   for x2, y2, r2 in zip(Hop_XCenter, Hop_YCenter, E_radius_entrain):
#       circle = Circle((x2, y2), r2)
#       patches.append(circle)
#   colors = 1000*np.random.rand(len(CenterCoord[0]))
#   p1 = (PatchCollection(patches, cmap=cm.Blues,
#                        alpha=0.9, linewidths=(0, )))
#   p1.set_array(colors)
#   p1.set_clim([5, 950])
#   ax.add_collection(p)
#   ax.add_collection(p1)
#   plt.show()
   # Need to insert some commands to write files to subdirectory.
   # os.makedirs('\\ScriptTest')
   #fig.savefig('.\ScriptTest\Entrained_Particles_Step_%d.pdf' % Step_2, format='pdf', dpi=2400)
   # Save initial results for archiving and plotting.
#   np.save('.\ScriptTest\XCenter_Entrain_Step_%d' % Step_2, XCenter_entrain)
#   np.save('.\ScriptTest\YCenter_Entrain_Step_%d' % Step_2, YCenter_entrain)
#   np.save('.\ScriptTest\Radius_Entrain_Step_%d' % Step_2, E_radius_entrain)
#   np.save('.\ScriptTest\XHop_Entrain_Step_%d' % Step_2, Hop_XCenter)
#   np.save('.\ScriptTest\YHop_Entrain_Step_%d' % Step_2, Hop_YCenter)
#   np.save('.\ScriptTest\Travel_Time_Step_%d' % Step_2, T_p)
#   np.save('.\ScriptTest\Hop_Distance_Step_%d' % Step_2, L_x)

    if Step_2 == Loops:

       my_array = np.concatenate([np.array(i) for i in L_x_Store],axis=1).reshape(-1)
       fig = plt.figure(10)
       bins = np.arange(0, 20, 1) # fixed bin size
       plt.hist(my_array, bins=bins, density=True, color="#3F5D7D", edgecolor='black')
       plt.title('Histogram of Hop Distance (n = 32,000)')
       plt.xlabel('Hop Distance (mm)')
       plt.ylabel('Count')
       plt.show()
       fig.savefig('./ScriptTest/HopDistanceHistogram.pdf', format='pdf', dpi=2400)

       range1 = 0 + (Bed_sampreg - 1)
       range2 = (Loops * Bed_sampreg)
       range3 = 1
       xx = np.arange(0, Loops, 1) #range1, range2, Bed_sampreg)
       for n in range(0, len(xx)):
           idx = xx[n]
           Flux = Nu_out_Store[(Bed_sampreg) - 1, idx]
           Nu_out_ds.extend([Flux])
           if n == (len(xx)-1):
               del Nu_out_ds [0]
       Nu_out_DSCumSum = np.cumsum(Nu_out_ds)
       del n

#       while range3 < (Bed_sampreg - 1):
#           xx = xx = np.arange(range3, range2, Bed_sampreg)
#           for n in range(range3, range2, Bed_sampreg):
#               idx = xx[n]
#               Flux = Nu_out_Store[idx]
#               Nu_out_ds.extend([Flux])
#               if n == (len(xx)-1):
#                  del Nu_out_ds [0]

       fig = plt.figure(11)
       ax = fig.add_subplot(1, 1, 1)
       bins = np.arange(0, 10, 1) # fixed bin size
       ax.set_xlim((-1, max(bins)+1))
       ax.set_xticks([-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

       hist, bin_edges = np.histogram(Nu_out_ds, bins = bins)
       EGrains_total = np.sum(hist * bin_edges[0:-1])
       plt.hist(Nu_out_ds, bins=bins, density=False, color="#3F5D7D", edgecolor='black')
       plt.title('Histogram of Particle Flux at Downstream Boundary')
       plt.xlabel('Flux (# of particles)')
       plt.ylabel('Fraction')
       plt.show()
       fig.savefig('./ScriptTest/FluxDownstreamBoundary.pdf', format='pdf', dpi=2400)

       Time = np.arange(1, Loops + 1, 1)

       fig = plt.figure(20)
       #fig = plt.figure()
       ax1 = fig.add_subplot(111)
       #fig, ax1 = plt.subplots()

       #ax.set_xlim((-1, Loops + 1))
       ax1.plot(Time, Nu_out_Store[(Bed_sampreg - 1), :])
       plt.title('Timeseries of Particle Flux at Downstream Boundary')
       ax1.set_xlabel('Numerical Step')
       ax1.set_ylabel('Particle Flux')
       ax1.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8])
       ax2 = ax1.twinx()

       ax2.plot(Time, Nu_out_DSCumSum, 'r.')
       ax2.set_ylabel('Particle Flux Cumulative Sum', color='r', rotation=270, labelpad=15)
       ax2.tick_params('y', colors='r')

       fig.tight_layout()
       plt.show()
       fig.savefig('./ScriptTest/FluxDownstreamBoundary_2YAx.png', format='png', dpi=1200)

#       for n in range(0, Bed_sampreg, 1):
#           #for nn in range(0, Loops, 1):
#           plt.plot(Nu_out_Store[n, :])
#               plt.title('A tale of 2 subplots')
#               plt.ylabel('Damped oscillation')

       print("--- loop %d out of %d total Loops ---" % (Step_2, Loops))
       print("--- %s total seconds ---" % (time.time() - Start_time))
       print("--- SIMULATION FINISHED ---")
       break

    else:
       #Particles_sort = np.place(Particles_sort[0,:], Rand_ELocs_Store, Hop_locx)
       print("--- loop %d out of %d total Loops ---" % (Step_2, Loops))
       print("--- %s total particles in domain ---" % (Step_1))
       print("--- %s total particles in domain after loop ---") % (len(Grains_sort[0]))
       Step_2 = Step_2 + 1
# Crcle area.
# END OF CODE