# -*- coding: utf-8 -*-
# """
# Created on Thu May 24 12:25:15 2018
# @author: schartrand
# """
import math
import random
import time
import os
import numpy as np
import matplotlib.cm as cm

from numba import jit
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

start_time = time.time()
#
# BEGIN CODE
# This code was developed by Shawn Chartrand and David Furbish, March 2018.
# The code randomly populates a rectangular area of specified dimension with
# circles of random and nonuniform size. Placement of circles at step S
# does not depend on prior placements. Placement of circles stops when a
# specified packing fraction (circle volume/total volume) is reached.
# The code is the basis of Markov birth-death simulations of bedload transport
# in rivers under relatively low transport rates. The code is split into two
# parts. Part 1 builds the random streambed surface. Part 2 entrains, deposits
# and moves particles along the streambed surface.
###############################################################################
###############################################################################
# PART 1
#
# STEPS ONE, TWO AND THREE: DEFINE VARIABLES AND BUILD THE MODEL DOMAIN.
#
# STEP ONE: DEFINE DATA STRUCTURES, INITIALIZE SOME VARIABLES.
# DATA LISTS
# Circle area.
Area = np.zeros([1, 1], dtype=int, order='F')
# Total area of placed circles. Initialize to 0.
AreaTotal = 0
# Coordinates of placed circles.
Coordinates = np.zeros([4, 1], dtype=int, order='F')
# Center coord. of placed circles: r1:x,r2:y.
CenterCoordinates = np.zeros([2, 1], dtype=int, order='F')
# Final placed circle diameters.
Diameter = np.zeros([1, 1], dtype=int, order='F')
# Open coordiantes within rectangular area.
OpenCLStore = np.array([])
# COUNTER AND VARAIBLES
# Defines size of colormap.
color_iter = 1000
# Milestones used to display calculation progress.
milestones = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 100]
# Initialize the step counter for part 1.
step_1 = 1
# Circle perimeter coordinate calculations.
theta = np.arange(0, (2 * math.pi), 0.01)
# Define origin for plotting.
XCoordinates_Orig = 0
YCoordinates_Orig = 0
# END SECTION
#
# STEP TWO: BUILD THE VIRTUAL DOMAIN WITH ALL NECESSARY DETAILS
# KEY INPUT PARAMETER VALUES
# Packing density of the bed; ranges from > 0 to <~0.70
Pack = 0.52
# Length of the domain in the streamwise direction in millimeters
x_max = 200
# Length of the domain in the cross-stream direction in millimeters
y_max = 100
# Spacing nodes
Grid_Spacing = 2.0
# Size of largest grain relative to the cross-stream length
Size_Frac = 0.06
# Time threshold to restart script
time_threshold = 10
# DEPENDENT CALCULATIONS
# Total virtual streambed area.
Total_VArea = x_max * y_max
# Grid spacing in x and y.
x_min_incre = Grid_Spacing
y_min_incre = Grid_Spacing
# Max grain size fraction of height.
CircDiam_Fac = 1 / Size_Frac
# Specify the max circle diameter.
Max_CircleDiameter = y_max / CircDiam_Fac
# END SECTION
#
# STEP THREE: BUILD MATRIX OF COORDINATE ADDRESSES FOR CENTER OF CIRCLES
# Build x array for grid.
x_array = np.arange(x_min_incre, (x_max+1), x_min_incre, dtype=float)
# Build y array for grid.
y_array = np.arange(y_min_incre, (y_max+1), y_min_incre, dtype=float)
# Length of array dim 1 and dim 2
length1 = int(len(y_array))
length2 = int(len(x_array))
# For logical mask.
XYCoordinateslogical = np.zeros((length1, length2), dtype=int)
# Length of matrix dim 1 and dim 2
Dim1 = len(XYCoordinateslogical)
Dim2 = len(XYCoordinateslogical[0])
# This defines the logical array length for the mask.
SIZE = Dim1 * Dim2
# This variable is used for the logical mask.
idx = np.linspace(0, SIZE-1, SIZE, dtype=int)
# BUILD A COORDINATE MATRICES.
YCoordinates = np.repeat(y_array[None], length2, axis=0).transpose()
XCoordinates = np.repeat(x_array[None], length1, axis=0)
# END SECTION
###############################################################################
###############################################################################
# STEP FOUR: DEFINE TWO FUNCTIONS USED TO PLACE and REVISE CIRCLE DIAMETERS


def fu_radcheck(x_Center, y_Center, radius):
    # These operations are used to check if the randomly placed
    # circle falls outside the rectangular area
    x_circ_maxcheck = x_max - (x_Center + radius)
    x_circ_mincheck = (x_Center - radius) - XCoordinates_Orig
    y_circ_maxcheck = y_max - (y_Center + radius)
    y_circ_mincheck = (y_Center - radius) - YCoordinates_Orig
    circ_array = np.arange(0, 4)
    # Build array to check the cumulative status of overlap
    circ_check = (np.reshape(np.asarray([x_circ_maxcheck, x_circ_mincheck,
                                         y_circ_maxcheck, y_circ_mincheck]),
                             (-1)))
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
    # Calculate the four circle edge coordinates.
    x1 = x_Center - radius
    x2 = x1 + (radius * 2)
    y1 = y_Center - radius
    y2 = y1 + (radius * 2)
    coord = [x1, x2, y1, y2]
    if step_1 == 1:
        newCircle_Found = 1 
    else:
        newCircle_Found = 0 
    return radius, coord, newCircle_Found

def fu_diamrevise(x_Center, y_Center, radius):
    # The following operations determines if the test circle overlaps
    # with existing circles, but only to the resolution of the
    # underlying grid. The outcome is a logical matrix. Values of 1
    # mean that grid point overalps between the test circle and an
    # existing circle.
    dist_check = ((CenterCoordinates[0, :] - x_Center) ** 2 +
                  (CenterCoordinates[1, :] - y_Center) ** 2) ** 0.5
    radius_check = (Diameter / 2) + radius
    radius_overlap = (radius_check >= dist_check)
    # Now actually check the overalp condition.
    if np.sum([radius_overlap]) == 0:
        # The new circle does not overlap so proceed.
        newCircle_Found = 1
    elif np.sum([radius_overlap]) == 1:
        # The new circle overlaps with one other circle
        overlap = (np.arange(0, len(radius_overlap[0]),
                             dtype=int).reshape(1, len(radius_overlap[0])))
        overlap_logix = (radius_overlap == 1)
        idx_true = overlap[overlap_logix]
        radius = dist_check[idx_true] - (Diameter[0, idx_true] / 2)
        # New circle radius can be negative when inside another circle.
        if radius > 0:
            newCircle_Found = 1
        else:
            newCircle_Found = 0
    elif np.sum([radius_overlap]) >= 2:
        # Overlap with two or more circles.
        # Set newCircle flag to zero so the step is repeated.
        newCircle_Found = 0
    return radius, newCircle_Found
 # END SECTION
###############################################################################
###############################################################################
# STEP FIVE: RANDOMLY PLACE CIRCLES IN RECTANGULAR AREA
# Circles are randomly placed until the total area of placed circles
# equals or is greater than the packing threshold.


while (AreaTotal / Total_VArea) < Pack:
    # Initial calculations for all steps.
    # Random number to assign a color.
    c_rand = random.randint(1, color_iter)
    # Flatten the logical matrix to an array
    coord_mat_transform = XYCoordinateslogical.flatten('F')
    # Random circle diameter
    diameter_Test = random.randint(Grid_Spacing, Max_CircleDiameter)
    # Index used to track open coordinate
    test = (coord_mat_transform == 0)
    # Index used to track open coord.
    idx2 = idx[test]
    # Array length.
    openCoordLength = len(idx2)
    # Random number.
    randCoord = random.randint(1, (openCoordLength-1))
    # Choose open coord. with random number.
    idx3 = idx2[randCoord]
    # Specify random x and y center coordinate.
    x_Center = XCoordinates[np.unravel_index([idx3], [length1, length2],
                                             order='F')]
    y_Center = YCoordinates[np.unravel_index([idx3], [length1, length2],
                                             order='F')]
    CenterPosition_data = [x_Center, y_Center]
    radius = diameter_Test / 2
    # CALL FUNCTION
    radius, coord, newCircle_Found = fu_radcheck(x_Center, y_Center, radius)
    # Check to see if radius needs revision for steps greater than the first
    if step_1 != 1:
        # CALL FUNCTION
        radius, newCircle_Found = fu_diamrevise(x_Center, y_Center, radius)    
    # Only write data if a circle was placed.
    if newCircle_Found == 1:
        # Write data to variables
        # coord = np.array([x1, x2, y1, y2])
        # This function moves circles in a rectangle to available grid points.
        mask = (((XCoordinates - x_Center) ** 2 + (YCoordinates - y_Center)
                 ** 2) ** 0.5 <= radius)
        # Track which grid points are available for circles and
        # which are not and pass to the next step in the loop.
        XYCoordinateslogical[np.where(mask == 1)] = 1
        # Write data to variables.
        Circle_Area = (np.reshape(np.asfortranarray(math.pi * (radius ** 2),
                                                    dtype=float), (1, +1)))
        Area = np.hstack((Area, Circle_Area))
        AreaTotal = np.sum(Area)
        Coordinates = np.hstack((Coordinates, coord))
        CenterCoordinates = np.hstack((CenterCoordinates, CenterPosition_data))
        Diameter = np.hstack((Diameter, np.reshape(radius * 2, (1, +1))))
        if step_1 == 1:
            Area = np.delete(Area, 0, 1)
            Coordinates = np.delete(Coordinates, 0, 1)
            CenterCoordinates = np.delete(CenterCoordinates, 0, 1)
            Diameter = np.delete(Diameter, 0, 1)
        # Advance counter.
        step_1 = step_1 + 1
        # Calculate the percent completion of simulation and print to console.
        percentage_complete = (100.0 * (AreaTotal / Total_VArea)/Pack)
        while len(milestones) > 0 and percentage_complete >= milestones[0]:
            print "{}% complete".format(milestones[0])
            # remove that milestone from the list
            milestones = milestones[1:]       
    # This operation restarts the script if the execution time is too long.
    # I think there is a bug but I cannot figure out what it is. Filling the
    # rectangular area stalls out at 95 or 99% about 1 out of every 4 runs.
    current_time = time.time()
    if (current_time - start_time) > time_threshold:
        start_time = time.time()
        step_1 = 1
        percentage_complete = 0
        # Circle area.
        Area = np.zeros([1, 1], dtype=int, order='F')
        # Total area of placed circles. Initialize to 0.
        AreaTotal = 0
        # Coordinates of placed circles.
        Coordinates = np.zeros([4, 1], dtype=int, order='F')
        # Center coord. of placed circles: r1:x,r2:y.
        CenterCoordinates = np.zeros([2, 1], dtype=int, order='F')
        # Final placed circle diameters.
        Diameter = np.zeros([1, 1], dtype=int, order='F')
        # For logical mask.
        XYCoordinateslogical = np.zeros((length1, length2), dtype=int)
        # Milestones used to display calculation progress.
        milestones = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 100]
        print("--- RESTARTING LOOP ---")
            
# END SECTION
###############################################################################
###############################################################################
# STEP SIX: PRINT RESULTS TO SCREEN AND WRITE PERTINENT RESULTS TO FILE
print("--- %s seconds ---" % (time.time() - start_time))
# STEP SIX: PLOT THE RESULTS
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlim((-1, x_max+1))
ax.set_ylim((-1, y_max+1))
resolution = 50  # the number of vertices
Radius_Array = np.asarray((Diameter / 2), dtype=float)
Radius_Array = Radius_Array.reshape(-1)
XCenter = np.reshape(CenterCoordinates[0, :], (1, len(CenterCoordinates[0])))
XCenter = XCenter.reshape(-1)
YCenter = np.reshape(CenterCoordinates[1, :], (1, len(CenterCoordinates[0])))
YCenter = YCenter.reshape(-1)
plt.rcParams['image.cmap'] = 'gray'
# This method of plotting circles comes from Stack Overflow questions/32444037
# Note that the patches won't be added to the axes, instead a collection will.
patches = []
for x1, y1, r in zip(XCenter, YCenter, Radius_Array):
    circle = Circle((x1, y1), r)
    patches.append(circle)
colors = 1000*np.random.rand(len(CenterCoordinates[0]))
p = (PatchCollection(patches, cmap=cm.binary,
                     alpha=0.9, linewidths=(0, )))
p.set_array(colors)
p.set_clim([5, 950])
ax.add_collection(p)
plt.show()
# Need to insert some commands to write files to subdirectory.
# os.makedirs('\\ScriptTest')
# fig.savefig('.\ScriptTest\Starting_Bed.pdf', format='pdf', dpi=2400)
fig.savefig('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/Starting_Bed.pdf', format='pdf', dpi=2400)
# Save initial results for archiving and plotting.
np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/XCenter_Initial', XCenter)
np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/YCenter_Initial', YCenter)
np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/Radius_Array_Initial', Radius_Array)
# END OF PART 1 VIRTUAL BED CREATION
###############################################################################
###############################################################################
# PART 2
# 
# STEP ONE: DEFINE DATA STRUCTURES, INITIALIZE SOME VARIABLES.
# DATA LISTS
# Upstream particle supply rate (particles/t).
Nu_in = np.zeros([1, 1], dtype=int, order='F')
# Downstream particle emigration rate (particles/t).
Nu_out = np.zeros([1, 1], dtype=int, order='F')
# Particle birth or entrainment rate within the control volume (particle/t).
# Set birth rate as a constant value for trial runs
Lambda_1 = 3
# Entrainment events per unit area of the bed and unit time (events/t)
k_events = np.zeros([1, 1], dtype=int, order='F')
# Particle death or deposition rate within the control volume (1/t).
Sigma = np.zeros([1, 1], dtype=int, order='F')
# Particle travel time (t).
T_p = np.zeros([1, 1], dtype=int, order='F')
# Particle travel time minimum (t).
T_pmin = 0
# Particle travel time maximum (t).
T_pmax = 1.0
# Particle hop distance (L).
L_x = np.zeros([1, 1], dtype=int, order='F')
L_x_Store = []
L_x_Store_temp = np.zeros([1, 1], dtype=int, order='F')
# Parameter to store coordinate and diameter data from step 1
Particle_details = np.vstack((CenterCoordinates, Diameter))
# Initialize the step counter for part 2.
step_2 = 0
# Temporary loop completion criteria
loops = 10
# Downstream boundary condition x_location. All x_locations downstream of 
# x_max-DS_Bound will have a fixed surface particle condition for all numerical
# time steps [L].
DS_Bound = 10
DS_Fix = x_max - DS_Bound
# END SECTION
#
# STEP TWO: DIVIDE THE BED INTO SAMPLING REGIONS AND BUILD SAMPLING ARRAY
# Number of bed sampling regions within the control volume V_c
Bed_sampreg = 10
# Index to hold boundaries between subsampling regions to facilitate random
# sampling from within each region.
SubSampl_idx = np.zeros([1, Bed_sampreg], dtype=int, order='F')
# Bed sampling boundaries in the x-direction
BS_boundaries = XCoordinates_Orig + (x_max / Bed_sampreg) * np.arange(1, Bed_sampreg + 1, dtype=int)
###############################################################################
###############################################################################
# STEP TWO: DEFINE A SEARCH FUNCTION AND ONE TO SPECIFY PARTICLE KINEMATICS


@jit(nopython=True)
def fu_search(item, vec):
    """return the index of the first occurence of item in vec"""
    for i in xrange(len(vec)):
        if vec[i] >= item:
            return i
    return -1

#def fu_kinematics(x_Center, y_Center, radius):
#    # k_events is the entrainment events per unit area of the bed
#    k_events = np.random.poisson(Lambda_1,None)
#    # Total number of entrained particles not counting immigrants across V_c
#    E_particles = k_events * Bed_sampreg
#    # Calculate the travel time as a randomly sampled variable from a uniform
#    # distribution constrained by Fathel et al., 2015.
#    T_p = (np.random.uniform(T_pmin,T_pmax,E_particles).reshape(1, E_particles))
#    # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
#    # For now I am multuplying by 10 so units are consistent ()
#    L_x = np.round((T_p ** 2) * 10, 2)
#    # Figure out where to randomly entrain particles within the control volume.
#    CenterCoord_sort = CenterCoordinates[:, CenterCoordinates[0].argsort()]
#    # Initialize a variable to store entrainment location information
#    Entrain_locs = np.zeros([1, k_events], dtype=int, order='F')
#        
#    # Loop to boundaries between subsampling regions
#    for n in range(0, Bed_sampreg):
#        # Find indices within sort coordinates array for subsampling
#        SubSampl_idx[0,n] = search(BS_boundaries[n],CenterCoord_sort[0,:])
#        # Specify random entrainment sampling indices
#        if n == 0:
#            locs = (np.random.randint(0,SubSampl_idx[0,n],k_events).reshape(1, k_events))
#        else:
#            locs = (np.random.randint(SubSampl_idx[0,n-1],SubSampl_idx[0,n],k_events).reshape(1, k_events))
#        # Write the entrainment location data to a variable
#        Entrain_locs = np.vstack((Entrain_locs, locs))
#        # Delete the initialization value b/c it is meaningless
#        if n == 0:
#            Entrain_locs = np.delete(Entrain_locs, (0), axis=0)
#    # This is the destination of partciles after the hop distance is applied.
#    # The result returns the destination in the streamwise coordinate only.        
#    Hop_loc = CenterCoord_sort[0,Entrain_locs] + L_x
#        #if Hop_loc > DS_Fix:
#            
#    return radius, coord, newCircle_Found, Nu_out
###############################################################################
###############################################################################
# STEP TWO: ENTRAIN, MOVE AND DEPOSIT PARTICLES IN CONTROL VOLUME
# Particles are entrained, moved and deposited until a steady-state condition 
# is met.


while step_2 < loops:
#    radius, coord, newCircle_Found, Nu_out = fu_kinematics(x_Center, y_Center, radius)
    
    # k_events is the entrainment events per unit area of the bed
    k_events = np.random.poisson(Lambda_1,None)
    if k_events == 0:
        k_events = np.random.poisson(Lambda_1,None)
    # Total number of entrained particles not counting immigrants across V_c
    E_particles = k_events * Bed_sampreg
#    # Calculate the travel time as a randomly sampled variable from a uniform
#    # distribution constrained by Fathel et al., 2015.
#    T_p = (np.random.uniform(T_pmin,T_pmax,E_particles).reshape(1, E_particles))
#    # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
#    # For now I am multuplying by 10 so units are consistent ()
#    L_x = np.round((T_p ** 2) * 10, 2)
    # Figure out where to randomly entrain particles within the control volume.
    Particle_details_sort = Particle_details[:, Particle_details[0].argsort()]
    # Initialize two variables to store entrainment location information
    Entrain_idx = np.zeros([1, k_events], dtype=int, order='F')
    Entrain_locsx = np.zeros([1, k_events], dtype=int, order='F')
    Entrain_locsy = np.zeros([1, k_events], dtype=int, order='F')
    T_p = np.zeros([1, k_events], dtype=int, order='F')
    L_x = np.zeros([1, k_events], dtype=int, order='F')
    E_diameter = np.zeros([1, k_events], dtype=int, order='F')
        
    # Loop to boundaries between subsampling regions
    for n in range(0, Bed_sampreg):
        # Find indices within sort coordinates array for subsampling
        SubSampl_idx[0,n] = fu_search(BS_boundaries[n],Particle_details_sort[0,:])
        # Calculate the travel time as a randomly sampled variable from a uniform
        # distribution constrained by Fathel et al., 2015.
        T_ptemp = (np.random.uniform(T_pmin,T_pmax,k_events).reshape(1, k_events))
        # https://stackoverflow.com/questions/2106503/
        T_ptemp2 = np.log(1-T_ptemp) / (-Lambda_1)
        # Calculate L_x per Fathel et al., 2015 and Furbish et al., 2017.
        # For now I am multuplying by 10 so units are consistent ()
        L_xtemp = np.round((T_ptemp2 ** 2) * 10, 2)
        # Specify random entrainment sampling indices
        if n == 0:
            locs = (np.random.randint(0,SubSampl_idx[0,n],k_events).reshape(1, k_events))
            # Actual entrainment x-coordinates
            locs_2 = Particle_details_sort[0,locs]
            # Actual entrainment y-coordinates
            locs_3 = Particle_details_sort[1,locs]
            # Associated entrained particle diameters
            E_diam = Particle_details_sort[2,locs]
        else:
            locs = (np.random.randint(SubSampl_idx[0,n-1],SubSampl_idx[0,n],k_events).reshape(1, k_events))
            # Actual entrainment x-coordinates
            locs_2 = Particle_details_sort[0,locs]
            # Actual entrainment y-coordinates
            locs_3 = Particle_details_sort[1,locs]
            # Associated entrained particle diameters
            E_diam = Particle_details_sort[2,locs]            
        # Write the entrainment location data to three variables
        Entrain_idx = np.vstack((Entrain_idx, locs))
        Entrain_locsx = np.vstack((Entrain_locsx, locs_2))
        Entrain_locsy = np.vstack((Entrain_locsy, locs_3))
        T_p = np.vstack((T_p, T_ptemp))
        L_x = np.vstack((L_x, L_xtemp))
        E_diameter = np.vstack((E_diameter, E_diam))
        # Delete the initialization value b/c it is meaningless
        if n == 0:
            Entrain_idx = np.delete(Entrain_idx, (0), axis=0)
            Entrain_locsx = np.delete(Entrain_locsx, (0), axis=0)
            Entrain_locsy = np.delete(Entrain_locsy, (0), axis=0)
            T_p = np.delete(T_p, (0), axis=0)
            L_x = np.delete(L_x, (0), axis=0)
            E_diameter = np.delete(E_diameter, (0), axis=0)
    # This is the destination of partciles after the hop distance is applied.
    # The result returns the destination in the streamwise coordinate only.        
    Hop_locx = Entrain_locsx + L_x
    Hop_locy = Entrain_locsy
    step_2 = step_2 + 1
    E_radius = E_diameter / 2
    L_x_Store_temp = np.asarray(L_x).ravel()
    L_x_Store.append(L_x_Store_temp)
    
    flat_list = [item for sublist in L_x_Store for item in sublist]
    for sublist in L_x_Store:
        for item in sublist:
            flat_list.append(item)
            
    bins = np.arange(0, 20, 0.1) # fixed bin size
    plt.xlim([min(L_x_Store)-5, max(L_x_Store)+5])

    plt.hist(L_x_Store, bins=bins, alpha=0.5)
    plt.title('Random Gaussian data (fixed bin size)')
    plt.xlabel('variable X (bin size = 5)')
    plt.ylabel('count')

    plt.show()
        
###############################################################################
## STEP SIX: PRINT RESULTS TO SCREEN AND WRITE PERTINENT RESULTS TO FILE
    print("--- %s seconds ---" % (time.time() - start_time))
    # STEP SIX: PLOT THE RESULTS
    fig = plt.figure(2)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_xlim((-1, x_max+1))
    ax.set_ylim((-1, y_max+1))
    resolution = 50  # the number of vertices
    Radius_Array = np.asarray((Diameter / 2), dtype=float)
    Radius_Array = Radius_Array.reshape(-1)
    XCenter = np.reshape(CenterCoordinates[0, :], (1, len(CenterCoordinates[0])))
    XCenter = XCenter.reshape(-1)
    YCenter = np.reshape(CenterCoordinates[1, :], (1, len(CenterCoordinates[0])))
    YCenter = YCenter.reshape(-1)
    plt.rcParams['image.cmap'] = 'gray'
    # This method of plotting circles comes from Stack Overflow questions/32444037
    # Note that the patches won't be added to the axes, instead a collection will.
    patches = []
    for x1, y1, r in zip(XCenter, YCenter, Radius_Array):
        circle = Circle((x1, y1), r)
        patches.append(circle)
    colors = 1000*np.random.rand(len(CenterCoordinates[0]))
    p = (PatchCollection(patches, cmap=cm.binary,
                         alpha=0.9, linewidths=(0, )))
    p.set_array(colors)
    p.set_clim([5, 950])
    ax.add_collection(p)
    plt.hold(True)
    #plt.show()
    # This method of plotting circles comes from Stack Overflow questions/32444037
    # Note that the patches won't be added to the axes, instead a collection will.
    #fig = plt.figure(3)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_xlim((-1, x_max+1))
    ax.set_ylim((-1, y_max+1))
    resolution = 50  # the number of vertices
    XCenter_entrain = np.asarray(Entrain_locsx).ravel()
    YCenter_entrain = np.asarray(Entrain_locsy).ravel()
    E_radius_entrain = np.asarray(E_radius).ravel()
    Hop_XCenter = np.asarray(Hop_locx).ravel()
    Hop_YCenter = np.asarray(Hop_locy).ravel()
    patches = []
    for x1, y1, r in zip(XCenter_entrain, YCenter_entrain, E_radius_entrain):
        circle = Circle((x1, y1), r)
        patches.append(circle)
    colors = 1000*np.random.rand(len(CenterCoordinates[0]))
    p = (PatchCollection(patches, cmap=cm.Reds,
                         alpha=0.9, linewidths=(0, )))
    p.set_array(colors)
    p.set_clim([5, 950])
    plt.hold(True)
    patches = []
    for x2, y2, r2 in zip(Hop_XCenter, Hop_YCenter, E_radius_entrain):
        circle = Circle((x2, y2), r2)
        patches.append(circle)
    colors = 1000*np.random.rand(len(CenterCoordinates[0]))
    p1 = (PatchCollection(patches, cmap=cm.Blues,
                         alpha=0.9, linewidths=(0, )))
    p1.set_array(colors)
    p1.set_clim([5, 950])
    ax.add_collection(p)
    plt.hold(True)
    ax.add_collection(p1)
    plt.show()
    fig.savefig('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/Entrained_TotalBed_Step_%d.pdf' % step_2, format='pdf', dpi=2400)
    ##
    ##
    # This method of plotting circles comes from Stack Overflow questions/32444037
    # Note that the patches won't be added to the axes, instead a collection will.
    fig = plt.figure(3)
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    ax.set_xlim((-1, x_max+1))
    ax.set_ylim((-1, y_max+1))
    resolution = 50  # the number of vertices
    XCenter_entrain = np.asarray(Entrain_locsx).ravel()
    YCenter_entrain = np.asarray(Entrain_locsy).ravel()
    E_radius_entrain = np.asarray(E_radius).ravel()
    Hop_XCenter = np.asarray(Hop_locx).ravel()
    Hop_YCenter = np.asarray(Hop_locy).ravel()
    patches = []
    for x1, y1, r in zip(XCenter_entrain, YCenter_entrain, E_radius_entrain):
        circle = Circle((x1, y1), r)
        patches.append(circle)
    colors = 1000*np.random.rand(len(CenterCoordinates[0]))
    p = (PatchCollection(patches, cmap=cm.Reds,
                         alpha=0.9, linewidths=(0, )))
    p.set_array(colors)
    p.set_clim([5, 950])
    plt.hold(True)
    patches = []
    for x2, y2, r2 in zip(Hop_XCenter, Hop_YCenter, E_radius_entrain):
        circle = Circle((x2, y2), r2)
        patches.append(circle)
    colors = 1000*np.random.rand(len(CenterCoordinates[0]))
    p1 = (PatchCollection(patches, cmap=cm.Blues,
                         alpha=0.9, linewidths=(0, )))
    p1.set_array(colors)
    p1.set_clim([5, 950])
    ax.add_collection(p)
    plt.hold(True)
    ax.add_collection(p1)
    plt.show()
    # Need to insert some commands to write files to subdirectory.
    # os.makedirs('\\ScriptTest')
    fig.savefig('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/Entrained_Particles_Step_%d.pdf' % step_2, format='pdf', dpi=2400)
    # Save initial results for archiving and plotting.
    np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/XCenter_Entrain_Step_%d' % step_2, XCenter_entrain)
    np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/YCenter_Entrain_Step_%d' % step_2, YCenter_entrain)
    np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/Radius_Entrain_Step_%d' % step_2, E_radius_entrain)
    np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/XHop_Entrain_Step_%d' % step_2, Hop_XCenter)
    np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/YHop_Entrain_Step_%d' % step_2, Hop_YCenter)
    np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/Travel_Time_Step_%d' % step_2, T_p)
    np.save('/Volumes/ESD-USB/GitHub/MarkovProcess_Bedload/ScriptTest/Hop_Distance_Step_%d' % step_2, L_x)
# Circle area.
# END OF CODE
