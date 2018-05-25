# -*- coding: utf-8 -*-
"""
Created on Thu May 24 12:25:15 2018
@author: schartrand
"""
# First two lines clear variables from memory each time the script is executed.
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import math
import random
# import matplotlib.pyplot as plt
# BEGIN CODE
# This code was developed by Shawn M. Chartrand, March 2018. The code
# randomly populates a rectangular area of specified dimension with 
# circles of random and nonuniform size. The range of circle sizes is 
# is specified througha user interface. Placement of circles at step S 
# does not depend on prior placements. Placement of circles stops when a
# specified packing fraction (circle volume/total volume) is reached.
# A maximum packing fraction of 70% is generally specified. The code is the
# basis of Markov birth-death simulations of bedload transport in rivers
# under low transport rates. I write comments above and to the side of 
# variable names depending on how much documentation is needed. I like to 
# document operations thoroughly. 
#
# STEP ONE: DEFINE DATA STRUCTURES, INITIALIZE SOME VARIABLES.                      # COMMENTS
# DATA LISTS
DebugArray = []                                                                     # De-bugging values.
Area = []                                                                           # Circle area.
AreaTotal = 0                                                                       # Total area of placed circles. Initialize to 0.
OpenCLStore = []                                                                    # Open coordiantes within rectangular area.
Coordinates = []                                                                    # Coordinates of placed circles.
CenterCoordinates = []                                                              # Center coord. of placed circles: r1:x,r2:y.
Diameter = []                                                                       # Final placed circle diameters.
# COUNTER AND VARAIBLES
step = 1;                                                                           # Initialize the step counter.
theta = np.arange(0,(2 * math.pi),0.01)                                             # Circle perimeter coordinate calculations.
color_iter = 1000                                                                   # Defines size of colormap.
XCoordinates_Orig = 0                                                               # Define origin for plotting.
YCoordinates_Orig = 0                                                               # Define origin for plotting.
bound_buffer = 1                                                                    # Buffer for particle placement near boundaries.
milestones = [5, 15, 30, 45, 60, 75, 90, 100]                                       # Milestones used to display calculation progress.
# END SECTION
#
# STEP TWO: BUILD THE VIRTUAL DOMAIN WITH ALL NECESSARY DETAILS 
# PROMPT THE USER
Pack  = 0.7 #float(raw_input('Enter the packing fraction of the bed (0.70 is recommended maximum): '))# Packing fraction of virtual bed
x_max = 64 #float(raw_input('Enter the downstream (x) length of the virtual streambed as an integer (mm):'))# Length of virtual streambed.
y_max = 64 #float(raw_input('Enter the crossstream (y) height of the virtual streambed as an integer(mm):'))# Height of virtual streambed.
Grid_Spacing = 1 #float(raw_input('Enter the grid spacing (mm - minimum circle size set by spacing and 1 mm is recommended):'))    
Size_Frac = 0.5 #float(raw_input('Enter the max grain size as a fraction of the virtual streambed height (0.50 is recommended)):'))        
# DEPENDENT CALCULATIONS
Total_VArea = x_max * y_max                                                         # Total virtual streambed area.
x_min_incre = Grid_Spacing                                                          # Grid spacing in x.
y_min_incre = Grid_Spacing                                                          # Max grain size fraction of height.   
CircDiam_Fac = 1 / Size_Frac                                                        # Factor for max circle diameter.
Max_CircleDiameter = y_max / CircDiam_Fac                                           # Specify the max circle diameter. 
# END SECTION
#
# STEP THREE: BUILD MATRIX OF COORDINATE ADDRESSES FOR CENTER OF CIRCLES
x_array = np.arange(x_min_incre,(x_max+1),x_min_incre)                              # Build x array for grid.
y_array = np.arange(y_min_incre,(y_max+1),y_min_incre)                              # Build y array for grid.
length1 = int(len(y_array))                                                         # Length of array dim 1
length2 = int(len(x_array))                                                         # Length of array dim 2
XYCoordinateslogical = np.zeros((length1,length2),dtype=int)                        # For mask.
Dim1 = len(XYCoordinateslogical)                                                    # Length of matrix dim 1
Dim2 = len(XYCoordinateslogical[0])                                                 # Length of matrix dim 2
SIZE = Dim1 * Dim2                                                                  # This defines the logical array length for the mask.
idx = np.linspace(1, SIZE, SIZE, dtype=int)                                         # This variable is used for the logical mask.
# BUILD A COORDINATE MATRICES.
YCoordinates = np.repeat(y_array[None],length2,axis=0).transpose()
XCoordinates = np.repeat(x_array[None],length1,axis=0)
# END SECTION
#
# STEP FOUR: RANDOMLY PLACE CIRCLES IN RECTANGULAR AREA
# Circles are randomly placed until the total area of placed circles
# equals or is greater than the packing threshold.
while (AreaTotal / Total_VArea) < Pack:
    # Initial calculations for all steps.
    diameter_Test = random.randint(Grid_Spacing,Max_CircleDiameter)                 # Random circle dia.
    c_rand = random.randint(1,color_iter)                                           # Random number to assign a color.
    coord_mat_transform = XYCoordinateslogical.flatten('F')                         # Flatten the matrix to an array
    test = (coord_mat_transform == 0)                                               # Index used to track open coord.
    idx2 = idx[test]                                                                # Index used to track open coord.
    openCoordLength = len(idx2)                                                     # Array length.
    randCoord = random.randint(1,openCoordLength)                                   # Random number.
    idx3 = idx2[randCoord]                                                          # Choose open coord. with random number.
    x_Center = XCoordinates[np.unravel_index([idx3],[length1,length2],order='F')]   # Specify random x coordinate.
    y_Center = YCoordinates[np.unravel_index([idx3],[length1,length2],order='F')]   # Specify random y coordinate.
    CenterPosition_data = ([x_Center,y_Center])                                     # Combine center data for appending.
    radius = diameter_Test / 2                                                      # Assign random circle radius
    # For the first step
    if step == 1:
        # These operations are used to check if the randomly placed
        # circle falls outside the rectangular area
        x_circ_maxcheck = x_max - (x_Center + radius)
        x_circ_mincheck = (x_Center - radius) - XCoordinates_Orig
        y_circ_maxcheck = y_max - (y_Center + radius)
        y_circ_mincheck = (y_Center - radius) - YCoordinates_Orig
        circ_array = [1, 2, 3, 4]
        circ_check = [x_circ_maxcheck, x_circ_mincheck, y_circ_maxcheck, y_circ_mincheck]
        circ_check_logical = np.where(circ_check <= 0)
        # Check to see if the circle radius needs to be changed.
        if sum(circ_check_logical) > 0:                                             # Revise circle radius
            # Use logical indexing to speed things up.
            circ_logidx = np.where(circ_check_logical == 1)
            circ_check_idx = circ_array[circ_logidx]
            circ_check_radrev = min(circ_check(circ_check_idx))
            radrev = radius + circ_check_radrev
            radius = radrev                                                         # Revise circle radius.                     
        # Calculate the four circle coordinates.
        x1 = x_Center - radius                                                      # Edge 1 of vert. tangent line.
        x2 = x1 + (radius * 2)                                                      # Edge 2 of vert. tangent line.
        y1 = y_Center - radius                                                      # Edge 1 of hor. tangent line.
        y2 = y1 + (radius * 2)                                                      # Edge 2 of hor. tangent line.
        coord = [x1, x2, y1, y2]                                                    # Build an array.
        debug_value = 1
    # AFTER THE INITIAL STEP  
    else:
        # Same operations as used above
        x_circ_maxcheck = x_max - (x_Center + radius)
        x_circ_mincheck = (x_Center - radius) - XCoordinates_Orig
        y_circ_maxcheck = y_max - (y_Center + radius)
        y_circ_mincheck = (y_Center - radius) - YCoordinates_Orig
        circ_array = [1, 2, 3, 4]
        circ_check = [x_circ_maxcheck, x_circ_mincheck, y_circ_maxcheck, y_circ_mincheck]
        circ_check_logical = np.where(circ_check <= 0)
        # Check to see if the circle radius needs to be changed.
        if np.sum(circ_check_logical) > 0:                                          # Revise circle radius
            # Use logical indexing to speed things up.
            circ_logidx = np.where(circ_check_logical == 1)
            circ_check_idx = circ_array[circ_logidx]
            circ_check_radrev = min(circ_check(circ_check_idx))
            radrev = radius + circ_check_radrev
            radius = radrev                                                         # Revise circle radius.                      
        # Calculate the four circle coordinates.
        x1 = x_Center - radius                                                      # Edge 1 of vert. tangent line.
        x2 = x1 + (radius * 2)                                                      # Edge 2 of vert. tangent line.
        y1 = y_Center - radius                                                      # Edge 1 of hor. tangent line.
        y2 = y1 + (radius * 2)                                                      # Edge 2 of hor. tangent line.
        coord = [x1, x2, y1, y2]                                                    # Build a vector.
        #clear(mstring('mask'))                                                     # The mask keeps track of newly taken grid cells.
        # The following operations determines if the test circle overlaps
        # with existing circles, but only to the resolution of the
        # underlying grid. The outcome is a logical matrix. Values of 1 
        # mean that grid point overalps between the test circle and an 
        # existing circle.
        CenterCoord_Array = np.asarray(CenterCoordinates)
        Diameter_Array = np.asarray(Diameter)
        dist_check = ((CenterCoord_Array[:,0] - x_Center) ** 2 + (CenterCoord_Array[:,1] - y_Center) ** 2) ** 0.5
        radius_check = (Diameter_Array / 2) + radius
        radius_check_update = np.reshape(radius_check,(len(radius_check),1))
        radius_overlap = (radius_check_update >= dist_check)
        # Now actually check the overalp condition.
        if np.sum([radius_overlap]) == 0:
            # The new circle does not overlap so proceed.
            newCircle_Found = 1
            debug_value = 2
        elif np.sum([radius_overlap]) == 1:
            # The new circle overlaps with one other circle
            overlap = np.arange(0,len(radius_overlap), dtype=int)
            overlap_update = np.reshape(overlap,(len(overlap),1))
            overlap_logix = (radius_overlap == 1)
            idx_true = overlap_update[overlap_logix]
            radius = dist_check(idx_true,1) - (Diameter(idx_true,1) / 2)
            # New circle radius can be negative when inside another circle.
            if radius > 0:              
                newCircle_Found = 1
                debug_value = 3
            else:
                newCircle_Found = 0
                debug_value = 4
        elif np.sum([radius_overlap]) >= 2:      
            # Overlap with two or more circles.
            # Set newCircle flag to zero so the step is repeated.
            newCircle_Found = 0
            debug_value = 5 
    # Only write data if a circle was placed.
    if step == 1 or newCircle_Found == 1:
        # Write data to variables
        coord = [x1, x2, y1, y2]
        # This function moves circles in a rectangle to available grid points
        mask = ((XCoordinates - x_Center) ** 2 + (YCoordinates - y_Center) ** 2) ** 0.5 <= radius
        # Track which grid points are available for circles and
        # which are not and pass to the next step in the loop
        XYCoordinateslogical[np.where(mask == 1)] = 1
        # Write data to variables
        Coordinates.append(coord)
        CenterCoordinates.append(CenterPosition_data)
        Circle_Area = math.pi * (radius ** 2)         
        Area.append(Circle_Area)
        AreaTotal = np.sum(Area)
        Diameter.append(radius * 2)
        DebugArray.append(debug_value)
        step = step + 1                                                              # Advance counter
        # Calculate the percent completion of the simulation
        percentage_complete = (100.0 * AreaTotal / Total_VArea)
        while len(milestones) > 0 and percentage_complete >= milestones[0]:
            print "{}% complete".format(milestones[0])
            #remove that milestone from the list
            milestones = milestones[1:]  
#
# END OF CODE