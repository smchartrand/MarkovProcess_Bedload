# -*- coding: utf-8 -*-
# """
# Created on Thu May 24 12:25:15 2018
# @author: schartrand
# """
import math
import random
import time

import matplotlib
import matplotlib.collections
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

# First two lines clear variables from memory each time the script is executed.
# Import the necessary libraries
# from IPython import get_ipython
# get_ipython().magic('reset -sf')

start_time = time.time()
# BEGIN CODE
# This code was developed by Shawn Chartrand and David Furbish, March 2018.
# The code randomly populates a rectangular area of specified dimension with
# circles of random and nonuniform size. The range of circle sizes is
# is specified through a user interface. Placement of circles at step S
# does not depend on prior placements. Placement of circles stops when a
# specified packing fraction (circle volume/total volume) is reached.
# A maximum packing fraction of 70% is generally specified. The code is the
# basis of Markov birth-death simulations of bedload transport in rivers
# under low transport rates.
###############################################################################
###############################################################################
# STEPS ONE, TWO AND THREE DEFINE VARIABLES AND BUILD THE MODEL DOMAIN.
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
# Initialize the step counter.
step = 1
# Circle perimeter coordinate calculations.
theta = np.arange(0, (2 * math.pi), 0.01)
# Define origin for plotting.
XCoordinates_Orig = 0
YCoordinates_Orig = 0
# END SECTION
# STEP TWO: BUILD THE VIRTUAL DOMAIN WITH ALL NECESSARY DETAILS
# KEY INPUT PARAMETER VALUES
# Packing density of the bed; ranges from > 0 to <~0.70
Pack = 0.52
# Length of the domain in the streamwise direction in millimeters
x_max = 2000
# Length of the domain in the cross-stream direction in millimeters
y_max = 100
# Spacing nodes
Grid_Spacing = 2.0
# Size of largest grain relative to the cross-stream length
Size_Frac = 0.06
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
# STEP FOUR: DEFINE FUNCTION USED TO PLACE CIRCLES


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
    return radius, coord
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
    radius, coord = fu_radcheck(x_Center, y_Center, radius)
    # Check to see if radius needs revision for steps greater than the first
    if step != 1:
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
    # Only write data if a circle was placed.
    if step == 1 or newCircle_Found == 1:
        # Write data to variables
        # coord = np.array([x1, x2, y1, y2])
        # This function moves circles in a rectangle to available grid points
        mask = (((XCoordinates - x_Center) ** 2 + (YCoordinates - y_Center)
                 ** 2) ** 0.5 <= radius)
        # Track which grid points are available for circles and
        # which are not and pass to the next step in the loop
        XYCoordinateslogical[np.where(mask == 1)] = 1
        # Write data to variables
        Circle_Area = (np.reshape(np.asfortranarray(math.pi * (radius ** 2),
                                                    dtype=float), (1, +1)))
        Area = np.hstack((Area, Circle_Area))
        AreaTotal = np.sum(Area)
        Coordinates = np.hstack((Coordinates, coord))
        CenterCoordinates = np.hstack((CenterCoordinates, CenterPosition_data))
        Diameter = np.hstack((Diameter, np.reshape(radius * 2, (1, +1))))
        if step == 1:
            Area = np.delete(Area, 0, 1)
            Coordinates = np.delete(Coordinates, 0, 1)
            CenterCoordinates = np.delete(CenterCoordinates, 0, 1)
            Diameter = np.delete(Diameter, 0, 1)
        # Advance counter
        step = step + 1
        # Calculate the percent completion of simulation and print to console
        percentage_complete = (100.0 * (AreaTotal / Total_VArea)/Pack)
        while len(milestones) > 0 and percentage_complete >= milestones[0]:
            print "{}% complete".format(milestones[0])
            # remove that milestone from the list
            milestones = milestones[1:]
# END SECTION
###############################################################################
###############################################################################
# STEP SIX PRINT RESULTS TO SCREEN AND WRITE PERTINENT RESULTS TO FILE
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
# Note that the patches won't be added to the axes, instead a collection will
patches = []
for x1, y1, r in zip(XCenter, YCenter, Radius_Array):
    circle = Circle((x1, y1), r)
    patches.append(circle)
colors = 1000*np.random.rand(len(CenterCoordinates[0]))
p = (PatchCollection(patches, cmap=matplotlib.cm.gist_gray,
                     alpha=0.9, linewidths=(0, )))
p.set_array(colors)
p.set_clim([5, 950])
ax.add_collection(p)
plt.show()
fig.savefig('Circles_test.pdf', format='pdf', dpi=2400)
# Save files for 3D plotting
np.save('XCenter', XCenter)
np.save('YCenter', YCenter)
np.save('Radius_Array', Radius_Array)
# END SECTION
###############################################################################
###############################################################################

# END OF CODE
