#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 13:38:54 2018

@author: smchartrand
"""

#import math
#import random
import time
import numpy as np
# import matplotlib.cm as cm
#
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
Start_time = time.time()

def sample_trig(npoints, r):
    theta = 2 * np.pi * np.random.rand(npoints)
    phi = np.pi - (np.pi * (1 - np.sqrt(np.random.rand(npoints))))
    #phi = np.arccos(2 * np.random.rand(npoints) - 1)
    if np.any(phi > 1.57):
        phi0 = (phi > 1.57)
        phi[phi0] = phi[phi0] - 1.57
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return np.array([x, y, z, theta, phi])

r = 1
#phi2 = np.linspace(0, np.pi, 20)
#theta2 = np.linspace(0, 2 * np.pi, 20)
#xi = r * np.outer(np.sin(theta2),  np.cos(phi2))
#yi = r * np.outer(np.sin(theta2),  np.sin(phi2))
#zi = r * np.outer(np.cos(theta2), np.ones_like(phi2))

npoints = 500
x, y, z, theta, phi = sample_trig(npoints, r)
z0 = (z < 0)

x = x[z0==0]
y = y[z0==0]
z = z[z0==0]


fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'3d', 'aspect':'equal'})
ax.scatter3D(x, y, z, s=10, c='r', zorder=10)
#ax.plot_wireframe(xi, yi, zi, color = 'tab:gray', rstride=1, cstride=1)

print("--- %s seconds ---" % (time.time() - Start_time))
# rotate the axes and update
#for angle in range(0, 360):
#    ax.view_init(60, angle)
#    plt.draw()
#    plt.pause(.1)
