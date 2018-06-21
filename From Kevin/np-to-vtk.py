import pyevtk
from pyevtk import hl 
import numpy as np


R_np = np.load('Radius_Array.npy')
X_np = np.load('XCenter.npy')
Y_np = np.load('YCenter.npy')
Z_np = np.zeros(len(R_np)) # z coordinate needs specification

hl.pointsToVTK('XYZR',X_np,Y_np,Z_np,{'radius':R_np})
