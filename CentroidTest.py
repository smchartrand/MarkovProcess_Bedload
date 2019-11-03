#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:28:49 2018

@author: smchartrand
"""

import numpy as np

data = np.random.randint(0, 10, size=(100000, 2))
def centeroidpython(data):
    x, y = zip(*data)
    l = len(x)
    return sum(x) / l, sum(y) / l
#take the data conversion out to be fair!
data = list(tuple(i) for i in data)