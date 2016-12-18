#!/usr/bin/env python
import numpy as np
a = np.array([12,14,18,23,27,28,34,37,39,40])
print (a-12)/28.0

b = np.array([3,5,10,20,35,40,43,60,25,27])*100.0
print (b-min(b))/(max(b)-min(b))
