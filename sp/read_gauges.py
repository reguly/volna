#!/usr/bin/env python
import sys
import re
import datetime
import h5py
import numpy as np

f = h5py.File("gauges.h5","r")
dims = f['/dims']
data = f['/gauges']
print 'timestamp','gauge 10', 'gauge 23'
for i in range(0,dims[1]):
  print data[i*dims[0]], data[i*dims[0]+10], data[i*dims[0]+23]
f.close()
