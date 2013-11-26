#!/usr/bin/env python
import sys
import re
import datetime
import h5py
import numpy as np

#By default, open gauges.h5
file = 'gauges.h5'
#But if a filename is passed in the command line, open that
if len(sys.argv)>1:
  file = str(sys.argv[1])
f = h5py.File(file,"r")

#Read the dimensions of the data
#dims[0] is the number of gauges+1 (for timestamps)
#dims[1] is the number of timesteps
dims = f['/dims']

#Read in gauge data
data = f['/gauges']

#Print some data
#to access gauge N (indexed 1...M) at time T:
#data[T*dims[0]+N]
print 'timestamp','gauge 10', 'gauge 23'
for i in range(0,dims[1]):
  print data[i*dims[0]], data[i*dims[0]+11], data[i*dims[0]+24]

#Close data file
f.close()
