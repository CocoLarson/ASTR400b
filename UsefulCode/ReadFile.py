#!/usr/bin/env python
# coding: utf-8

# In[11]:


# import needed modules
import numpy as np
import astropy.units as u
#define read function that takes file and returns data
def Read(filename):
    #open file
    file = open(filename,'r' )
    #read first line and store time in units of Myr
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    #read the second line and store number of particles 
    line2 = file.readline()
    labe2, value = line2.split()
    num = float(value)
    #close file
    file.close()
    #store rest of data as array, keeping the column header informatio
    data = np.genfromtxt(filename,dtype=None,names=True, skip_header=3)
    #Return time, number of particles and data
    return time, num, data

