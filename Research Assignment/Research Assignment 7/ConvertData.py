#!/usr/bin/env python
# coding: utf-8

# In[15]:


#import numpy
import numpy as np


# In[17]:


def ConvertData(par,x,y):
    """ Conver array of x,y,z postion data to counts 
    given the choosen plane x=0, y=1, z=2,
     returns counts, xedges, yedges, x_binsize, y_binsize
    ----------
    Input
    par : array
        array of particles x,y,z position
    x : int
        first axis of choosen dimenshion
    y : int
        second axis of choosen dimenshion
    OUtput
    counts : numpy.ndarray
        array with the number of particles in each bin 
    xedges : numpy.ndarray
        array of x edges
    yedges : numpy.ndarray
        array of y edges
    x_binsize : numpy.float64
        bin size in the x direction
    y_binsize : numpy.float64
        bin size on the y direction """

    #all particles xy plane
    counts, xedges, yedges = np.histogram2d(par[:,x],par[:,y], bins = 100, range = ((-2000, 2000), (-2000, 2000)))
    #our bin size is:
    x_binsize = xedges[1] - xedges[0]
    y_binsize = yedges[1] - yedges[0]

    return counts, xedges, yedges, x_binsize, y_binsize 


# In[ ]:




