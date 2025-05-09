#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import numpy
import numpy as np


# In[3]:


def ConcatenateHalo(name1, name2, filename):
    """ Combine data of two diffrent txt files and 
    return as a new text file with with a 1 or 8 in the last 
    column indicating the orginal file.
    Parameters
    ----------
    name1 : str
        first file name 
    name2 : str
        second file name
    filename: str 
        name of the file you would like data returned as
          
    """
    #Get data for snap for given files
    file1=np.genfromtxt(name1,skip_header=3)
    file2=np.genfromtxt(name2,skip_header=3)

    #create index the will remove all non halo particles
    index = np.where(file1[:,0] == 1)
    #Create array full of ones with one more columb then the first set of data
    a = np.full((len(file1),len(file1[0])+1),1,dtype='float')
    #add the data to new array, columb 8 keeps its 1
    a[:,:-1] = file1
    #remove none halo particles
    a=a[index]
    #create index the will remove all non halo particles
    index = np.where(file2[:,0] == 1)
    #Create array full of twos with one more columb then the second set of data
    b = np.full((len(file2),len(file2[0])+1),2,dtype='float')
    #add new data to new array, columb 8 keeps its 2
    b[:,:-1] = file2
    #remove none halo particles
    b=b[index]
    #combine  arrays
    halo=np.concatenate((a,b))

    #create txt file with combined halo paricles so we can use center of mass class
    np.savetxt(filename+".txt", halo,  comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format( 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz','gal'))


# In[ ]:




