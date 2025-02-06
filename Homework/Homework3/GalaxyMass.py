#!/usr/bin/env python
# coding: utf-8

# In[73]:


#import numpy, previously written function Read, and astropy units
import numpy as np
from ReadFile import Read
import astropy.units as u
"""
creates function that takes a file and particle type and retunrs the total mass of of particle trye in file
input : filename speciese the file you wish to use
    partype specifes the type of particel on wants to file 
    'halo', 'disk, 'buldge' ans 1,2,3 are valid enties for partype
output : total mass of specifed partype in 10**12 Msun units
"""
def ComponentMass(filename, partype):
    #use read function to find the data for the given file
    time, num, data = Read(filename)
    #identify chossen particle type by name or number and create mask encluding only this types
    if partype=='disk' or partype==2:
        #create index asking for only disk particles
        index = np.where(data['type'] == 2)
    if partype=='halo' or partype==1:
        #create index asking for only dark/halo matter particles
        index = np.where(data['type'] == 1)
    if partype=='buldge' or partype==3:
        #create index asking for only buldge particles
        index = np.where(data['type'] == 3)
    #apply index to the mass data from file, finding only masses from sepcifed type
    masslist=data['m'][index]
    #sum list of masses and convert to Msun units from 10**10 Msun to sinf total mass
    totalmass=sum(masslist)*10**10
    #converst mass units from Msun to Msun 10**12 and raound to third decimal place
    totalmass=np.round(totalmass/10**12,3)
    #retunr found mass
    print(totalmass)
    return totalmass


# In[75]:


#Find compunents of mass of MW
ComponentMass("MW_000.txt",1)
ComponentMass("MW_000.txt",2)
ComponentMass("MW_000.txt",3)


# In[76]:


#Find compunents of mass of M31
ComponentMass("M31_000.txt",1)
ComponentMass("M31_000.txt",2)
ComponentMass("M31_000.txt",3)


# In[77]:


#Find compunents of mass of M33
ComponentMass("M33_000.txt",1)
ComponentMass("M33_000.txt",2)
ComponentMass("M33_000.txt",3)


# In[ ]:




