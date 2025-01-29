#!/usr/bin/env python
# coding: utf-8

# In[44]:


"""funciton to take a file, find a type of particle, 
find the nth number of this type of paritcle
return position,velocity and mass
input are file
type of particle must be written as 'disk', 'dark', or 'buldge' exactly
number of this type of particle, count start at 1
"""
import numpy as np
import astropy.units as u
from ReadFile import Read
def ParticleInfo(filename,par_type,n):#the function
    #use read function to get data from file
    time, par_num, data=Read(filename)
    #idneifity chossen particle type and create mask encluding only this types
    if par_type=='disk':
        #create index asking for only disk particles
        index = np.where(data['type'] == 2)
    if par_type=='dark':
        #create index asking for only dark/halo matter particles
        index = np.where(data['type'] == 1)
    if par_type=='buldge':
        #create index asking for only buldge particles
        index = np.where(data['type'] == 3)
    #aply index to x, y, z positions, keeping only the desired types of particles 
    x=data['x'][index]
    y=data['y'][index]
    z=data['z'][index]
    #Find position of the nth number of chossen type of particles
    x=x[int(n-1)]*u.kpc
    y=y[int(n-1)]*u.kpc
    z=z[int(n-1)]*u.kpc
    #convert x,y,z position to radial distance
    r=(x**2+y**2+z**2)**(1/2)
    #round radial distancce
    r=np.around(r,3)
    #aply index to vx, vy, vz velocity, keeping only the desired types of particles 
    vx=data['vx'][index]*u.km/u.s
    vy=data['vy'][index]*u.km/u.s
    vz=data['vz'][index]*u.km/u.s
    #Find velocity of the nth number of chossen type of particles
    vx=vx[int(n-1)]
    vy=vy[int(n-1)]
    vz=vz[int(n-1)]
    #find the magnitude of the velcoity
    vr=(vx**2+vy**2+vz**2)**(1/2)
    #round velcoity to 3rd decimal
    vr=np.around(vr,3)
    #find the mass of the chossen type then find nth number particle
    m=data['m'][index][int(n-1)]
    #convert 10**10 solar mass to 1 solar amss by dividing by 10**10
    mass=m/10**10*u.M_sun
    #print calculated values of distance, velocity and mass
    print('Magnitude of the distance '+str(r))
    print('Magnitude of the velocity '+str(vr))
    print('Mass '+str(mass))
    #return values of distance, velocity and mass
    return r, vr, mass
    

