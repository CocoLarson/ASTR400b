#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Homework 6 Template
# G. Besla & R. Li


# In[3]:


# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2 
#also removed astronomy units for ease of use
from CenterOfMass2 import CenterOfMass



# In[77]:


def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy: str
            name of galaxy
        start: int
            the inital snap we start usgin
        end: int
            the final snap we will be using
        n: int
            the number of snaps we will skip when doing the calcuatinos
            i.e we will reocrd every nth snap
          
    outputs: 
        Orbit_glaxyname.txt: textfile
            file containing the Com pos and vel as well as the time
    """
    
    # compose the filename for output

    fileout = "Orbit_"+"%s"%(galaxy)+".txt"
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
        #currently set for M33 
        #volDec=2 for M31 and MW
    delta=0.1
    volDec=4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    #if the total number of stans is not perfectly divisable by the number of skipped 
    #snaps end function and print warning
    if ((end-start)%n) ==1:
        print("this is not a valid n, please try a new integer")
        return
    #else create array with the snap ids we will use
    else:
        snapId=np.arange(start, end, n)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM

    orbit=np.zeros([len(snapId), 7])
    
    
    # a for loop 
    for  i in range(len(snapId)):# loop over files
        
        # compose the data filename (be careful about the folder)
        #add a string of the filenumber to the value “000”
        ilbl = '000' + str(snapId[i])
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        #create file name in self
        #this work for me on Jupiter notebook, but may need to bbe modified 
        #depending on where you store your files
        filename="VLowRes/" + "%s"%(galaxy) + "/" + "%s_"%(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        Loop=CenterOfMass(filename)
        #FInd the COM of pos
        COM_P=Loop.COM_P(delta,volDec)
        #Use COM of pos to find COM of vel
        COM_V=Loop.COM_V(COM_P[0],COM_P[1],COM_P[2])
        #find time in gyr and remove units
        time=(Loop.time)/1000/u.Myr
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        orbit[i] = time, *tuple(COM_P),*tuple(COM_V)
    
        # print snap_id to see the progress
        print(snapId[i])
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


# In[7]:


# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
#create files for MW and M31 using velDec=2
OrbitCOM('MW', 0, 800, 5)
OrbitCOM('M31', 0, 800, 5)
# Note: This might take a little while - test your code with a smaller number of snapshots first! 


# In[79]:


#Create files for M33 using velDec=4
OrbitCOM('M33', 0, 800, 5)


# In[81]:


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MW=np.genfromtxt('Orbit_MW.txt',skip_header=1)
M31=np.genfromtxt('Orbit_M31.txt',skip_header=1)
M33=np.genfromtxt('Orbit_M33.txt',skip_header=1)


# In[95]:


# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def mag_dif(a,b):
    #create empty array to store all magnitudes of multiple vectors
    Mag=[0]*len(a)
    #star looping of each set of vectors
    for i in range(len(a)):
        #find diffrence in vectors
        c=a[i]-b[i]
        #find magnitude of diffrence and add to array
        Mag[i]=np.sqrt(sum(c*c))
    #return magnitus
    return Mag


# In[127]:


# Determine the magnitude of the relative position and velocities 
#find time, this will be the same for all calutations
time=MW[:,0]
#isolate postions veclotys from data
MW_pos=MW[:,1:3]
M31_pos=M31[:,1:3]
M33_pos=M33[:,1:3]
#isolate velocity vectors from data
MW_vel=MW[:,4:6]
M31_vel=M31[:,4:6]
M33_vel=M33[:,4:6]
# of MW and M31
#use mag_dif to find postiona and velcity seperation for MW and M31
Mag_MW_M31_pos=mag_dif(MW_pos,M31_pos)
Mag_MW_M31_vel=mag_dif(MW_vel,M31_vel)
# of M33 and M31
#use mag_dif to find postiona and velcity seperation for M33 and M31
Mag_M33_M31_vel=mag_dif(M33_vel,M31_vel)
Mag_M33_M31_pos=mag_dif(M33_pos,M31_pos)


# In[133]:


# Plot the Orbit of the galaxies 
#################################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plotting time vs. Difference in seperation
plt.plot(time, Mag_MW_M31_pos, color='blue', 
         linewidth=3, label='MW-M31')

# Add axis labels
plt.xlabel('Time Gyr', fontsize=22)
plt.ylabel('Seperation kpc', fontsize=22)


#adjust tick label font size
label_size = 15
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.title( 'Sepertaion of MW and M31 over time', fontsize=22)

plt.show()



# In[135]:


# Plot the Orbit of the galaxies 
#################################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plotting time vs. Difference in seperation
plt.plot(time, Mag_M33_M31_pos, color='blue', 
         linewidth=3, label='M33-M31')

# Add axis labels
plt.xlabel('Time Gyr', fontsize=22)
plt.ylabel('Seperation kpc', fontsize=22)


#adjust tick label font size
label_size = 15
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.title( 'Sepertaion of M33 and M31 over time', fontsize=22)

plt.show()



# In[137]:


# Plot the orbital velocities of the galaxies 
#################################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plotting time vs. Difference in seperation
plt.plot(time, Mag_MW_M31_vel, color='blue', 
         linewidth=3, label='MW-M31')

# Add axis labels
plt.xlabel('Time Gyr', fontsize=22)
plt.ylabel('Seperation km/s', fontsize=22)


#adjust tick label font size
label_size = 15
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.title( 'Sepertaion of MW and M31 velocity over time', fontsize=22)

plt.show()



# In[139]:


# Plot the orbital velocities of the galaxies 
#################################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plotting time vs. Difference in seperation
plt.plot(time, Mag_M33_M31_vel, color='blue', 
         linewidth=3, label='M33-M31')

# Add axis labels
plt.xlabel('Time Gyr', fontsize=22)
plt.ylabel('Seperation km/s', fontsize=22)


#adjust tick label font size
label_size = 15
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.title( 'Sepertaion of M33 and M31 velocity over time', fontsize=22)

plt.show()



# In[ ]:


#Q.1
#There will be two close encounters between MW and M31 in the futcer, not encluding the final collsion. 
#Q.2
#When the galaxies are closeest to each other, the smallest distance seperation, they have the largest diffrence in velcoites.
#in short the separation and relative velocity have an inverse realtion ship, when one goes up the other goes down
#Q.3
#MW and M31 colide around 6.25 Gyr. This is then the seperation becomes realivelt flat compared to the rest of time
#after this time M33 continues to orbits the combines MW and M31. It oscilated between 20-110 kpc, 
#gradually getting closer to the remant as time goes on. These oscialtions are more frequent then beofre the merge
#Q.4
#If we take the first apocenter after 6 Gyr (103 kpc speration, at 7.643 Gyr) and the last (62.4 kpc, 10.929 Gyr) we can find the slope
#the the slope is rought the rate of decay. We find that 103.4-62.4/10.929-7.643~-12.5kpc/Gyr or considering this recording is over about
#3 oscilations ~-13.7 kpc/oscilation. If the galaxy is 75 kpc away it will collide in 75 kpc/12.5 kpc/Gyr or 6 Gyr or 5.47 oscilatoins. 
#Assuming constant rate of decay and consistant oscialtions.


# In[247]:


#find positions of apocenters after 6 Gyr
print(Mag_M33_M31_pos[107])
print(Mag_M33_M31_pos[121])
print(Mag_M33_M31_pos[142])
print(Mag_M33_M31_pos[153])
#times of just first and last apocenter
print(time[107])
print(time[153])


# In[243]:


# Plot the Orbit of the galaxies zoom in on area of merge MW and M31
#################################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plotting time vs. Difference in seperation
#plot MW and M31
plt.plot(time[84:], Mag_MW_M31_pos[84:], color='blue', 
         linewidth=3, label='MW-M31')
#plot apocenters
plt.plot(time[107], Mag_M33_M31_pos[107], 'go', 
         linewidth=3, label='MW-M31')
plt.plot(time[121], Mag_M33_M31_pos[121], 'go', 
         linewidth=3, label='MW-M31')

plt.plot(time[142], Mag_M33_M31_pos[142], 'go', 
         linewidth=3, label='MW-M31')
plt.plot(time[153], Mag_M33_M31_pos[153], 'go', 
         linewidth=3, label='MW-M31')
#plot M33 M31
plt.plot(time[84:], Mag_M33_M31_pos[84:], color='red', 
         linewidth=3, label='M33-M31')

# Add axis labels
plt.xlabel('Time Gyr', fontsize=22)
plt.ylabel('Seperation kpc', fontsize=22)


#adjust tick label font size
label_size = 15
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.title( 'Sepertaion of MW,M33 and M31 over time', fontsize=22)

plt.show()



# In[ ]:




