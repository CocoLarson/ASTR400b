#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import moduels
import numpy as np
import astropy.units as u
from astropy.constants import G
#redifine graviation constant before functions so it will always be avalibl
g = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

#importing own functions
from ReadFile import Read
from CenterOfMassCode import CenterOfMass

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')


# In[21]:


#intiale class
class MassProfile:
    def __init__(self, galaxy, Snap):
        """ Function that initates on creatinv a class object
        creates file name as well as sotres mass and position 
        data from file
        Inputs:
            self: self
               class used to store data in so we dont have to keep calling things
            galaxy: string
                name of galaxy name in file
            Snap: float
                snap number of file
            
        """
        #add a string of the filenumber to the value “000”
        ilbl = '000' + str(Snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        #create file name in self
        self.filename="%s_"%(galaxy) + ilbl + '.txt'
        #get data from fiel using read
        self.time, self.total, self.data = Read(self.filename)
        #isolate mass data and save to self
        self.m = self.data['m']
        #isolate position data and save to self
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        #save galaxy name in self
        self.gname = galaxy
        
    def MassEnclosed(self, radii, ptype):
        """ Function that find the mass in Msun enclosed 
        by given radii kpc of a given type of particle
        Inputs:
            self: self
                the class we store data in
            radii: array of astropy quantity
                array of radii we want to use to find mass enclosed in kpc
            ptype: int
                integral representing particle type 
            
        Ouputs:
             MassEnclosed:  array of astropy quantity
                 array with the masses in Msun units
                
        """
        Class=CenterOfMass(self.filename)
        COM_P = Class.COM_P(0.1)
        #create index limiting particel to given ptype
        self.index = np.where(self.data['type'] == ptype)
        #position all particel centerd on COM and apply index
        x=self.x[self.index]-COM_P[0]
        y=self.y[self.index]-COM_P[1]
        z=self.z[self.index]-COM_P[2]
        #find magnitude of radii for all particles
        r_par=np.sqrt(x**2 + y**2 + z**2)
        #set up empty array same lengh as raddi array to store mass enclosed
        MassEnclosed=[0] * len(radii)
        #set up loop gogin over each radii in array
        for i in range(len(radii)):
            #find specific radii in array and units
            rmax=radii[i]*u.kpc
            #create index getting all particel inside specifc radii and apply to masses
            index=np.where(r_par<=rmax)
            m=self.m[index]
            #sum all masses to find total
            MassSum=sum(m)
            #add total mass to array and adjust units
            MassEnclosed[i]=MassSum*1e10
        #add units
        MassEnclosed=MassEnclosed*u.Msun
        #return list of masses
        return MassEnclosed
        
        
    def MassEnclosedTotal(self,radii):
        """ Function that find the mass in Msun enclosed 
        by given radii in kpc of all particles
        Inputs:
            self: self
                the class we store data in
            radii: array of astropy quantity
                array of radii we want to use to find mass enclosed in kpc
            
        Ouputs:
             TotalMass:  array of astropy quantity
                 array with the masses enclosed by the given radii in
                 solar amss units
                
        """       
        #test if galaxy if M33 if so do not find buldge mass
        if self.gname=='M33':
            #use MassEnclosed to find disk and halo mass
            HaloMass = self.MassEnclosed(radii,1)
            DiskMass = self.MassEnclosed(radii,2)
            #sum mass to find total
            TotalMass = HaloMass+DiskMass
        #all other galaxies
        else:
            #use MassEnclosed to find disk, halo buldge mass
            HaloMass= self.MassEnclosed(radii,1)
            DiskMass= self.MassEnclosed(radii,2)
            BuldgeMass= self.MassEnclosed(radii,3)
            #sum mass to find total
            TotalMass=HaloMass+DiskMass+BuldgeMass
        #return total mass
        return TotalMass
        

    #note I toke this function from lab 4 verbaitum
    def hernquist_mass(self,r,h_a,MHalo): 
        """ Function that defines the Hernquist 1990 
        dark matter mass profile 
        Inputs:
            r: astropy quantity
                Galactocentric distance in kpc
            h_a: astropy quantity
                scale radius of the Hernquist profile in kpc
            MHalo: float
                total halo mass in units of 1e12 Msun 
            
        Ouputs:
            mass:  astropy quantity
                total mass within the input radius r in Msun
        """
        #correct units of halo mass, constants
        a=MHalo*1e12*u.Msun
        #find not constant part of mass
        b=r**2/(h_a+r)**2
        #muliply mass parts
        mass = a*b
        
        return mass
        
# cir v = sqrt(MG/r)
    def CircularVelocity(self, radii, partype):
        """ Function that find the circular velocity in km/s 
        of a given radii kpc of a given type of particle
        Inputs:
            self: self
                the class we store data in
            radii: array of astropy quantity
                array of radii we want to use to find velcoity in kpc
            ptype: int
                integer representing particle type 
            
        Ouputs:
             VCirc:  array of astropy quantity
                 array with the circular velocity of the given radii in
                 km/s
                
        """
        #find masses using MassEnclosed
        masses=self.MassEnclosed(radii, partype)
        #find circular velocity of each mass
        VCirc=np.sqrt(masses*g/(radii*u.kpc))
        #return value
        return VCirc

    def CircularVelocityTotal(self,radii):
        """ Function that find the circular velocity in km/s 
        of a given radii kpc of all particle
        Inputs:
            self: self
                the class we store data in
            radii: array of astropy quantity
                array of radii we want to use to find velcoity in kpc
            
        Ouputs:
             VCircTotal:  array of astropy quantity
                 array with the circular velocity of the given radii in
                 km/s
                
        """
        #find masses using MassEnclosedTotal
        masses=self.MassEnclosedTotal(radii)
        #find circular velocity of each mass
        VCircTotal=np.sqrt(masses*g/(radii*u.kpc))
        #return value
        return VCircTotal

    def HernquistVCirc(self,radii,a,HaloMass):
        """ Function that find the circular velocity in km/s 
        of a given radii using the heriquist mass
        Inputs:
            self: self
                the class we store data in
            r: astropy quantity
                Galactocentric distance in kpc
            h_a: astropy quantity
                scale radius of the Hernquist profile in kpc
            MHalo: float
                total halo mass in units of 1e12 Msun 
            
        Ouputs:
             HernVCirc:  array of astropy quantity
                 array with the circular velocity of the given radii in
                 km/s
                
        """
        #create empty array the same leng as the radii array to store values
        #find mass using hernquist_mass
        mass=self.hernquist_mass(r,a,HaloMass)
        #calculate cirular velocity
        HernVCirc=np.sqrt(mass*g/(r*u.kpc))
        #save circular velocity to array
        #return array
        return HernVCirc
            


# In[23]:


#create class for MW
MW = MassProfile('MW',000)
#create array of radii 
r = np.arange(0.1, 30.5, 1.5)
#find specifc partype mass using MassEnclosed
MWHaloMass = MW.MassEnclosed(r,1)
MWDiskMass = MW.MassEnclosed(r,2)
MWBuldgeMass = MW.MassEnclosed(r,3)
#find total mass using MassEnclosedTotal
MWTotalMass = MW.MassEnclosedTotal(r)


# In[ ]:


#find herniquist mass using hernquist_mass
#best fit a=60*kpc
MWHernMass=MW.hernquist_mass(r, 60, 1.975)


# In[ ]:


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot curve for MW

# 
# Plotting radius vs. massenclosed
plt.semilogy(r, MWHaloMass,'o', color='blue', 
         linewidth=5, label='HaloMass')

plt.semilogy(r, MWDiskMass,'o', color='red', 
         linewidth=5, label='DiskMass')

plt.semilogy(r, MWBuldgeMass,'o', color='green', 
         linewidth=5, label='BuldgeMass')

plt.semilogy(r, MWTotalMass,'o', color='black', 
         linewidth=5, label='TotalMass')

plt.semilogy(r, MWHernMass,'o', color='yellow', 
         linewidth=5, label='HernMass')

# Add axis labels
plt.xlabel('Radii kpc', fontsize=22)
plt.ylabel('Log of Mass Enclosed in Solar Mass', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

#add figure text
plt.figtext(0.5, 0.15, 'Mass Profile MW', fontsize=22)
plt.show()


# In[ ]:


#creat object for M31
M31 = MassProfile('M31',000)
#create masses for all types using M31
M31HaloMass = M31.MassEnclosed(r,1)
M31DiskMass = M31.MassEnclosed(r,2)
M31BuldgeMass = M31.MassEnclosed(r,3)
M31TotalMass = M31.MassEnclosedTotal(r)


# In[ ]:


#cerate list of Herniquist mass for M31
#best fit a=60*kpc
M31HernMass=M31.hernquist_mass(r,60,1.921)


# In[ ]:


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot curve for M31

# 
# Plotting radius vs. massenclosed
plt.semilogy(r, M31HaloMass,'o', color='blue', 
         linewidth=5, label='HaloMass')

plt.semilogy(r, M31DiskMass,'o', color='red', 
         linewidth=5, label='DiskMass')

plt.semilogy(r, M31BuldgeMass,'o', color='green', 
         linewidth=5, label='BuldgeMass')

plt.semilogy(r, M31TotalMass,'o', color='black', 
         linewidth=5, label='TotalMass')

plt.semilogy(r, M31HernMass,'o', color='yellow', 
         linewidth=5, label='HernMass')

# Add axis labels
plt.xlabel('Radii kpc', fontsize=22)
plt.ylabel('Log of Mass Enclosed in Solar Mass', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.5, 0.15, 'Mass Profile M31', fontsize=22)
plt.show()


# In[ ]:


#create object for M33
M33 = MassProfile('M33',000)
#create masses for all types except buldge using M33
M33HaloMass = M33.MassEnclosed(r,1)
M33DiskMass = M33.MassEnclosed(r,2)
M33TotalMass = M33.MassEnclosedTotal(r)


# In[ ]:


#cerate list of Herniquist mass for M33
#best fit a=25*kpc
M33HernMass=M33.hernquist_mass(r,25,0.187)


# In[ ]:


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot curve for M33

# 
# Plotting radius vs. massenclosed
plt.semilogy(r, M33HaloMass,'o', color='blue', 
         linewidth=5, label='HaloMass')

plt.semilogy(r, M33DiskMass,'o', color='red', 
         linewidth=5, label='DiskMass')

plt.semilogy(r, M33TotalMass,'o', color='black', 
         linewidth=5, label='TotalMass')

plt.semilogy(r, M33HernMass,'o', color='yellow', 
         linewidth=5, label='HernMass')

# Add axis labels
plt.xlabel('Radii kpc', fontsize=22)
plt.ylabel('Log of Mass Enclosed in Solar Mass', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.5, 0.15, 'Mass Profile M33', fontsize=22)
plt.show()


# In[ ]:


#find variouse veloicites for MW
MWHaloCirV = MW.CircularVelocity(r,1)
MWDiskCirV = MW.CircularVelocity(r,2)
MWBuldgeCirV = MW.CircularVelocity(r,3)
MWTotalCirV = MW.CircularVelocityTotal(r)


# In[ ]:


#find herniquist velcoity using h_a found previousely
MWHernquistCirV=MW.HernquistVCirc(r,60,1.975)


# In[ ]:


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot curve for MW

# 
# Plotting radius vs. massenclosed
plt.plot(r, MWHaloCirV,'o', color='blue', 
         linewidth=5, label='HaloMass')

plt.plot(r, MWDiskCirV,'o', color='red', 
         linewidth=5, label='DiskMass')

plt.plot(r, MWBuldgeCirV,'o', color='green', 
         linewidth=5, label='BuldgeMass')

plt.plot(r, MWTotalCirV,'o', color='black', 
         linewidth=5, label='TotalMass')

plt.plot(r, MWHernquistCirV,'o', color='yellow', 
         linewidth=5, label='HernMassCir')

# Add axis labels
plt.xlabel('Radii kpc', fontsize=22)
plt.ylabel('Circular Velocity km/s', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

#add figure text
plt.figtext(0.5, 0.15, 'Mass Profile MW', fontsize=22)
plt.show()


# In[ ]:


#find the circular velocity for the raddi for all partypres and total mass for M31
M31HaloCirV = M31.CircularVelocity(r,1)
M31DiskCirV = M31.CircularVelocity(r,2)
M31BuldgeCirV = M31.CircularVelocity(r,3)
M31TotalCirV = M31.CircularVelocityTotal(r)


# In[ ]:


#find herniquist velcoity using h_a found previousely
M31HernquistCirV=M31.HernquistVCirc(r,60,1.921)


# In[ ]:


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot curve for M31

# 
# Plotting radius vs. massenclosed
plt.plot(r, M31HaloCirV,'o', color='blue', 
         linewidth=5, label='HaloCirV')

plt.plot(r, M31DiskCirV,'o', color='red', 
         linewidth=5, label='DiskCirV')

plt.plot(r, M31BuldgeCirV,'o', color='green', 
         linewidth=5, label='BuldgeCirV')

plt.plot(r, M31TotalCirV,'o', color='black', 
         linewidth=5, label='TotalCirV')

plt.plot(r, M31HernquistCirV,'o', color='yellow', 
         linewidth=5, label='HernCirV')

# Add axis labels
plt.xlabel('Radii kpc', fontsize=22)
plt.ylabel('Circular Velocity km/s', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

#add figure text
plt.figtext(0.5, 0.15, 'Mass Profile M31', fontsize=22)
plt.show()


# In[ ]:


#find the circular velocity for the raddi for all partypres and total mass for M33
M33HaloCirV = M33.CircularVelocity(r,1)
M33DiskCirV = M33.CircularVelocity(r,2)
M33TotalCirV = M33.CircularVelocityTotal(r)
#find herniquist velcoity using h_a found previousely
M33HernquistCircV=M33.HernquistVCirc(r,25,0.187)


# In[ ]:


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot curve for M33

# 
# Plotting radius vs. circular velocity
plt.plot(r, M33HaloCirV,'o', color='blue', 
         linewidth=5, label='HaloCirV')

plt.plot(r, M33DiskCirV,'o', color='red', 
         linewidth=5, label='DiskCirV')

plt.plot(r, M33TotalCirV,'o', color='black', 
         linewidth=5, label='TotalCirV')

plt.plot(r, M33HernquistCircV,'o', color='yellow', 
         linewidth=5, label='HernCirV')

# Add axis labels
plt.xlabel('Radii kpc', fontsize=22)
plt.ylabel('Circular Velocity km/s', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

#add figure text
plt.figtext(0.5, 0.15, 'Mass Profile M33', fontsize=22)
plt.show()


# In[ ]:





# In[ ]:




