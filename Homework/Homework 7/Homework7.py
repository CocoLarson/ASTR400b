#!/usr/bin/env python
# coding: utf-8

# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# 

# In[2]:


# Make edits where instructed - look for "****", which indicates where you need to add code. 


# In[3]:


# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex
get_ipython().run_line_magic('matplotlib', 'inline')

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMassCode import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component

from GalaxyMass import ComponentMass


# # M33AnalyticOrbit

# In[5]:


class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self,filename): # **** add inputs
        """ function that intiates when class is called
        generates positon and velocity vecloty between M33 and M31
        Parameters
        input
        filename: string
            name of file analytical orbit of M33 around M31 will be stored
        """
        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename= filename+".txt"
        
        # **** create an instance of the  CenterOfMass class for M33 
        COM_33=CenterOfMass("M33_000.txt")
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        
        self.r33=COM_33.COM_P(0.1)

        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)

        self.v33=COM_33.COM_V( self.r33[0], self.r33[1], self.r33[2])
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        COM_31=CenterOfMass("M31_000.txt")
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        self.r31=COM_31.COM_P(0.1)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        
        self.v31=COM_31.COM_V( self.r31[0], self.r31[1], self.r31[2])
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        
        self.r0=(self.r33-self.r31).value
        
        self.v0=(self.v33-self.v31).value
       
        ### get the mass of each component in M31 
        #mass from HW 3
        ### disk
        self.rdisk = 5

        self.Mdisk= 0.12*1e12#Msun
        
        ### bulge
        self.rbulge =1 

        self.Mbulge=0.019*1e12
        
        # Halo
        self.rhalo = 62 #kpc

        self.Mhalo=1.921*1e12
     
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        """ 
        Function that calculates the Herniquist accleration assoceated with a given compenent
        found using a=−GM/rmag(ra + rmag)**2*r_vector 
        Parameters
        Input
            M: float
            Mass of the component being measured
            r_a: float
            radius of compnent being measured
            r: array
            position vector of component
        Output
            Hern: array
            Harniquist acceleration assoceated with a specific compnent
        """
        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(sum(r**2))
        
        ### *** Store the Acceleration
        Hern =  -self.G*M/(rmag *(r_a + rmag)**2) * r
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        """Function that calculates the Miyamoto Nagai profile's accleration assoceated with a given compenent
        found using a=-G*M/((R**2+B**2)**(1.5))*zcomp*r R=sqrt(r[0]**2+r[1]**2) 
        B= r_d+sqrt(r[2]**2+( r_d/5)**2)  zcomp=([1,1,B/sqrt(r[2]**2+(r_d/5)**2)]) 
        Parameters
        Input
            M: float
            Mass of the component being measured
            r_d: float
            radius of compnent being measured
            r: array
            position vector of component
        Output
            Miya: array
            Miyamoto Nagai profile's acceleration assoceated with a specific compnent

        """
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        R=np.sqrt(r[0]**2+r[1]**2)

        B= r_d+np.sqrt(r[2]**2+( r_d/5)**2)
        
        zcomp=np.array([1,1,B/np.sqrt(r[2]**2+( r_d/5)**2)]) 
        # where ZSTUFF are the terms associated with the z direction
        
        Miya=-self.G*M/((R**2+B**2)**(1.5))*zcomp*r
       
        return Miya
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self,r): # input should include the position vector, r
        """ Function that uses the previously defined function for acceleration
        to find the acceleration of an objects components and sum them
        Parameters
        Input
        r: array
            posiion vector of object
        Output
        Acc: array
            sumed accelation of object from all masses
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        Halo=self.HernquistAccel(self.Mhalo, self.rhalo, r)
        Bulge=self.HernquistAccel(self.Mbulge, self.rbulge, r)
        Disk=self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 

        Acc=Halo+Bulge+Disk
        return Acc
    
    
    
    def LeapFrog(self,dt,r,v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """
        Functionat that progress the simuation by 1 timestep using a half step methode and 
        the accelatrion suming fucnctio,  rhalf = r+v*dt/2, vnew = v+a(rhalf)*dt, rnew =rhalf+vnew*dt/2
        Parameter
        Input
        Output
        rnew: array
            vector of position after leaf froging
        vnew: array
            vector of velcoity after leaf froging"""
        
        # predict the position at the next half timestep
        rhalf = r+v*dt/2
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v+self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew =rhalf+vnew*dt/2
        
        return rnew, vnew
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ 
        Fucniton that uses LeapFrog  to solve the equations of motion and compute the 
        future orbit of M33 from time t0 to time tmax in dt timesteps and outputs the result as a txt file
        Parameter
        Input
        t0: float
            start time of simulation
        tmax: float
            end time of simulation
        dt: float
            timestep
        Output
        filename.txt: text file
            out puts a text file with the caluculated positon and veloicty after 
            integrating over the full time
        """

        # initialize the time to the input starting time
        t = np.arange(t0,tmax+0.01,dt)
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((len(t),7))
        
        # initialize the first row of the orbit
        #orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        orbit[0,:] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (i<len(t)):  # as long as t has not exceeded the maximal time 

            # **** advance the time by one timestep, dt
           
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t[i]
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)

            pos, vel = self.LeapFrog(dt,orbit[i-1,1:4],orbit[i-1,4:7])
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            orbit[i,1:4]=pos
            orbit[i,4:7]=vel;
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i+=1
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function
        
        


# In[6]:


#create class object for M33AnalyticOrbit
A=M33AnalyticOrbit('M33AnalyticOrbit')


# In[7]:


#use object to create obbit of M33 from 0-10 Gyr using 0.1Gyr intervals
A.OrbitIntegration( 0, 0.1, 10)


# In[8]:


#get data of M33 and M31  orbit from files generted in Hw6 
Sim33=np.genfromtxt('Orbit_M33.txt',skip_header=0)
Sim31=np.genfromtxt('Orbit_M31.txt',skip_header=0)
#get data of M33 and M31 from the anaylisi we just did
Analy33=np.genfromtxt('M33AnalyticOrbit.txt',skip_header=0)


# In[9]:


#find the seperation between M33 and M31 from the simulations 
SimRMag=np.sqrt(np.sum((Sim33[:,1:4]-Sim31[:,1:4])**2,axis=1))
#find the seperation between M33 and M31 from the analsis
AnaRMag=np.sqrt(np.sum(Analy33[:,1:4]**2,axis=1))
#find the seperation between M33 and M31 velocityes from the simulations 
SimVMag=np.sqrt(np.sum((Sim33[:,4:7]-Sim31[:,4:7])**2,axis=1))
#find the seperation between M33 and M31 velocities from the analsis
AnaVMag=np.sqrt(np.sum(Analy33[:,4:7]**2,axis=1))
#get time for simuation and analysis
timeS=Sim33[:,0]
timeA=Analy33[:,0]


# In[25]:


# Plot the Orbit of the galaxies 
#################################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plotting time vs. velocity
plt.plot(timeS, SimVMag, color='blue', 
         linewidth=3, label='Simulation M33-M31 Velocity')
plt.plot(timeA, AnaVMag, color='red', linewidth=3, label='Analytical M33-M31 Velocity')

# Add axis labels
plt.xlabel('Time Gyr', fontsize=22)
plt.ylabel('Seperation Velocity km/s', fontsize=22)


#adjust tick label font size
label_size = 15
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.title( 'Velocity of M33-M31', fontsize=22)

plt.show()


# In[23]:


# Plot the Orbit of the galaxies 
#################################

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plotting time vs. Difference in seperation
plt.plot(timeS, SimRMag, color='blue', 
         linewidth=3, label='Simulation M33-M31 Position')
plt.plot(timeA, AnaRMag, color='red', linewidth=3, label='Analytical M31-M31 Position')

# Add axis labels
plt.xlabel('Time Gyr', fontsize=22)
plt.ylabel('Seperation kpc', fontsize=22)


#adjust tick label font size
label_size = 15
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.title( 'Seperation of M33 and M31 over time', fontsize=22)

plt.show()


# In[ ]:


"""
2.How do the plots compare?
the plots are accurate till about 2 GYR, at which point the data becomes wildly different
3. What missing physics could make the difference?
Dynamic friction, interactions between individual particles (code assumes galaxies act as point masses)
centrifugal force, conservation of angular momentum
4. The MW is missing in these calculations. How might you include its effects?
MW is a massive body with its own gravitational effects on acceleration. To add it in Is effects,
I would add in a third object acceleration using the previous defined function for acceleration
but using the radius and velocity vectors between it and M33 and M31. this would be used to find how
MW affects the change in position/velocity using leapfrog, then integrate.
"""

