#!/usr/bin/env python
# coding: utf-8

# # In Class Lab 1
# 
# ### Due by midnight, thursday in your github repository 'Labs/Lab1' folder
# 

# In[1]:


# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants


# ### Astropy Units:  https://docs.astropy.org/en/stable/units/index.html

# 
# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 

# In[5]:


""" function that will use 4.74*u*R=VLSR-v_s
    to find the VLSR 
    input: the solar radius as R in KPC
    output: the VLSR in units km/s
"""
def VLSR(R):#function
    #mu is the proper motion of sgr a* mas/yr
    mu=6.379
    #v_s is peculiar motin of sun km/s
    v_s=12.24
    #calculate VLSR using equaotin above
    v_lsr=4.74*mu*(R/u.kpc)*(u.km/u.s)-v_s*(u.km/u.s)
    #return vlsr unit km/s
    return v_lsr


# In[ ]:





# In[8]:


#run VLSR for water maser distance and print
WMd_V=VLSR(8.43*u.kpc)
print('Water Maser distance = '+str(WMd_V))
#run VLSR for GRAVITY Collaboration distance and print
Grav_V=VLSR(8.178*u.kpc)
print('Gravity Collaboration = '+str(Grav_V))
#run VLSR for Sparke & Gallagher distance and print
sparke_V=VLSR(7.9*u.kpc)
print('Sparke & Gallager = '+str(sparke_V))


# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr

# In[56]:


#define GRAVITY Collaboration distance
#orbial period=2piR/V
def TOrbSun(R,v):
    """funcitons that computes the orbial period T= 2 pi R/v
    input R : the radius to the galactric center kpc
    v velocity of the sun in km/s
    out put: period"""
#calcute distance traveled by sun in one rotation/ the circumfrance in kpc
    cir=2*np.pi*R
    #convert velocity from km//s to kpc/gyr
    VkpcGyr= v.to(u.kpc/u.Gyr)
    #use VLSR funciton to find VLSR to sun at the given radius
    Grav_V=VLSR(R)
    #convert VLSR to proper units
    Grav_V.to(u.kpc/u.Gyr)
    #add peculiar and VLSR together to get total velocity of sun
    Vsun=Grav_V+VkpcGyr
#preiod is distance traveled divided by speed
    T=cir/Vsun
    return T
#Show result


# In[58]:


#Run TOrbSUn function with r=8.178*u.kpc from the GRAVITY Collaboration and v=12.24*u.kpc/u.Gyr 
period=TOrbSun(8.178*u.kpc, 12.24*u.km/u.s)
print(str(period))


# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

# In[62]:


#Find period using TOrbSUn functions
period=TOrbSun(8.178*u.kpc, 12.24*u.km/u.s)
#define the total time of the universe
Time=13.8*u.Gyr
#divide total time by the time it takes to complete one rotation to get total number of rotations
num_rot=Time/period
#show results
print(str(num_rot)+' Rotations')


# In[64]:


NumRot(8.178*u.kpc, 12.24*u.km/u.s,13.8*u.Gyr)


# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 

# In[68]:


#find units of astopy G
print(const.G)


# In[74]:


#convert m3/ kg s2 to kpc/gyr/msun
Grav=const.G.to(u.kpc**3/u.Gyr**2/u.Msun)
print(Grav)


# In[84]:


#p=VLSR**2/(4*pi*G*R**2)
#mass = integrate p dv
#     = p 4*pi*r**2*dr
#intergrate VLSR**2/(4*pi*G*R**2)* p 4*pi*r**2*dr
#Mass = int VLSR**2/G dr
#VLSR**2/G*r
def massIso(R, VLSR):
    """function that cp the dark matter mass eloced a given distance R, assuming isothermal sphere
    M=VLSR**2/G*r
    input: r-distance from GC (in kpc
        VLSR: velocity of local standard of rest km/s
    output : M: the dark matter mass enclosed by a sphere of r in Msun"""
    #convert velocty from km/s to kpc/gyr
    VLSRkpcGyr=VLSR.to(u.kpc/u.Gyr)
    #isohermal sphere mass profile
    M= VLSRkpcGyr**2/Grav*R
    return M


# In[92]:


#compute mass from Grav Collad radius and knonw VLSR (235.033 km/s from part 1)
M=massIso(8.178*u.kpc, Grav_V)
print(f"{M:.2e}")


# In[94]:


#find mass at 260 kpc assuming VLSr does not change
M=massIso(260*u.kpc, Grav_V)
print(f"{M:.2e}")


# In[96]:


#find mass at 260 kpc assuming VLSr is 220/km/s
M=massIso(260*u.kpc, 220*u.km/u.s)
print(f"{M:.2e}")


# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

# In[102]:


#poteinal for a Hern Sphere
#phi= -G*M/(r+a)
#vesc**2=2*G*M/(r+a)
#rearange M=vesc**2/2/G*(r+a)
def HernMass(Vesc, r, a=30*u.kpc):
    """ function calulates mass using Hernquist sphere poteinal
    input: Vesc:the escape velocity of an object in km/s
        R: the distaince of the object from the galactic center in kpc
        a: Hernquist scale radius in kpc, set to 30 kpc
    output: M: Mass with in units Msun"""
    #convert vesc to astopy units
    Vesckpcgyr=Vesc.to(u.kpc/u.Gyr)
    #calcuate mass
    M= Vesckpcgyr**2/2/Grav*(r+a)
    #retunr mass
    return M


# In[106]:


#run Hern Mass with raidus of Leo 1 r=260pkc and velcoity is 196km/s
HM=HernMass(196*u.km/u.s,260*u.kpc)
print(f"{HM:.2e}")


# In[ ]:




