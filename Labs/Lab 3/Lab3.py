
#In Class Lab 3 Template
# G Besla ASTR 400B

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


# The Figure illustrates the color magnitude diagram (CMD) for the Carina Dwarf along with the interpreted 
# star formation history from isochrone fitting to the CMD.
# The image is from Tolstoy+2009 ARA&A 47 review paper about dwarf galaxies
# 
# ![Iso](./Lab3_Isochrones.png)
# 

# # This Lab:
# 
# Modify the template file of your choice to plot isochrones that correspond to the inferred star formation episodes (right panel of Figure 1) to recreate the dominant features of the CMD of Carina (left panel of Figure 1). 



# Some Notes about the Isochrone Data
# DATA From   http://stellar.dartmouth.edu/models/isolf_new.html
# files have been modified from download.  ( M/Mo --> M;   Log L/Lo --> L)
# removed #'s from all lines except column heading
# NOTE SETTINGS USED:  Y = 0.245 default   [Fe/H] = -2.0  alpha/Fe = -0.2
# These could all be changed and it would generate a different isochrone




# Filename for data with Isochrone fit for 1 Gyr
# These files are located in the folder IsochroneData
filename1="./IsochroneData/Isochrone1.txt"


#get file name for 6&7 gyr
filename6="./IsochroneData/Isochrone6.txt"
filename7="./IsochroneData/Isochrone7.txt"

#filw name for 11 &13 gyr
filename11="./IsochroneData/Isochrone11.txt"
filename13="./IsochroneData/Isochrone13.txt"

# READ IN DATA
# "dtype=None" means line is split using white spaces
# "skip_header=8"  skipping the first 8 lines 
# the flag "names=True" creates arrays to store the date
#       with the column headers given in line 8 

# Read in data for an isochrone corresponding to 1 Gyr

data1 = np.genfromtxt(filename1,dtype=None,names=True,skip_header=8)


# get data for 6&7 gyr
data6 = np.genfromtxt(filename6,dtype=None,names=True,skip_header=8)
data7 = np.genfromtxt(filename7,dtype=None,names=True,skip_header=8)

#get data for 11&13 gyr
data11 = np.genfromtxt(filename11,dtype=None,names=True,skip_header=8)
data13 = np.genfromtxt(filename13,dtype=None,names=True,skip_header=8)

# Plot Isochrones 
# For Carina

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot Isochrones

# Isochrone for 1 Gyr
# Plotting Color vs. Difference in Color 


plt.plot(data1['B']-data1['R'], data1['R'], color='blue', linewidth=5, label='1 Gyr')
###EDIT Here, following the same format as the line above 

# plat 6&7 gyr
plt.plot(data6['B']-data6['R'], data6['R'], color='green', linewidth=5, label='6 Gyr')
plt.plot(data7['B']-data7['R'], data7['R'], color='yellow', linewidth=5, label='7 Gyr')

#plot 11 & 13 gyr
plt.plot(data11['B']-data11['R'], data11['R'], color='orange', linewidth=5, label='11 Gyr')
plt.plot(data13['B']-data13['R'], data13['R'], color='red', linewidth=5, label='13 Gyr')

###EDIT Here, following the same format as the line above 

# Add axis labels
plt.xlabel('B-R', fontsize=22)
plt.ylabel('M$_R$', fontsize=22)

#set axis limits
plt.xlim(-0.5,2)
plt.ylim(5,-2.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.6, 0.15, 'CMD for Carina dSph', fontsize=22)

plt.savefig('IsochroneCarina5Iso.png')

plt.show()

# # Q2
# 
# Could there be younger ages than suggested in the Tolstoy plot?
# Try adding younger isochrones to the above plot.
# 
# # Q3
# 
# What do you think might cause the bursts of star formation?
# Galaxy Merger, temp changes-gas cooling, stars die and provide new materials
#Supernova




