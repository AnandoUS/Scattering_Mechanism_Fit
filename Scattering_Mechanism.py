import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
from fdint import *

#This script fits experimental Seebeck and Conductivity data to 
#[Acoustic Phonon Scattering, Point Defect Scattering, Non-polar Optical Scattering] (given by l=0 below), 
#Polar Optical Scattering (given by l=1 below), and 
#ionized impurity scattering (given by l=2 below)
#CAVEAT: These models are only valid for single band transport in single phase semiconductors. So beware of possible impurity phases in your sample. 
#CAVEAT: Also transport theory often treats grain boundaries as a second "phase". So the physics of grain boundary scattering is not explicitly captured.

#####     Universal Constants      ######

kB = 1.38064852 * 10**-23          					#Boltzman Constant SI unit
e = 1.6*10**-19                        					#Charge on electron SI unit
h = 6.626*10**-34                       				#Planks Constant SI unit
pi = np.pi                              				#pi
m_e = 9.10938356 * 10**-31              				#mass of an electron in SI unit


#####     Other Constants         ########

T = 300									#Temperature at which properties were measured
Efs = -8								#Starting Fermi-level
Efe = 90								#Ending Fermi-level
npts = 10000                          					#Number of points between starting and ending Ef


##########################################

df=pd.read_csv('Co_2009.csv')											#Example experimental data
len = len(df.index)

cond = df['Conductivity']
seeb = df['Seebeck']

Efermi=[Ef for Ef in np.linspace(Efs,Efe,npts)]

l=0
S = np.array([10**6*(kB/e)*(((l+2)*fdk(l+1,Ef)/((l+1)*fdk(l+0,Ef)))-Ef) for Ef in np.linspace(Efs,Efe,npts)])	#Seebeck coefficient for scattering mechanism l = 0 
SigmaEo = np.empty([len])

l=1
S1 = np.array([10**6*(kB/e)*(((l+2)*fdk(l+1,Ef)/((l+1)*fdk(l+0,Ef)))-Ef) for Ef in np.linspace(Efs,Efe,npts)])	#Seebeck coefficient for scattering mechanism l = 1
SigmaEo_1 = np.empty([len])

l=2
S2 = np.array([10**6*(kB/e)*(((l+2)*fdk(l+1,Ef)/((l+1)*fdk(l+0,Ef)))-Ef) for Ef in np.linspace(Efs,Efe,npts)])	#Seebeck coefficient for scattering mechanism l = 2
SigmaEo_2 = np.empty([len])

for i in range(len):  
	index = np.argmin(np.abs(S-np.abs(seeb[i])))								
	SigmaEo[i] = (cond[i]/(fdk(0,Efermi[index])))								#SigmaEo for each data point with scattering mechanism l = 0 
	index1 = np.argmin(np.abs(S1-np.abs(seeb[i])))
	SigmaEo_1[i] = (cond[i]/(2*fdk(1,Efermi[index1])))							#SigmaEo for each data point with scattering mechanism l = 1
	index2 = np.argmin(np.abs(S2-np.abs(seeb[i])))
	SigmaEo_2[i] = (cond[i]/(3*fdk(2,Efermi[index2])))							#SigmaEo for each data point with scattering mechanism l = 2

l=0
sig = np.array([np.mean(SigmaEo)*(l+1)*fdk(l+0,Ef) for Ef in np.linspace(Efs,Efe,npts)])

l=1
sig1 = np.array([np.mean(SigmaEo_1)*(l+1)*fdk(l+0,Ef) for Ef in np.linspace(Efs,Efe,npts)])

l=2
sig2 = np.array([np.mean(SigmaEo_2)*(l+1)*fdk(l+0,Ef) for Ef in np.linspace(Efs,Efe,npts)])




######### Plotting the "best fit" for each model along with the experimental data ############

plt.rcParams['figure.figsize'] = 8, 6

plt.scatter(cond, seeb,s=50)						
plt.plot(sig,S, c='steelblue',linewidth=2)
plt.plot(sig1,S1, c='darkgreen',linewidth=1)
plt.plot(sig2,S2, c='firebrick',linewidth=1)

plt.xlabel('Conductivity $\sigma$ (S/m)', fontsize= 20)
plt.ylabel('Thermopower |S| ($\mu$V/K)', fontsize =20)

plt.xticks(fontsize= 18)
plt.yticks(fontsize= 18)

plt.xscale('log')
plt.yscale('log')

plt.xlim(10**3,10**7)
plt.ylim(10,1000)

plt.show()

###############################################################################################