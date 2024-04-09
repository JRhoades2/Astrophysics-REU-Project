# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:37:40 2019

@author: jorda
"""

from astropy.io import fits
import numpy as np
import math
import decimal

decimal.getcontext().prec = 100

#opening 12CO .fits file
with fits.open('12co_new.fits') as hdu12:
    header12 = hdu12[0].header
    data12 = hdu12[0].data

print(header12)


#opening 13CO .fits file  
with fits.open('13co_allregions.fits') as hdu13:
   header13 = hdu13[0].header
   data13 = hdu13[0].data
   
#test be check if I can read into the .fits files
#print(header12)
#print("\n\n\n")
#print(header13)

#showing shape of numpy arrays
#print(np.shape(data12[0]))
#print(np.shape(data13[0]))

#showing how to look at the length
#print(len(data12[0]))#shows length of velocity channels 
#print(len(data12[0][0]))#show length of y channels
#print(len(data12[0][0][0]))#shows length of x channels

#showing the specific intensity @ vel=38, y=43, x=178
#print(data12[0][38,43,178])


with fits.open('13co_new2.fits') as hdu132:
    header132 = hdu132[0].header
    data132 = hdu132[0].data

#print(header132)
#print(data132[0,78,53,174])
#print(data13[0,78,53,174])

#print(np.shape(data132[0]))
#print(header130['CTYPE1'])

#to do reverse list of a range in decending order:
#x=list(range(20,-1,-1))
#print(x)
#Ico = data12[0]
#finding the exitation Temp
#peak = max(Ico[:,:,:])

#testing going backwards with a range of numbers
#x = []
#for i in range(124,-1,-1):   
#    x.append(i)
#print(x)    
#x = list(range(124,-1,-1))

#time to calculate
#constants
Tex = 30#k #exitation temp for 13CO
T0 = 5.29 #for 13CO
#
#momheader = header12
#
#del momheader['CTYPE3']
#del momheader['CRVAL3']
#del momheader['CDELT3']
#del momheader['CRPIX3']
##del momheader['CUNIT3']
##del momheader['CTYPE4']
#del momheader['CDELT4']
#del momheader['CRPIX4']
##del momheader['CUNIT4']
#del momheader['CRVAL4']
##del momheader['PC1_3']
##del momheader['PC2_3']
##del momheader['PC3_3']
##del momheader['PC4_3']
##del momheader['PC1_4']
##del momheader['PC2_4']
##del momheader['PC3_4']
##del momheader['PC4_4']
##del momheader['PC1_1']
##del momheader['PC2_1']
##del momheader['PC3_1']
##del momheader['PC4_1']
##del momheader['PC1_2']
##del momheader['PC2_2']
##del momheader['PC3_2']
##del momheader['PC4_2']
maheader = momheader
enheader = momheader

#channel width
chanwidth = .340#km/s
#red vel for 13CO
#vel1 = np.arange(6.29106,.170029,-chanwidth)
vel1 = np.arange(6.29106,2.55043,-chanwidth)
#blue vel for 13CO
#vel2 = np.arange(-.17029,-1.1902,-chanwidth)
vel2 = np.arange(-.17029,-1.1902,-chanwidth)
#distance
dist=230.0
#size of a pixel
pixel =  2.777777779420e-03 * 3600
#pixel area
A=(dist*pixel*1.496e13)**2

#placing data that's less than 2sigma to zero for 12CO
LT2sig1=(data12[0]< 0.2).nonzero()
data12[0][LT2sig1]=0.0

#placing data that's less than 2sigma to zero for 13CO
LT2sig2=(data132[0]< 0.2).nonzero()
data132[0][LT2sig2]=0.0

#R_{12/13} ~ 62 Ratio of 12CO/13CO
Ratio = 62
#Technically 12CO made to approximate 13CO

I_new12 = data12[0]/Ratio
I_old13 = data132[0]
##how to copy np.arrays##
momentum = Mint.copy()
momentum[:,:] = decimal.Decimal(0.0)
I13CO = I_old13.copy()
I13CO = decimal.Decimal(0.0)

#Now to create my Frankenstien!
I13CO = I_new12[15:43,0:87,0:330]#28 red vel channels
I13CO = np.concatenate((I13CO,I_old13[44:55,0:87,0:330]))#12 red vel channels
I13CO = np.concatenate((I13CO,I_old13[44:62,0:87,0:330]))#19 red vel channels
I13CO = np.concatenate((I13CO,I_old13[63:66,0:87,0:330]))#4 blue vel channels
I13CO = np.concatenate((I13CO,I_new12[78:98,0:87,0:330]))#21 blue vel channels

LT2SIG = (I13CO < 0.2).nonzero()
I13CO[LT2SIG]=0.0

print(np.shape(I13CO))


#optical depth
#Tau12=-1.0 * np.log(1.0 -((IntensityData)/(T0*(((math.exp(T0/Tex)-1)**(-1.0)) - 0.16))))
Tau13=-1.0 * np.log(1.0 -((I13CO)/(T0*(((math.exp(T0/Tex)-1)**(-1.0)) - 0.16))))
#print(np.shape(Tau13))

#column density
N13=2.5e14 * Tex * Tau13 / (1.0 - math.exp(((-1.0) * T0 / Tex))) * chanwidth

#column density in terms of mole. Hydrogen
NH2=7.0e5 * N13

#conversion
con = A * 3.34e-24/2.0e33

#mass
M = NH2 * con

#summed mass
Mint=np.sum(M,axis=0)

#vel of the mole. cloud
vcloud=4.5#km/s

#creating the momentum np.arrays
momentum = Mint.copy()
momentum[:,:] = decimal.Decimal(0.0)

#creating the energy np.arrays
energy = (momentum.copy())

vel_red = np.arange(16.0925,6.8967,-0.3251)#km/s
#print(len(vel_blue))
vel_blue = np.arange(-4.38886,-10.8909,-0.3251)#km/s
#momentum and energy of entire outflows
for i in range(0,69):
   if(i<29 and i >= 0):
       momentum=momentum+M[i,:,:]*np.sqrt((vel_red[i]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel_red[i]-vcloud)**2*(100000.0)**2)
   elif(i>=29 and i<47):
       momentum=momentum+M[i,:,:]*np.sqrt((vel1[i-29]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel1[i-29]-vcloud)**2*(100000.0)**2)
   elif(i>=47 and i<50):        
       momentum=momentum+M[i,:,:]*np.sqrt((vel2[i-47]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel2[i-47]-vcloud)**2*(100000.0)**2)
   elif(i>=51):
       momentum=momentum+M[i,:,:]*np.sqrt((vel_blue[i-51]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel_blue[i-51]-vcloud)**2*(100000.0)**2)
   
for i in range(0,59):
   if(i<29 and i >= 0):
       momentum=momentum+M[i,:,:]*np.sqrt((vel_red[i]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel_red[i]-vcloud)**2*(100000.0)**2)
   elif(i>=29 and i<41):
       momentum=momentum+M[i,:,:]*np.sqrt((vel1[i-29]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel1[i-29]-vcloud)**2*(100000.0)**2)
   elif(i>=41 and i<44):        
       momentum=momentum+M[i,:,:]*np.sqrt((vel2[i-41]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel2[i-41]-vcloud)**2*(100000.0)**2)
   elif(i>=41):
       momentum=momentum+M[i,:,:]*np.sqrt((vel_blue[i-41]-vcloud)**2)
       energy=energy+(0.5*2.0 *M[i,:,:]*(vel_blue[i-41]-vcloud)**2*(100000.0)**2)
       
   
   

fits.writeto('13coMass3.fits',Mint,maheader)
fits.writeto('13coMomentum3.fits', momentum,momheader)
fits.writeto('13coEnergy3.fits', energy, enheader)


fits.writeto('FrankTest.fits',I13CO,maheader)



