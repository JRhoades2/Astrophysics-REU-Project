# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 16:30:40 2018

@author: jorda
"""

from astropy.io import fits
import numpy as np
import math
import decimal

decimal.getcontext().prec = 100

with fits.open('12co_new.fits') as hdu:
    header = hdu[0].header
    data   = hdu[0].data
    
momheader = header

del momheader['CTYPE3']
del momheader['CRVAL3']
del momheader['CDELT3']
del momheader['CRPIX3']
#del momheader['CUNIT3']
#del momheader['CTYPE4']
del momheader['CDELT4']
del momheader['CRPIX4']
#del momheader['CUNIT4']
del momheader['CRVAL4']
#del momheader['PC1_3']
#del momheader['PC2_3']
#del momheader['PC3_3']
#del momheader['PC4_3']
#del momheader['PC1_4']
#del momheader['PC2_4']
#del momheader['PC3_4']
#del momheader['PC4_4']
#del momheader['PC1_1']
#del momheader['PC2_1']
#del momheader['PC3_1']
#del momheader['PC4_1']
#del momheader['PC1_2']
#del momheader['PC2_2']
#del momheader['PC3_2']
#del momheader['PC4_2']
maheader = momheader
enheader = momheader
#momheader['BUNIT'] = 'Msun * km * s^-1'
#enheader['BUNIT'] = 'E43 erg'
#print(data[0][38,43,178])
#constants
T0=6.626e-27*115.271203e9 * 2/1.38e-16#K
Tex=30#K
chanwidth=0.3251#km/s
vel=np.arange(20.644,-20.644,-0.3251,)#km/s
vel_red = np.arange(16.0925,6.8967,-0.3251)#km/s
#print(len(vel_blue))
vel_blue = np.arange(-4.38886,-10.8909,-chanwidth)#km/s
#print(len(vel_red))
#print(len(vel))
#print(header)
dist=230.0
pixel =  2.777777779420e-03 * 3600

#placing data that's less than 2sigma to zero
LT2sig=(data[0]< 0.2).nonzero()

data[0][LT2sig]=0.0

#intensity of 12CO
Ico=data[0]/ 70
#optical depth
Tau12=-1.0 * np.log(1.0 -((Ico)/(T0*(((math.exp(T0/Tex)-1)**(-1.0)) - 0.16))))
Tau12_blue=-1.0 * np.log(1.0 -((Ico[78:98,:,:])/(T0*(((math.exp(T0/Tex)-1)**(-1.0)) - 0.16))))
Tau12_red=-1.0 * np.log(1.0 -((Ico[15:43,:,:])/(T0*(((math.exp(T0/Tex)-1)**(-1.0)) - 0.16))))
#column density
N12=2.5e14 * Tex * Tau12 / (1.0 - math.exp(((-1.0) * T0 / Tex))) * chanwidth
N12_blue=2.5e14 * Tex * Tau12_blue / (1.0 - math.exp(((-1.0) * T0 / Tex))) * chanwidth
N12_red=2.5e14 * Tex * Tau12_red / (1.0 - math.exp(((-1.0) * T0 / Tex))) * chanwidth
#print(N12[38,43,178])
#pixel area
A=(dist*pixel*1.496e13)**2
#print(A)
NH2=7.0e5 * N12 
NH2_blue=7.0e5 * N12_blue
NH2_red=7.0e5 * N12_red 
con = A * 3.34e-24/2.0e33
#mass
M = NH2 * con
M_blue=NH2_blue * con
M_red=NH2_red * con
#summed mass
Mint=np.sum(M,axis=0)
Mint_blue=np.sum(M_blue,axis=0)
Mint_red=np.sum(M_red,axis=0)
#print(Mint.shape)
#column density summed over all velocity channels
#print(Mint[43,178])
vcloud=4.5#km/s
momentum = Mint.copy()
momentum[:,:] = decimal.Decimal(0.0)
momentum_blue = Mint.copy()
momentum_blue[:,:] = decimal.Decimal(0.0)
momentum_red = Mint.copy()
momentum_red[:,:] = decimal.Decimal(0.0)
#print(momentum.shape)
energy = (momentum.copy())
energy_blue = (momentum_blue.copy())
energy_red = (momentum_red.copy())
#print(Mint[43,178])

#momentum and energy of entire outflows
for i in range(0,128):
    momentum=momentum+M[i,:,:]*np.sqrt((vel[i]-vcloud)**2)
    energy=energy+(0.5*2.0 *M[i,:,:]*(vel[i]-vcloud)**2*(100000.0)**2)
#momentum and energy of blue-doppler shifted outflows
for i in range(0,21):
    momentum_blue=momentum_blue+M[78+i,:,:]*np.sqrt((vel_blue[i]-vcloud)**2)
    energy_blue=energy_blue+(0.5*2.0 *M[78+i,:,:]*(vel_blue[i]-vcloud)**2*(100000.0)**2)#multiply by e33
#momentum and energy of red-doppler shifted outflows
for i in range(0,29):
    momentum_red=momentum_red+M[15+i,:,:]*np.sqrt((vel_red[i]-vcloud)**2)
    energy_red=energy_red+(0.5*2.0 *M[i,:,:]*(vel_red[i]-vcloud)**2*(100000.0)**2)#multiply by e33
    
fits.writeto('12coMass.fits',Mint,maheader)
fits.writeto('12coMomentum.fits', momentum,momheader)
fits.writeto('12coEnergy.fits', energy, enheader)

fits.writeto('12coMassBlue.fits', Mint_blue,maheader)
fits.writeto('12coMassRed.fits',Mint_red,maheader)
fits.writeto('12COMBLUE.fits',momentum_blue,momheader)
fits.writeto('12COMRED.fits',momentum_red,momheader)
fits.writeto('12COEBLUE.fits',energy_blue,enheader)
fits.writeto('12COERED.fits',energy_red,enheader)
    




