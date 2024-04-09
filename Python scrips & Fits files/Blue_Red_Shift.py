#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 10:19:54 2018

@author: Jordan
"""

#libraries
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
from sys import version_info,exit
#from collections import Iterable
assert version_info[0] == 3


#opening .fits files and setting up variables

with fits.open('12co_allregions.fits') as hdu:
    header = hdu[0].header
    data   = hdu[0].data


with fits.open('comb.fits') as hdu:
    header = hdu[0].header
    data1   = hdu[0].data
data1 = np.log(data1[0])
    
#hdul = fits.open('allregions.fits')
#data = hdul[0].data
#header = hdul[0].header

def pix_transform(pixcrd,header):
    """Assumes pix_array is a list or numpy array of pixel values to convert
       Header must be the header object or a dictionary of the header vals"""
    w = wcs.WCS(naxis=2) 
    w.wcs.crpix=[ header['CRPIX{}'.format(x+1)] for x in range(2)]
    w.wcs.cdelt=[ header['CDELT{}'.format(x+1)] for x in range(2)]
    w.wcs.crval=[ header['CRVAL{}'.format(x+1)] for x in range(2)]
    w.wcs.ctype=[ header['CTYPE{}'.format(x+1)] for x in range(2)]
    world = w.wcs_pix2world(pixcrd,1)
    return world

def immoment(data, moment=0,lval=-1,uval=-1):
    if lval == -1: lval = 0
    if uval == -1: uval = data.shape[1]
    assert uval > lval
    if   moment == -1:
        momentdata = np.mean(data[0,lval:uval,:,:],axis=0)
    elif moment == 0 :
        momentdata = np.sum(data[0,lval:uval,:,:],axis=0)
#    elif moment == 1: 
#        momentdata = 
    else:
        print('Other moments not supported')
        exit(0)
    return momentdata        

#plotting in general

red = immoment(data, moment=0,lval=1,uval=42)
blue = immoment(data, moment=0,lval=65,uval=-1)
file = plt.figure()

#plotting protostars
#plt.plot(163.13196,52.315086,'g+',label='L1448C-S', markersize = 7)
#plt.plot(163.8624,52.719186,'y+',label='Per-emb-26',markersize = 7)
#plt.plot(163.1496,52.328021,'m+',label='Per-emb-42',markersize = 7)
#plt.plot(170.6764,56.360426,'g+',label='Per-emb-33-B',markersize = 7)
#plt.plot(170.58033,56.378071,'y+',label='Per-emb-33-A',markersize = 7)
#plt.plot(170.65278,56.348071,'m+',label='Per-emb-33-C',markersize = 7)
#plt.plot(172.52837,57.477571,'g+',label='L1448IRS3',markersize = 7)
#plt.plot(172.47182,57.388412,'m+',label='L1448NW-A',markersize = 7)
#plt.plot(172.4824,57.411319,'y+',label='L1448NW-B',markersize = 7)
#plt.plot(180.60324,46.695253,'g+',label='L1448IRS2E',markersize = 7)
#plt.plot(184.95471,45.529484,'g+',label='Per-emb-22-A',markersize = 7)
#plt.plot(185.00863,45.476675,'g+',label='Per-emb-22-B',markersize = 7)
#plt.plot(185.061,45.456394,'c+',label='L1448IRS2',markersize = 7)
#plt.plot(202.46489,41.098359,'g+',label='L1448IRS1-A',markersize = 7)
#plt.plot(202.42485,40.966362,'c+',label='L1448IRS1-B',markersize = 7)
#plt.legend(loc='upper left')

#making red and blue shifted graphs w/ disabling channels inbetween
plt.contour(blue,colors='blue', levels=list(range(8,76,6)))
plt.contour(red, colors='red',levels=list(range(8,76,6)))
plt.imshow(data1, origin = 'lower')
# plt.tight_layout()

#creating x and y axis into galactic coord.
xtick = range(0,340,50)
ytick = range(0,97, 32)
xarray = np.full([len(xtick),2],0)
yarray = np.full([len(ytick),2],0)
for i,x in enumerate(xtick):
    xarray[i,0] = x
for i,y in enumerate(ytick):
    yarray[i,1] = y  
newx1 = pix_transform(xarray, header)[:,0]
newy1 = pix_transform(yarray, header)[:,1]
newx = newx1
newy = newy1
for i,x in enumerate(newx1):
    newx[i] = round(x,1)
for i,y in enumerate(newy1):
    newy[i] = round(y,2)
    
#plt.title('Mean Flux of L1448')
plt.xlabel('Galactic Longitude ($^\circ$)')
plt.xticks(xtick, newx)
plt.yticks(ytick, newy)
plt.ylabel('Galactic Latitude ($^\circ$)')
plt.draw()
  
file.savefig('TestingOverlaycombMM0.png', dpi=300)

plt.show()




