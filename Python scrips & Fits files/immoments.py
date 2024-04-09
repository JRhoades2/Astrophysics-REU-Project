#Code made by Nick Reynolds
#used by Jordan Rhoades

#importing all of the libraries needed
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from sys import version_info,exit
from collections import Iterable
assert version_info[0] == 3
from astropy import wcs


hdul = open.fits('12co_allregions.fits')
header = hdul[0].header
data = hdul[0].data
'''
#
def typecheck(obj): 
    return not isinstance(obj, str) and isinstance(obj, Iterable)

#opening .fits file and making the data a glorified array
def reading(fname):
    with fits.open(fname) as hdu:
        header = hdu[0].header
        data   = hdu[0].data
    return data,header

#editing the fits file but not used in this case
def saving(fname,data,header):
    fits.writeto(fname,data,header)
    pass

#having it where dependig on the data, we could be getting a -1 moment or a 0 moment
def immoment(data,header=None, moment=-1,axis='channel',lval=-1,uval=-1):
    if lval == -1: lval = 0
    if uval == -1: uval = data.shape[1]
    assert uval > lval
    if   moment == -1:
        momentdata = np.mean(data[0,lval:uval,:,:],axis=0)
    elif moment == 0 :
        momentdata = np.sum(data[0,lval:uval,:,:],axis=0)
    else:
        print('Other moments not supported')
        exit(0)
    return momentdata
#round to specific decimal
    
'''
'''def changecoordx(pixel):
    degreex = round(((-0.002777)*((pixel) - 142)) + 158.127 , 1)
    return degreex
#same as x, fix rounding to like ~2 - 3 decimal places
def changecoordy(pixel):
    degreey = ((0.002777)*((pixel) - 41.999)) - 21.45
    return degreey
'''
# assumes header
'''
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

#essentially what it says, you are plotting the data
def plotter(data,header,moment=0,cmap='viridis',save=False):
    title,xlabel,ylabel,zlabel=header['OBJECT'],header['CTYPE1'],header['CTYPE2'],header['BUNIT']
    #locs, labels = plt.xticks()
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
   
   
    #xax = ['158$^\circ$.5','','158$^\circ$.3','','158^\$circ$.1','','157$^\circ$.9','','157$^\circ$.7','',]
    
    plt.imshow(data)
    plt.title('Mean Flux of L1448')
    plt.xlabel('Galactic Longitude ($^\circ$)')
    plt.xticks(xtick, newx, rotation = 45)
    plt.yticks(ytick, newy, rotation = 45)
    plt.ylabel('Galactic Latitude ($^\circ$)')
    plt.gca().invert_yaxis()
    cbar = plt.colorbar()
    plt.tight_layout()
    plt.draw()
    if moment == -1: cbar.ax.set_ylabel('K')
    else: cbar.ax.set_ylabel(r'K*(km*s$^{-1}$)')
    if save: plt.savefig(save,dpi=400)
    plt.show()
    
#

#main function that calls to the other used functions
def main(args):
    data,header = reading(args.input)
    #data = np.log10(data1)
    if (type(args.moments) == int) or (type(args.moments) == float):
        moment=immoment(data,header,moment=args.moments,lval=args.lval,uval=args.uval)
        plotter(moment,header,moment=moment,save=args.save)
    else:
        args.moments = [int(x) for x in args.moments.strip('[').strip(']').split(',')]
        for x in args.moments:
            moment=immoment(data,header,moment=x,lval=args.lval,uval=args.uval)
            plotter(moment,header,moment=x,save=args.save)

#looks like it basically a glorified helper in the sense that you can use all functions
if __name__ == "__main__":
    from argparse import ArgumentParser as ap
    parser = ap()
    parser.add_argument('-i', '--input',help='Input fits file name',default='image.fits',required=True)
    parser.add_argument('-s', '--save',help='Output pdf file name')
    parser.add_argument('-o', '--output',help='Output fits file name',default='cores.fits')
    parser.add_argument('-m', '--moments',help='Integer or list of moment maps to compute',default=-1)
    parser.add_argument('--lval',help='Lower channel value',default=-1)
    parser.add_argument('--uval',help='Upper channel value',default=-1)
    args = parser.parse_args()

    
    main(args)
'''
# end of code
