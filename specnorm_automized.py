#!/usr/bin/env python
import sys
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cte

from scipy.interpolate import splrep,splev
from astropy.io import fits


# from http://python4esac.github.io/plotting/specnorm.html

def read_fits(fname):
    """
    Reads fits file
    :param fname: filename
    :return: wave   |   Numpy array containing wavelengths
             flux   |   Numpy array containing fluxes
    """
    with fits.open(fname) as hdu:
        hdr = hdu[0].header
        flux = hdu[0].data
        wave = np.arange(hdr['naxis1']) * hdr['cdelt1'] + hdr['crval1']
        bvcor = hdr['BVCOR']
        return bvcor, wave, flux
    
def onclick(event):
    # when none of the toolbar buttons is activated and the user clicks in the
    # plot somewhere, compute the median value of the spectrum in a .5 angstrom
    # window around the x-coordinate of the clicked point. The y coordinate
    # of the clicked point is not important. Make sure the continuum points
    # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
    toolbar = plt.get_current_fig_manager().toolbar
    if event.button==1 and toolbar.mode=='':
        window = ((event.xdata-.25)<=wave) & (wave<=(event.xdata+.25))
        y = np.median(flux[window])
        plt.plot(event.xdata,y,'rs',ms=10,picker=5,label='cont_pnt')
    plt.draw()

def onpick(event):
    # when the user clicks right on a continuum point, remove it
    if event.mouseevent.button==3:
        if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
            event.artist.remove()

def ontype(event):
    # when the user hits enter:
    # 1. Cycle through the artists in the current axes. If it is a continuum
    #    point, remember its coordinates. If it is the fitted continuum from the
    #    previous step, remove it
    # 2. sort the continuum-point-array according to the x-values
    # 3. fit a spline and evaluate it in the wavelength points
    # 4. plot the continuum
    if event.key=='enter':
        cont_pnt_coord = []
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                cont_pnt_coord.append(artist.get_data())
            elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                artist.remove()
        cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
        sort_array = np.argsort(cont_pnt_coord[:,0])
        x,y = cont_pnt_coord[sort_array].T
        spline = splrep(x,y,k=3)
        continuum = splev(wave,spline)
        plt.plot(wave,continuum,'r-',lw=2,label='continuum')

    # when the user hits 'n' and a spline-continuum is fitted, normalise the
    # spectrum
    elif event.key=='n':
        global flux, wave    
        continuum = None
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                continuum = artist.get_data()[1]
                break
        if continuum is not None:
            plt.cla()
            plt.plot(wave,flux/continuum,'k-',label='normalised')
            
            plt.axhline(1.01)
            plt.axhline(0.99)
            flux = flux/continuum

    # when the user hits 'r': clear the axes and plot the original spectrum
    elif event.key=='r':
        plt.cla()
        plt.plot(wave,flux,'k-')

    # when the user hits 'w': if the normalised spectrum exists, write it to a
    # file.
    elif event.key=='w':
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                data = np.array(artist.get_data())
                np.savetxt(os.path.splitext(filename)[0]+'_norm.dat',data.T)
                print('Saved to file')
                break
    
    # when the user hits 'a': automatically divide spectrum in bins and search for continuum point in each bin
    elif event.key=='a':        
        global cont_wave_list, cont_flux_list
        automized_search(wave, flux, 20)
        for i in range(0, len(cont_flux_list)):
            plt.plot(cont_wave_list[i],cont_flux_list[i],'rs',ms=10,picker=5,label='cont_pnt')    
        cont_wave_list = []
        cont_flux_list = []

    elif event.key=='q':
        global cont_wave_list, cont_flux_list
        automized_search(wave, flux, 100)
        for i in range(0, len(cont_flux_list)):
            plt.plot(cont_wave_list[i],cont_flux_list[i],'rs',ms=10,picker=5,label='cont_pnt')   
        cont_wave_list = []
        cont_flux_list = []            
    plt.draw()

# Find the index of the value, closest to the specified value in a given array.
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def automized_search(wave, flux, bin_width):
    i = len(wave)-1
    while i > 0:
        low_edge = wave[i]-bin_width
        low_index= find_nearest(wave, low_edge)
        
        bin_wave = []
        bin_flux = []
        for j in range(low_index, i):
            bin_flux.append(flux[j])
            bin_wave.append(wave[j])
        bin_flux = np.asarray(bin_flux)
        
        cont_flux = np.percentile(bin_flux, 85)
        index = find_nearest(bin_flux, cont_flux)
        cont_wave = bin_wave[index]
        
        cont_flux_list.append(cont_flux)
        cont_wave_list.append(cont_wave)     
            
        start_new_bin = cont_wave - 10 # leave 10 angstrom between a continuum point and the start of the next bin
        i = find_nearest(wave, start_new_bin)
        
#def ignore_lines(wave, flux)
    #nr_points = 300
    ##for i in range(0, nr_points):
        
    ##for i in range(nr_points, len(wave)-nr_points):
    
    #for i in range(0, len(wave)):
        #if i < nr_points:
            
        #elif i>len(wave)-nr_points:
        
        #else:
def barycentric_correction(bvcor, wvl, flx):
    print "Applying barycentric correction."
    print "Value found in fits file: BVCOR =", bvcor, " km/s"
    deltawvl = wvl * bvcor / cte.c  # The first two lines execute the barycentric correction
    wvl_cor = wvl + deltawvl

    # Create an evenly spaced wavelength vector and evaluate the fluxes on those wavelengths.
    evenlyspacedwvl = np.arange(min(wvl_cor), max(wvl_cor), 0.0156)
    f = interpolate.interp1d(wvl_cor, flx, kind='linear')
    evenlyspacedflux = f(evenlyspacedwvl)
    return evenlyspacedwvl, evenlyspacedflux


if __name__ == "__main__":
    bary_corr = True
    
    global flux, wave, cont_wave_list, cont_flux_list
    cont_flux_list = []
    cont_wave_list = []    
    # Get the filename of the spectrum from the command line, and plot it
    filename = sys.argv[1]
    if filename.endswith('fits') or filename.endswith('fit'):
        bvcor, wave, flux = read_fits(filename)
    else:
        wave, flux = np.loadtxt(filename).T[:2]
    spectrum, = plt.plot(wave,flux,'k-',label='spectrum')
    plt.title(filename)
    
    if bary_corr:
        wave, flux = barycentric_correction(bvcor, wave, flux)
        
    wave_nanremoved = []
    flux_nanremoved = []
    for i in range(0, len(wave)):
        if not np.isnan(flux[i]):
            flux_nanremoved.append(flux[i])
            wave_nanremoved.append(wave[i])
            
    flux = np.asarray(flux_nanremoved)
    wave = np.asarray(wave_nanremoved)

        
    # Connect the different functions to the different events
    plt.gcf().canvas.mpl_connect('key_press_event',ontype)
    plt.gcf().canvas.mpl_connect('button_press_event',onclick)
    plt.gcf().canvas.mpl_connect('pick_event',onpick)
    plt.show() # show the window
