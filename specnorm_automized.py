#!/usr/bin/env python
import sys
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cte

from scipy.interpolate import splrep, splev
from astropy.io import fits
from scipy import interpolate

# from http://python4esac.github.io/plotting/specnorm.html

def make_plot(continuum):
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(wave,flux,'k-',label='spectrum')
    axarr[0].plot(wave,continuum,'r-',lw=2,label='Continuum')
    axarr[0].set_title('Continuum fit')
    axarr[0].legend(loc='best')
    axarr[1].plot(wave,flux/continuum,'k-',label='Normalised spectrum')
    axarr[1].axhline(1.01, label='1% margin')
    axarr[1].axhline(0.99)
    axarr[1].set_ylim([-0.1,1.1])
    axarr[1].set_title('Normalized')
    axarr[1].legend(loc='best')
    plt.show()
    
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
        SNR = hdr['SNR50']
        return SNR, bvcor, wave, flux

def onclick(event):
    """
    When no toolbar buttons are activated and the user clicks in the plot
    somewhere, compute 85th percentile value of spectrum in a window 
    (with width of 1/10 bin width at both sides) around the x-coordinate of the 
    clicked point. The y-coordinate of the clicked point is not important. 
    The feel-radius (picker) is set to 5 point.
    :param event: User event
    :return: Continuum point added to graph when clicking left.
    """
    toolbar = plt.get_current_fig_manager().toolbar
    if event.button==1 and toolbar.mode=='':
        window = ((event.xdata-bin_width/10.)<=wave) & (wave<=(event.xdata+bin_width/10.))  # on click, use 1/10 of bin width as window
        y = np.nanpercentile(flux[window], 85)
        plt.plot(event.xdata,y,'ro',ms=5,picker=5,label='cont_pnt')
    plt.draw()

def drawContinuumPoint(pointList, markerColor):
    """
    Draws a continuum point at a specific location.
    :param pointList: List of continuum points containing tuples [wavelength, flux]
    :param markerColor: Color of marker to draw
    :return: Continuum points drawn on the current graph.
    """
    for pt in pointList:
        plt.plot(pt[0], pt[1], 'o', color=markerColor, ms=5, picker=5, label='cont_pnt')
        plt.draw()

def find_nearest(array,value):
    """
    Finds nearest value in a given array.
    :param array: The corresponding Numpy array.
    :param value: The value to search for.
    :return: The index of the nearest corresponding value in the given array.
    """
    idx = (np.abs(array-value)).argmin()
    return idx

def automized_search():
    rIdx = len(wave)-1
    nbBins = 1
    pointList = []

    while rIdx > 0:
        lIdx = find_nearest(wave, wave[rIdx] - bin_width)
        bin_flux = flux[lIdx:rIdx+1]
        bin_wave = wave[lIdx:rIdx+1]
        if (np.isnan(bin_flux).any() or np.isnan(bin_wave).any()):
            print "Nan value encountered in bin %s: interval skipped and continuing." % nbBins
            rIdx = lIdx
            nbBins += 1
            continue
        ctn_flux = np.nanpercentile(bin_flux, 85)
        ctn_wave = bin_wave[find_nearest(bin_flux, ctn_flux)]
        pointList.append([ctn_wave, ctn_flux])
        rIdx = find_nearest(wave, ctn_wave - bin_width/2.) # leave half a bin width angstrom between a continuum point and the start of the next bin
        nbBins += 1

    drawContinuumPoint(pointList, 'red')
    return pointList
    
def onpick(event):
    # when the user clicks right on a continuum point, remove it
    if event.mouseevent.button==3:
        if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
            event.artist.remove()

def ontype(event):
    # When the user hits a:
    # 1. Start on the right side of the spectrum, as HERMES is more accurate there
    # 2. Determine ideal continuum line
    # 3. Draw points
    if event.key=='a':
        automized_search()

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
        continuum = None
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                continuum = artist.get_data()[1]
                break
        if continuum is not None:
            #plt.cla()
            plt.plot(wave,flux/continuum,'k-',label='normalised')
            plt.axhline(1.01)
            plt.axhline(0.99)
            make_plot(continuum)

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
                #np.savetxt(os.path.splitext(filename)[0]+'_norm.dat',data.T)
                np.savetxt(saveName + '_norm.dat', data.T)
                write_fitsfile(data)
                print('Saved to file')
                break
    plt.draw()

def write_fitsfile(data):
    fits_header = fits.getheader(filename)
    hdu = fits.PrimaryHDU()
    hdu.header =fits_header
    hdu.data = data
    hdu.writeto(saveName + '_norm.fits', overwrite=True)


def add_exposures(SNR1, input_wvl_1, input_flux_1, SNR2, input_wvl_2, input_flux_2):
    flux_composed_spectrum = []
    wavelength_composed_spectrum = []
    for w in range (0, min(len(input_wvl_1), len(input_wvl_2))):
        cutoff = 4000
        if log_scale:
            cutoff = np.log(cutoff)
        if w > cutoff:
            if abs(input_wvl_1[w]-input_wvl_2[w]) > 1E-5:
                print 'add_exposures: wavelengths of input fluxes to be combined are not equal'
                sys.exit()
            flux_composed_spectrum.append((input_flux_1[w] + input_flux_2[w])/2)
            wavelength_composed_spectrum.append(input_wvl_1[w])
    flux_composed_spectrum = np.asarray(flux_composed_spectrum)
    #flux_composed_spectrum = (input_flux_1 + input_flux_2)/2
    #wavelength_composed_spectrum = input_wvl_1
    SNR = (SNR1+SNR2)/2 * np.sqrt(2)
    
    return SNR, wavelength_composed_spectrum, flux_composed_spectrum


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
    print "--- Running automated_normalisation.py ---"
    compose = False
    log_scale = False
    bin_width = 0.003       # bin width for the log scale
    # Get the filename of the spectrum and plot it
    if not compose:
        saveName = "00851884"
        filename = "data/00851884_HRF_OBJ_ext_CosmicsRemoved_wavelength_merged_c.fits"
        if 'log' in filename:
            log_scale = True
            saveName += "_log"
            
        SNR, bvcor, wave, flux = read_fits(filename)
        print "SNR of spectrum:", SNR

        if log_scale is False:
            wave, flux = barycentric_correction(bvcor, wave, flux)
            bin_width = 20          # set bin width in angstrom for linear scale
            
        spectrum, = plt.plot(wave,flux,'k-',label='spectrum')
        plt.title(saveName)
        
    else:
        saveName = "SAVENAME"
        filename = "PATH1"
        filename2 = "PATH2"
        if 'log' in filename:
            log_scale = True
            saveName += "_log"
            
        SNR1, bvcor1, wave1, flux1 = read_fits(filename)
        SNR2, bvcor2, wave2, flux2 = read_fits(filename2)

        if log_scale is False:
            bin_width = 20          # set bin width in angstrom for linear scale
            wave1_corr, flux1_corr = barycentric_correction(bvcor1, wave1, flux1)
            wave2_corr, flux2_corr = barycentric_correction(bvcor2, wave2, flux2)
            SNR, wave, flux = add_exposures(SNR1, wave1_corr, flux1_corr, SNR2, wave2_corr, flux2_corr)
        else:
            SNR, wave, flux = add_exposures(SNR1, wave1, flux1, SNR2, wave2, flux2)

        print "SNR of composed spectrum:", SNR

        spectrum, = plt.plot(wave, flux, 'k-', label='spectrum')
        plt.title("Composed spectrum")

    # Connect the different functions to the different events
    plt.gcf().canvas.mpl_connect('key_press_event',ontype)
    plt.gcf().canvas.mpl_connect('button_press_event',onclick)
    plt.gcf().canvas.mpl_connect('pick_event',onpick)
    
plt.show() # show the window
