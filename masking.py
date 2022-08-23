import numpy as np
import numpy.ma as ma
import os
from astropy.io import fits
from photutils.segmentation import make_source_mask
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from photutils.background import Background2D
from photutils.segmentation import detect_threshold
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils.segmentation import detect_sources
from photutils.segmentation import deblend_sources

import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import LogStretch
#--------
from glob import glob

from unmask_id import *

def import_data():

    input_file = np.genfromtxt('/home/maria/Documents/projects/Stripe82/sample_to_check.txt')
    RA = input_file[:,0]
    DEC = input_file[:,1]
    return(RA, DEC)

def get_image_data(image_file):
    data = fits.getdata(image_file, ext=0)
    return data

def initial_mask(data):
    mask = make_source_mask(data, nsigma=2, npixels=5, dilate_size=30)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    print('mean, median, std:', (mean, median, std))
    return(mask)

def background_calculations(data):
    coverage_mask = (data == 0)
    sigma_clip = SigmaClip(sigma=3.)
    mask = initial_mask(data)
    bkg = Background2D(data, (25, 25), filter_size=(3, 3),exclude_percentile=10.0,  coverage_mask = coverage_mask, sigma_clip=sigma_clip, mask = mask,  fill_value = 0.0)
    return(bkg)

def segmentation_map(data):
    #threshold = bkg.background + (2.0 * bkg.background_rms)
    threshold = detect_threshold(data, nsigma=2.)
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold, npixels=5,  filter_kernel=kernel)
    segm_deblend = deblend_sources(data, segm, npixels=5, filter_kernel=kernel, nlevels=32, contrast=0.001)
    return(segm_deblend)

def get_mask_index(RA, DEC):
    image_file = '/home/maria/Documents/projects/Stripe82/data_additional/EG_%.5f_%.5f/rdeep/segm.fits' % (RA,DEC)
    data = fits.getdata(image_file, ext=0)
    mask_index = data[int(data.shape[0]/2), int(data.shape[1]/2)]
    return mask_index

def put_a_mask_on_the_image(RA, DEC):
    input_dir = '/home/maria/Documents/projects/Stripe82/data_additional/EG_%.5f_%.5f/rdeep/' % (RA,DEC)
    image_file = input_dir +'EG_%.5f_%.5f_rdeep.fits' % (RA, DEC)
    data = get_image_data(image_file)
    mask_file = input_dir +'mask.fits'
    mask = fits.getdata(mask_file, ext = 0)
    data_with_mask = ma.array(data, mask=mask, fill_value = 0)
    #data_with_mask.filled()
    #norm = ImageNormalize(stretch=LogStretch())
    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
    # ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
    # ax1.set_title('Data')
    # ax2.imshow(data_with_mask.filled(), origin='lower', cmap='Greys_r', interpolation='nearest', norm=norm)
    # ax2.set_title('Image with a mask')
    #plt.show()
    fits.writeto(input_dir + 'data_with_mask.fits', data_with_mask.filled(), overwrite=True)


  

def my_main():
    RA, DEC = import_data()
    list_of_indexes = [3, 6, 23, 29, 30, 32, 33, 35, 36, 40, 43, 47,  65, 70, 102, 117]
    #list_of_indexes = np.genfromtxt('/home/maria/Documents/projects/Stripe82/golden_sample.txt')
    for i in list_of_indexes:
        i = int(i)
    
        input_dir = '/home/maria/Documents/projects/Stripe82/data_additional/EG_%.5f_%.5f/rdeep/' % (RA[i],DEC[i])
        image_file = input_dir+ "EG_%.5f_%.5f_rdeep.fits" % (RA[i], DEC[i])
        
        data = get_image_data(image_file)
        data[np.isnan(data)] = 0
        #bkg = background_calculations(data)

        norm = ImageNormalize(stretch=LogStretch())
        segm = segmentation_map(data)

        fits.writeto(input_dir + 'segm.fits', segm.data, overwrite=True)

        mask_index = get_mask_index(RA[i], DEC[i])
        mask = main(input_dir+'segm.fits', mask_index, input_dir + 'mask.fits')
        put_a_mask_on_the_image(RA[i], DEC[i])

        print(i, RA[i], DEC[i])

my_main()
