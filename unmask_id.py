#!/usr/bin/env python

from math import sin
from math import cos
from math import pi
from math import modf
from math import sqrt
import numpy as np
import argparse
from types import SimpleNamespace

from astropy.io import fits as pyfits
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from scipy.signal import find_peaks, find_peaks_cwt
from scipy.signal import argrelextrema

import argparse


def main(mask_image, ID, output_mask, operation='=='):
    hdulist = pyfits.open(mask_image)
    mask = hdulist[0].data
    header = hdulist[0].header

    if operation=='==':
        mask[mask==ID] = 0.
    elif operation=='>=':
        mask[mask>=ID] = 0.
    elif operation=='<=':
        mask[mask<=ID] = 0.        
        
    outHDU = pyfits.PrimaryHDU(mask, header=header)
    outHDU.writeto(output_mask, overwrite=True)
    print('Done!')
    return mask
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("mask", help="Name of the input mask",
                        type=str)
    parser.add_argument("--o", help="Name of the output mask",
                        type=str, default='unmask.fits')
    parser.add_argument("ID", help="ID to unmask",
                        type=int)

    
    args = parser.parse_args()
    mask_image = args.mask
    output_mask = args.o
    ID = args.ID

    
    main(mask_image, ID, output_mask)
            
