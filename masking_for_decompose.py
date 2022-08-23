import numpy as np
from astropy.io import fits
import subprocess
import numpy.ma as ma  
import os  



def read_from_file():

    input_file = np.genfromtxt('/home/maria/Documents/projects/Stripe82/super_final_sample1.txt')
    #input_file = np.genfromtxt('/home/maria/Documents/projects/Stripe82/params.cat')

    return (input_file)

def get_image_data(image_file):
    data = fits.getdata(image_file, ext=0)
    return data


def make_mask():
    """
    Use sextractor to make a mask
    """
    # Check if we have special sextractor parameters

    deblend_nthresh = 5
    deblend_mincont = 0.0002
    magzpt = 28.5
    target_x = image_data.shape[1] // 2
    target_y = image_data.shape[0] // 2

    callString = f"sex {image_file} -c ./libs/default.sex -MAG_ZEROPOINT {magzpt} "
    callString += f"-DEBLEND_NTHRESH {deblend_nthresh} -DEBLEND_MINCONT {deblend_mincont} "
    callString += f"-CHECKIMAGE_TYPE  SEGMENTATION -CHECKIMAGE_NAME  {input_dir}/segm.fits -VERBOSE_TYPE QUIET"
    subprocess.call(callString, shell=True)
    mask = fits.getdata(f"{input_dir}segm.fits")
    # Remove the target object from the mask
    i = mask[target_y, target_x]
    inds = np.where(mask == i)
    galaxy_segments = np.zeros_like(mask)
    galaxy_segments[inds] = 1
    mask[inds] = 0
    mask[np.where(mask != 0)] = 1
    data_with_mask = ma.array(image_data, mask=mask, fill_value = 0)
    fits.writeto(input_dir + 'data_with_mask.fits', data_with_mask.filled(), overwrite=True)


i = 0
params =  read_from_file()
RA, DEC = params[:,0], params[:,1]
input_dir = '/home/maria/Documents/projects/Stripe82/rdeep/EG_%.5f_%.5f/r/' % (RA[i],DEC[i])
os.chdir(input_dir)
image_file = input_dir + '/EG_%.5f_%.5f_r_.fits' % (RA[i],DEC[i])
image_data = get_image_data(image_file)

make_mask()