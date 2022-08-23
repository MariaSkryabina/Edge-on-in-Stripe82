import astropy.io.fits as pyfits 
import numpy as np  
from glob import glob   
from astropy.wcs import WCS
from math import sqrt
import download_Stripe82 as dwnld
import os

def find_new_coord(index, xc, yc):
    dist = 2000*0.396/3600
    if index == 0:
        x_new = xc
        y_new = yc+dist
    if index == 1:
        x_new = xc
        y_new =  yc-dist 
    if index == 2:
        x_new = xc+dist
        y_new = yc
    if index == 3:
        x_new = xc-dist
        y_new = yc
    return x_new, y_new

def find_edge_distant(x2, y2, x1, y1):
    d = sqrt((x2-x1)**2+(y2-y1)**2)
    return d

def find_center_coord(RA, DEC, band, hdu_inp=0):

    image_dir = '/media/maria/TOSHIBA/images_final/EG_%.5f_%.5f/' % (RA,DEC) + band +'/'
    os.chdir(image_dir)
    if len(glob('*rec.fits')) !=1:
        print(glob('*rec.fits'))
        print(glob('*psf.fits'))
        i = float(input('which rec image to use?(index has to start with 0)'))
        input_image = glob('*rec.fits')[i]
    else: input_image = glob('*rec.fits')[0]

    
    hdu = pyfits.open(input_image)
    data = hdu[hdu_inp].data
    header = hdu[hdu_inp].header

    ny,nx = np.shape(data)
    print(ny,nx)
    w = WCS(hdu[0].header)
    print(RA, DEC)
    xc_pix, yc_pix = w.wcs_world2pix(RA, DEC, 1)
    os.chdir('/media/maria/TOSHIBA/images_final/')
    return(xc_pix, yc_pix, nx, ny )

def find_RA_DEC(): 
    file = open('/home/maria/Documents/projects/Stripe82/final_sample_coord.dat', 'r')
    text = file.readlines()
    text_float = np.zeros((2,len(text)))

    for string in range(len(text)):
        text_float[:,string] = [float(x) for x in text[string].split()[0:2]]
    file.close()
    RA = text_float[0,:]
    DEC = text_float[1,:]
    return(RA, DEC)


def main():
    RA, DEC = find_RA_DEC()
    bands = ['g', 'r', 'i']
    #bands = ['r']
#    for i in range(73,len(RA)):
    i = 1048
    for band in bands:

        xc_pix, yc_pix, nx, ny = find_center_coord(RA[i], DEC[i], band)
        xc_pix, yc_pix= int(xc_pix), int(yc_pix)

        print(xc_pix, yc_pix)

        x_min = 0;x_max = nx
        y_min = 0 ; y_max = ny
        d_up = find_edge_distant(xc_pix, yc_pix, xc_pix, y_max)
        d_down = find_edge_distant(xc_pix, yc_pix, xc_pix, y_min)
        d_left = find_edge_distant(xc_pix, yc_pix, x_min, yc_pix)
        d_right = find_edge_distant(xc_pix, yc_pix, x_max, yc_pix)
        print(d_up, d_down, d_left, d_right)

        d_array = [d_up, d_down, d_left, d_right]
        index_min = min(range(len(d_array)), key=d_array.__getitem__)
        print(index_min)

        RA_new, DEC_new = find_new_coord(index_min, RA[i], DEC[i])
        print(RA_new, DEC_new)
        dwnld.main(ra=[RA_new], dec=[DEC_new], name=['EG_%.5f_%.5f'  % (RA[i],DEC[i])],  bands=[band], frames=['rec'])

main()

