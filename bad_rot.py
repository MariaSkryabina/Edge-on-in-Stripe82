import numpy as np 
from glob import glob 
import os
import shutil
import subprocess

def find_RA_DEC(): 
    file = open('/home/maria/Documents/projects/Stripe82/final_sample_coord.dat', 'r')
    text = file.readlines()
    text_float = np.zeros((2,len(text)))

    for string in range(0,len(text)):
        text_float[:,string] = [float(x) for x in text[string].split()[0:2]]
    file.close()
    RA = text_float[0,:]
    DEC = text_float[1,:]
    return(RA, DEC)

def open_fits(RA, DEC, i):
    image_dir = '/media/maria/TOSHIBA/images_final/EG_%.5f_%.5f/' % (RA,DEC) +'r/'
    png_dir = '/media/maria/TOSHIBA/png/'
    os.chdir(image_dir)
    print('I am in the dirrectory of EG_%.5f_%.5f' % (RA,DEC))
    output_image = 'EG_%.5f_%.5f_r_crop_rot.fits' % (RA,DEC)
    subprocess.call("ds9 %s -scale histequ -zoom 4" % (output_image), shell=True)

    answer = str(input('fine rotation?') or 'y')

    file_num = open('/home/maria/Documents/projects/Stripe82/rot_again_num.dat', 'a')

    if answer == 'n':
        file_num.write(str(i)+'\n')

    file_num.close()

RA, DEC = find_RA_DEC()
for i in range (1045, len (RA)):
    open_fits(RA[i], DEC[i], i)
