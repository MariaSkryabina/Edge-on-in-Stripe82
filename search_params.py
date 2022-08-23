import os
import numpy as np
import glob
from decimal import Decimal
from math import sqrt
text_float = np.genfromtxt('/home/maria/Documents/projects/Stripe82/super_final_sample.txt')
RA = text_float[:,0]
DEC = text_float[:,1]

def return_cat_data(cat_file):
	file = open(cat_file, 'r')
	text = file.readlines()
	text_float = np.zeros((75,len(text)-75))
	for string in range(75,len(text)):
		text_float[:,string-75] = [float(x) for x in text[string].split()]
		file.close()

	return (text_float)


for i in range(len(RA)):

	path_to_data = '/home/maria/Documents/projects/Stripe82/rdeep/EG_%.5f_%.5f/rdeep/' % (RA[i],DEC[i])
	os.chdir(path_to_data)
	cat_file = glob.glob('*.cat')[0]

	cat_data = return_cat_data(cat_file)
	for string in range(np.array(cat_data).shape[1]):
		#print(float(round(Decimal(cat_data[0,4201]), 5)),RA[i])
		distance_ra = cat_data[0,string] - RA[i] 
		distance_dec = cat_data[1,string] - DEC[i]
		deviation = sqrt(distance_ra**2 +distance_dec**2)
		if deviation < 0.001:
		
			f = open('/home/maria/Documents/projects/Stripe82/params.cat', 'a', newline = '\n')
			f.write(str(i))
			f.write(' ')
			f.write(str(cat_data[:,string]).replace('\n','').replace('[','').replace(']','') +'\n')
			f.close()
			print(i,'is found', "%.6f" % (deviation))
			break
		else: 
			continue



