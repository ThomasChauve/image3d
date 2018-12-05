import image3d.image3d as image3d
import numpy as np
import os
import re
from PIL import Image
import time
import matplotlib.pyplot as plt

def sorted_aphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

def load_image_from_tiff(adr_folder,resolution):
	'''
	Load image from a Folder containing N tiff of N*N dimention
	'''
	output=sorted_aphanumeric(os.listdir(adr_folder))
	
	s3=len(output)


	ss=np.shape(Image.open(adr_folder+output[0]))
	
	im=np.zeros([ss[0],ss[1],s3])	

	for i in list(range(s3)):
		im[:,:,i]=Image.open(adr_folder+output[i])

	return image3d.image3d(im,resolution)

def make_autocorrelation(adr_folder,output,resolution,Cinf_coeff):
    # Open image
    print('Begining construction image', time.asctime())
    image=load_image_from_tiff(adr_folder,resolution)
    print('End construction image', time.asctime())
    
    # Make autocorrelation
    print ('Begining autocorrelation', time.asctime())
    Autocor=image.xcorr3d()
    Autocor.Cinf=Cinf_coeff*Autocor.Cinf
    print ('End autocorrelation', time.asctime())
    
    print ('Begining stereo proj', time.asctime())
    plt.figure(figsize=(20,20),dpi=160)
    Autocor.stereographic_corr_length()
    plt.xlabel('Stereoprojection - radius length')
    plt.savefig(output + '.png')
    print ('End stereo proj', time.asctime())  