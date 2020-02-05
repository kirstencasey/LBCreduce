'''
Functions for:
	- fetching images
	- calculating basic statistics about images
	- displaying images in a notebook
'''
from ccdproc import ImageFileCollection
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Get specific images (input image type, object name, level of reduction, etc.)
def get_images(directory, imtype='*', objectname='*', filt='*', redlevel=None, display_options=False):
	'''
	This function returns a list of images given input parameters specifying the type of images needed.

	Parameters
	----------
	directory : str
		The directory in which getImages should look for .fits files. 

	imtype : str
		The IMAGETYP keyword in the header of the .fits files to be returned. Used to indicate whether 
		the returned images should be flats, biases, objects, etc.

	objectname : str
		The OBJECT keyword in the header of the .fits files to be returned. Used to indicate which
		target (eg. UGC1171, NGC628, etc.) is desired.
	
	filt : str
		The FILTER keyword in the header of the .fits files to be returned. Used to indicate which
		telescope filter (eg. B-BESSEL, R-BESSEL, etc.) is desired.

	redlevel : str
		The level of image reduction desired in the returned images. 

	display_options : bool
		If True, get_images will not return a list of images but will instead return a dict of available
		options for each other input parameter (imtype, objectname, redlevel) based on the contents of 
		the .fits headers in the given directory.

	Returns
	-------
	dict (when display_options == True)
		Returns dict containing information about available options for the input parameters given the .fits 
		files in the given directory

	list (when display_options == False)
		Returns a list of image filenames consistent with the input parameters	
	
	'''

	path = Path(directory)
	keys = ['imagetyp','object','filter','propid','exptime']
	all_images = ImageFileCollection(path, keywords=keys)		

	if display_options:
		images = all_images.filter(imagetyp=imtype, object=objectname, filter=filt)

		options = {}

		options['imtype'] = images.values('imagetyp', unique=True)
		options['objectname'] = images.values('object', unique=True)
		options['filt'] = images.values('filter', unique=True)

		print('Options','-------',sep='\n')		
		print('\n','imtype:', options['imtype'], sep='')
		print('\n', 'objectname:', options['objectname'], sep='')
		print('\n', 'filt:', options['filt'], sep='')

		return images
	
	else:

		images = all_images.files_filtered(include_path=True, imagetyp=imtype, object=objectname, filter=filt)

		return images	

# Calculate given image stats
def calc_stats(images):
	'''
	This function accepts a list of images and returns basic statistics for the image counts
	
	Parameters
	----------
	images : array of strings
		This is the list of filenames for which to calculate basic statistics

	Returns
	-------
	arr
		Returns information about the mean, median, std dev, etc. count for each input image
	
	'''
	stats = [] # Do stuff
	return stats

# Display given image (for use in notebooks)
def show_image(image):
	'''
	This function displays a given image in a notebook

	Parameters
	----------
	image : arr
		This is a 2D array representing the image to be displayed

	Returns
	-------
	None
	'''





