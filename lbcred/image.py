'''
Functions for:
	- fetching images
	- calculating basic statistics about images
	- displaying images in a notebook
'''
from ccdproc import ImageFileCollection
import numpy as np
import matplotlib.pyplot as plt
from astropy import stats

# Get images to do the reduction
def get_images(options):
	'''
	This function returns a list of images given input parameters specifying the type of images needed.

	Parameters
	----------
	options : dict
		Options input for selecting images. Necessary items in dictionary are:
			- image_dir : str
				Directory where raw images are saved
			- filenames : str or list of str
				Filename(s) to be included in image processing (should include any flats, bias images, darks, etc. necessary for image reduction)
			- glob_include : str
				Unix-style filename segment used to include files (eg. \'lbcb*\')
			- glob_exclude : str
				Unix-style filename segment used to exclude files (eg. \'lbcr*\')
			- include_attributes : dict
				Dictionary describing elements in the header of FITS files with values to include from image_dir
			- exclude_attributes : dict
				Dictionary describing elements in the header of FITS files with values to exclude from image_dir  (does nothing if include_attributes is specified; cannot exclude flats, bias images, darks, etc.)
			- check_output : int
				Option specifying how much interaction with the program the user desires (0 is no interaction, 1 is some interaction, 2 is all interaction)

	Returns
	-------
	all_images : ImageFileCollection
		Returns an ImageFileCollection object containing information about the images desired for reduction
	options : dict
		Dictionary containing all the same information as the input dictionary but with any options that were changed from user interaction updated
	'''

	# Get necessary images
	all_images = ImageFileCollection(options['image_dir'], filenames=options['filenames'], glob_include=options['glob_include'])



	# Check that expected images are available
	return all_images, options


# Calculate given image stats
def calc_stats(images, sigma_clip=True, sigma=3.0, maxiters=5):
	'''
	This function accepts a list of images and returns basic statistics for the image counts

	Parameters
	----------
	images : ImageFileCollection
		This is an ImageFileCollection object of images for which to calculate basic statistics

	Returns
	-------
	arr
		Returns information about the mean, median, std dev, etc. count for each input image

	'''
	num_images = len(images.summary)

	mean = numpy.zeros(num_images)
	median = numpy.zeros(num_images)
	mad_std = numpy.zeros(num_images)
	std = numpy.zeros(num_images)
	counts = numpy.zeros(num_images)

	i = 0
	for im in images.data(do_not_scale_image_data=True):
		data_sigclipped = stats.sigma_clip(im, sigma=sigma, maxiters=maxiters)
		mean[i] = data_sigclipped.mean()
		median[i] = data_sigclipped.median()
		std[i] = data_sigclipped.std()
		mad_std = stats.mad_std(data_sigclipped)
		i+=1

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
	return
