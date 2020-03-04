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
from astropy.utils.exceptions import AstropyUserWarning
import interactive
import sys

affirmative = ['y','Y','yes','Yes']
negative = ['n','N','no','No']
keys = ['imagetyp', 'filter', 'object']

# Get images to do the reduction
def get_images(config):
	'''
	This function returns an ImageFileCollection of images given input parameters specifying the type of images needed.

	Parameters
	----------
	options : dict
		Options input for selecting images. Necessary items in dictionary are:
			- image_dir : str
				Directory where raw images are saved
			- ext : int
				Extension of FITS HDUList to read
			- filenames : list
				Filename(s) to be included in image processing (should include any flats, bias images, darks, etc. necessary for image reduction)
			- object : str
				OBJECT keyword in FITS header; filters given files to include only those with given keyword as well as any flats, darks, zero frames, etc.
			- glob_include : str
				Unix-style filename segment used to include files (eg. \'lbcb*\')
			- glob_exclude : str
				Unix-style filename segment used to exclude files (eg. \'lbcr*\')

	Returns
	-------
	all_images : ImageFileCollection
		Returns an ImageFileCollection object containing information about the images desired for reduction
	options : dict
		Dictionary containing all the same information as the input dictionary but with any options that were changed from user interaction updated
	'''
	all_images = ImageFileCollection(config['image_dir'], ext=config['ext'], keywords=keys, **config['file_selection_options'])

	# Get necessary files
	if config['object'] != None:
		all_files = all_images.files_filtered(object=config['object'])
		all_files += all_images.files_filtered(imagetyp='dark')
		all_files += all_images.files_filtered(imagetyp='flat')
		all_files += all_images.files_filtered(imagetyp='zero')
		all_images = ImageFileCollection(config['image_dir'], filenames=all_files, ext=config['ext'], keywords=keys)


	if type(all_images.summary) == type(None):
		warnings.warn('No files found in image_dir!',AstropyUserWarning)	########## WORK ON THIS - for when no files are found ###########
		sys.exit('lbcreduce stopped.')


	return all_images, config

def check_files(files, config):
	'''
	This function checks to make sure all image types needed to do the image analysis are available

	Parameters
	----------
	files : ImageFileCollection
		Files previously selected to be checked for validity given options.
	options : dict
		Options input for selecting images. Necessary items in dictionary are:
			- image_dir : str
				Directory where raw images are saved
			- ext : int
				Extension of FITS HDUList to read
			- do_zero : bool
				If True, check for bias images
			- do_dark : bool
				If True, check for dark frames
			- do_flat : bool
				If True, check for flat fields
			- object : str
				OBJECT keyword in FITS header; filters given files to include only those with given keyword as well as any flats, darks, zero frames, etc.
			- glob_include : str
				Unix-style filename segment used to include files (eg. \'lbcb*\')
			- glob_exclude : str
				Unix-style filename segment used to exclude files (eg. \'lbcr*\')
			- check_output : int
				Option specifying how much interaction with the program the user desires (0 is no interaction, 1 is some interaction, 2 is all interaction)

	Returns
	-------
	all_images : ImageFileCollection
		Returns an ImageFileCollection object containing information about the images desired for reduction
	options : dict
		Dictionary containing all the same information as the input dictionary but with any options that were changed from user interaction updated
	'''
	# Get list of all image types available
	config_changed = False
	imagetypes = files.values('imagetyp',unique=True)

	# Check zero frames
	if config['zero']:
		if 'zero' not in imagetypes:
			warnings.warn('No zero/bias frames detected in ImageFileCollection!',AstropyUserWarning)
			if check_output != 0:
				response = interactive.get_input('Change the file selection criteria before continuing? [y/n]: ')
				if response in affirmative:
					# Ask for new
					config = change_file_selection(config, zeros=True)
					config_changed = True
			else:
				sys.exit('lbcreduce stopped.')
	# Check for darks
	if config['dark']:
		if 'dark' not in imagetypes:
			warnings.warn('No dark frames detected in ImageFileCollection!',AstropyUserWarning)
			if check_output != 0:
				response = interactive.get_input('Change the file selection criteria before continuing? [y/n]: ')
				if response in affirmative:
					# Ask for new
					config = change_file_selection(options, darks=True)
					config_changed = True
			else:
				sys.exit('lbcreduce stopped.')
	# Check for flats
	if config['flat']:
		if 'flat' not in imagetypes:
			warnings.warn('No flat fields detected in ImageFileCollection!',AstropyUserWarning)
			if check_output != 0:
				response = interactive.get_input('Change the file selection criteria before continuing? [y/n]: ')
				if response in affirmative:
					# Ask for new
					config = change_file_selection(options, flats=True)
					config_changed = True
			else:
				sys.exit('lbcreduce stopped.')

	# Check for object files
	if config['reduce_objects']:
		if 'object' not in imagetypes:
			warnings.warn('No object files detected in ImageFileCollection!',AstropyUserWarning)
			if check_output != 0:
				response = interactive.get_input('Change the file selection criteria before continuing? [y/n]: ')
				if response in affirmative:
					# Ask for new
					config = change_file_selection(options, objects =True)
					config_changed = True
			else:
				sys.exit('lbcreduce stopped.')

	# Get new files from changed options
	if config_changed:
		all_images, config = get_images(config)

	return all_images, config

def get_ccd_section(image_header, section_name):
	'''
	This function takes a string in the format, '[xmin:xmax,ymin:ymax]', 
	and converts the appropriate str-type FITS indices to int-type python indices
	'''
	section = image_header[section_name]
	
	xmin = int(section.split('[')[1].split(':')[0]) - 1
	xmax = int(section.split(':')[1].split(',')[0])
	ymin = int(section.split(',')[1].split(':')[0]) - 1
	ymax = int(section.split(':')[-1].split(']')[0])

	if section_name == 'BIASSEC':	############ DO THIS BETTER!!!!!! ############
		xmin += 10
	
	return ymin, ymax, xmin, xmax

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
