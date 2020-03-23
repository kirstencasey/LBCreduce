'''
Functions for:
	- fetching images
	- calculating basic statistics about images
	- displaying images in a notebook
'''
from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy import stats
from astropy.utils.exceptions import AstropyUserWarning
import interactive
import sys, glob, warnings

affirmative = ['y','Y','yes','Yes']
negative = ['n','N','no','No']
keys = ['imagetyp', 'object', 'propid']
'''
class File_Information():
	def __init__(self, config, glob_include=None, directory=None, midproc=False, **keys):
		# make table
		if glob_include = None:
			glob_include = config['glob_include']
		if directory is None:
			directory = config['image_dir']

		all_files = glob.glob(directory+glob_include)

		if len(all_files) == 0:
			warnings.warn('No files found in ' + directory, AstropyUserWarning)	########## WORK ON THIS - for when no files are found ###########
			sys.exit('lbcreduce stopped.')

		if not midproc and config['exclude'] is not None:
			all_files = [ im for im in all_files if config['exclude'] not in im ]

		# Loop through files to get relevant information
		propids = []
		objects = []
		imagetyps = []					######################### DO THIS BETTER!!!!!########################
		filters = []
		for fi in all_files:
			hdulist = fits.open(fi)
			propids.append(hdulist[0].header['PROPID'])
			objects.append(hdulist[1].header['OBJECT'])
			imagetyps.append(hdulist[1].header['IMAGETYP'])
			filters.append(hdulist[1].header['FILTER'])
			hdulist.close()

		propids = np.asarray(propids)
		all_files = np.asarray(all_files)
		objects = np.asarray(objects)
		imagetyps = np.asarray(imagetyps)
		filters = np.asarray(filters)

		file_info = Table()
		file_info['filename'] = all_files
		file_info['imagetyp'] = imagetyps
		file_info['object'] = objects
		file_info['filter'] = filters
		file_info['propid'] = propids
		self.file_info = file_info

	def values(self, col_name, unique=False):
		vals = np.asarray(self.file_info[col_name])
		if unique:
			vals = np.unique(vals)
		return vals

'''


# Get images to do the reduction
def get_image_info(config, image_dir=None, filenames=None):
	'''
	This function returns an ImageFileCollection of images given input parameters specifying the type of images needed.

	Parameters
	----------
	options : dict
		Options input for selecting images. Necessary items in dictionary are:
			- image_dir : str
				Directory where raw images are saved
			- object : str
				OBJECT keyword in FITS header; filters given files to include only those with given keyword as well as any flats, darks, zero frames, etc.
			- glob_include : str
				Unix-style filename segment used to include files (eg. \'lbcb*\')
			- exclude : str
				Filename segment used to exclude files (eg. \'lbcr\')

	Returns
	-------
	all_images : ImageFileCollection
		Returns an ImageFileCollection object containing information about the images desired for reduction
	'''
	if image_dir is None:
		image_dir = config['image_dir']

	if filenames is None:
		all_images = glob.glob(image_dir+config['glob_include'])

		if len(all_images) == 0:
			warnings.warn('No files found in image_dir!', AstropyUserWarning)	########## WORK ON THIS - for when no files are found ###########
			sys.exit('lbcreduce stopped.')

		if config['exclude'] is not None:
			all_images = [ im for im in all_images if config['exclude'] not in im ]

	else:
		all_images = filenames
	all_images = [im.split('/')[-1] for im in all_images]

	# Loop through files to get relevant information
	propids = []
	objects = []
	imagetyps = []					######################### DO THIS BETTER!!!!!########################
	filters = []
	lbcinst = []
	for fi in all_images:
		hdulist = fits.open(image_dir + fi)
		propids.append(hdulist[0].header['PROPID'])
		objects.append(hdulist[1].header['OBJECT'])
		imagetyps.append(hdulist[1].header['IMAGETYP'])
		filters.append(hdulist[1].header['FILTER'])
		lbcinst.append(hdulist[0].header['INSTRUME'])
		hdulist.close()

	propids = np.asarray(propids)
	all_images = np.asarray(all_images)
	objects = np.asarray(objects)
	imagetyps = np.asarray(imagetyps)
	filters = np.asarray(filters)
	lbcinst = np.asarray(lbcinst)
	'''																########### THIS DOESN'T WORK!!!! NEED TO NOT EXCLUDE FLATS, BIAS IMAGES, ETC.
	if config['propid'] is not None:
		ims_to_keep = np.where(propids==config['propid'])
		propids = propids[ims_to_keep]
		all_images = all_images[ims_to_keep]
		objects = objects[ims_to_keep]
		imagetyps = imagetyps[ims_to_keep]

	if config['object'] is not None:
		ims_to_keep = np.where(objects==config['object'])
		propids = propids[ims_to_keep]
		all_images = all_images[ims_to_keep]
		objects = objects[ims_to_keep]
		imagetyps = imagetyps[ims_to_keep]
	'''
	image_info = Table()
	image_info['filename'] = all_images
	image_info['imagetyp'] = imagetyps
	image_info['object'] = objects
	image_info['propid'] = propids
	image_info['filter'] = filters
	image_info['instrument'] = lbcinst

	return image_info

def check_files(file_info, config):
	'''
	This function checks to make sure all image types needed to do the image analysis are available

	Parameters
	----------
	files : Astropy Table
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
	all_images : Astropy Table
		Returns an ImageFileCollection object containing information about the images desired for reduction
	options : dict
		Dictionary containing all the same information as the input dictionary but with any options that were changed from user interaction updated
	'''
	# Get list of all image types available
	imagetypes = file_info['imagetyp']

	# Check zero frames
	if config['zero']:
		if 'zero' not in imagetypes:
			warnings.warn('No zero/bias frames detected in ImageFileCollection!',AstropyUserWarning)
			sys.exit('lbcreduce stopped.')
	# Check for darks
	if config['dark']:
		if 'dark' not in imagetypes:
			warnings.warn('No dark frames detected in ImageFileCollection!',AstropyUserWarning)
			sys.exit('lbcreduce stopped.')
	# Check for flats
	if config['flat']:
		if 'flat' not in imagetypes:
			warnings.warn('No flat fields detected in ImageFileCollection!',AstropyUserWarning)
			sys.exit('lbcreduce stopped.')

	# Check for object files
	if config['reduce_objects']:
		if 'object' not in imagetypes:
			warnings.warn('No object files detected in ImageFileCollection!',AstropyUserWarning)
			sys.exit('lbcreduce stopped.')

	return

def separate_chips(image_info, config):
	save_dir = config['out_dir']+'midproc/'

	# Loop through files
	for fi in image_info['filename']:

		# Get current HDU data
		hdul = fits.open(config['image_dir'] + fi)
		primary_hdu = hdul[0]
		primary_hdu.header['NEXTEND'] = 2
		chip1_hdu = hdul[1]
		chip2_hdu = hdul[2]
		chip3_hdu = hdul[3]
		chip4_hdu = hdul[4]

		# Save each chip
		fits.HDUList([primary_hdu,chip1_hdu]).writeto(save_dir + fi.split('/')[-1].replace('.fits','-chip1.fits'))
		fits.HDUList([primary_hdu,chip2_hdu]).writeto(save_dir + fi.split('/')[-1].replace('.fits','-chip2.fits'))
		fits.HDUList([primary_hdu,chip3_hdu]).writeto(save_dir + fi.split('/')[-1].replace('.fits','-chip3.fits'))
		fits.HDUList([primary_hdu,chip4_hdu]).writeto(save_dir + fi.split('/')[-1].replace('.fits','-chip4.fits'))

		hdul.close()

	files_by_chip = get_image_info(config, image_dir=config['out_dir']+'midproc/')

	return files_by_chip

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

def check_flat_counts(flat_info, config):
	return flat_info.copy()


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
