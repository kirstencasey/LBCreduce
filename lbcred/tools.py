'''
Functions for reducing LBC images:
	- Overscan and trim
	- Constructing master biases
	- Calibrating flat fields, etc.
	- Stacking images
'''
import numpy as np
import image
from astropy.stats import sigma_clip, mad_std
import ccdproc, os, sys, time, shutil, warnings
from astropy.utils.exceptions import AstropyUserWarning

def textfile(initial_options, notes, end_notes, dir_overwritten, final_options = None):
	'''
	This function creates a textfile in the directory containing the processed images which details all the decisions made during processing

	Parameters
	----------
	initial_options : list
		Options input by user when the program was first called.

	final_options : list
		Options actually used. If the options were not changed during the course of the program then final_options is None.

	end_notes : str
		Notes input by user after all the image processing steps have complted to be saved in textfile.

	dir_overwritten : bool
		Information about whether a directory was overwritten in order to make room for the output directory

	Returns
	-------
	N/A
	'''

	# Make text file with notes
	if final_options != None:
		out_dir = final_options['out_dir']
		image_dir = final_options['image_dir']
	else:
		out_dir = initial_options['out_dir']
		image_dir = initial_options['image_dir']
	file = open(out_dir+'image-production-details.txt','w')
	datetime = time.strftime('%Y-%m-%d at %H:%M:%S',time.gmtime())
	lines = [f'This file describes the options used to process the images in {image_dir}. \nThe image processing was performed by lbcreduce on {datetime}.',
			'\n\nNotes input by user at runtime:\n', notes, '\n\nNotes input by user after all image processing steps:\n', end_notes]

	# Make sure inputs are what were actually used if the user changed anything at all while running the program interactively
	file.writelines(lines)

	print('\nImage processing complete! :)\n')
	return

def initialize_directories(options):
	'''
	This function creates a textfile in the directory containing the processed images which details all the decisions made during processing

	Parameters
	----------
	options : dict
		Dictionary containing information about how directories should be initialized. Necessary items in dictionary are:
			- image_dir : str
				Directory where raw images are saved
			- out_dir : str
				Directory where processed images will be saved
			- overwrite : bool
				Option specifying whether an existing directory should be overwritten if it has the same name/path as out_dir
			- check_output : int
				Option specifying how much interaction with the program the user desires (0 is no interaction, 1 is some interaction, 2 is all interaction)

	Returns
	-------
	options : dict
		Dictionary containing all the same information as the input dictionary but with any options that were changed from user interaction updated

	dir_overwritten : bool
		Information on whether it was necessary to overwrite an existing directory
	'''
	# Get necessary information from options
	image_dir = options['image_dir']
	out_dir = options['out_dir']
	overwrite = options['overwrite']
	check_output = options['check_output']

	# Fix directory formats
	image_dir = os.path.abspath(image_dir) + '/'
	if out_dir != None:
		out_dir = os.path.abspath(out_dir) + '/'
	else:
		datetime = time.strftime('%Y%m%d_%H%M%S',time.gmtime())
		out_dir = image_dir + f'lbcreduce_{datetime}/'
		out_dir = os.path.abspath(image_dir + f'../lbcreduce_{datetime}/') + '/'

	directory_exists = os.path.isdir(out_dir)
	dir_overwritten = False

	if directory_exists and overwrite:
		warnings.warn('The output directory already exists and will be overwritten.',AstropyUserWarning)
		if check_output != 0:
			response = input('Are you sure you want to proceed? [y/n]: ')
			if response == 'y' or response == 'yes':
				shutil.rmtree(out_dir)
				os.mkdir(out_dir)
				dir_overwritten = True
			else:
				print('Output directory already exists and lbcreduce was not given permission to overwrite it.')
				overwrite = False
				sys.exit('lbcreduce stopped.')
	elif directory_exists and not overwrite:
		print('Error: Output directory already exists and lbcreduce was not given permission to overwrite it.')
		print('Either specify a different out_dir or set overwrite to True.')
		sys.exit('lbcreduce stopped.')

	else:
		os.mkdir(out_dir)

	os.mkdir(out_dir+'midproc/')

	# Update options
	options['image_dir'] = image_dir
	options['out_dir'] = out_dir
	options['overwrite'] = overwrite
	options['check_output'] = check_output

	return options, dir_overwritten

def bias(bias_type, options):

	# Make 2D bias image
	if bias_type == '2Dbias':


	# Make zero-frame bias image
	if bias_type == 'zero':

	return biasimage, options

# Overscan and trim
def overscan(options):
	'''

	'''
	# Get biases (image.get_images)

	# Loop through files:

		# Subtract overscan

		# Trim image

		# Save the result
	return

# Flat fielding
def flat(options):
	'''
	'''
	# Get flats (image.get_images)

	# Loop through files:

		#
	return options
# Combine calibrated images
def combine(options):
	'''
	'''
	# Combine calibrated bias images

	# Or, combine calibrated flats

		# Rescale individual flats

		# Combine flats
	return options
# Calibrate dark frames
def dark(options):
	'''
	'''
	return options

# Process images
def process(options):

	# Do stuff

	return options

# Stacking
def stack(options):
	'''
	Note: talk to Chris and Johnny about this
	'''
	return options
