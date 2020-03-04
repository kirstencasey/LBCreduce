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
import ccdproc, os, sys, time, shutil, warnings, yaml
from astropy.utils.exceptions import AstropyUserWarning

def initialize_config(input_options, config_filename):
	# Open and read config file
	with open(config_filename, 'r') as filename:
		config = yaml.load(filename, Loader=yaml.FullLoader)

	# Replace options input via command line into config
	for key in input_options:
		if input_options[key] != None:
			config[key] = input_options[key]

	return config

def textfile(initial_config, end_notes, dir_overwritten, final_config = None):
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
	if final_config != None:
		out_dir = final_config['out_dir']
		image_dir = final_config['image_dir']
	else:
		out_dir = initial_config['out_dir']
		image_dir = initial_config['image_dir']
	file = open(out_dir+'image-production-details.txt','w')
	datetime = time.strftime('%Y-%m-%d at %H:%M:%S',time.gmtime())
	lines = [f'This file describes the options used to process the images in {image_dir}. \nThe image processing was performed by lbcreduce on {datetime}.',
			'\n\nNotes input by user at runtime:\n', initial_config['notes'], '\n\nNotes input by user after all image processing steps:\n', end_notes]

	# Make sure inputs are what were actually used if the user changed anything at all while running the program interactively
	file.writelines(lines)

	print('\nImage processing complete! :)\n')
	return

def bias(bias_type, options, all_files):

	# Get bias images (imagetyp: 'zero')
	bias_images = all_files.files_filtered(imagetyp='zero')

	# Make 2D bias image
	if bias_type == '2Dbias':
		# Trim overscan region


		# Combine images
		biasimage = ccdproc.combine(bias_images)

	# Make zero-frame bias image
	if bias_type == 'zero':
		# Subtract median (overscan) from all images

		# Trim overscan region

		# Combine images
		biasimage = ccdproc.combine(bias_images)

	return biasimage, options

# Overscan and trim
def overscan(options, all_files):
	'''

	'''
	# Get biases (image.get_images)

	# Loop through files:

		# Subtract overscan

		# Trim image

		# Save the result
	return

# Flat fielding
def flat(options, all_files):
	'''
	'''
	# Get flats (image.get_images)

	# Loop through files:

		#
	return options
# Combine calibrated images
def combine(options, all_files):
	'''
	'''
	# Combine calibrated bias images

	# Or, combine calibrated flats

		# Rescale individual flats

		# Combine flats
	return options
# Calibrate dark frames
def dark(options, all_files):
	'''
	'''
	return options

# Process images
def process(options, all_files):

	# Do stuff

	return options

# Stacking
def stack(options, all_files):
	'''
	Note: talk to Chris and Johnny about this
	'''
	return options
