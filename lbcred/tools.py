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
from astropy.nddata import CCDData

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

def bias(bias_type, config, all_files):
	
	raw_bias_images = all_files.files_filtered(include_path=True, imagetyp=config['bias_image_keyword'])
	processed_bias_images = []

	# Make 2D bias image
	if bias_type == '2Dbias':
		# Loop through bias images
		for bias_im in raw_bias_images:
			# Get data 
			data = CCDData.read(bias_im, unit=config['data_units'])
			
			# Trim overscan region
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data_t = ccdproc.trim_image(data[xmin:xmax,ymin:ymax], **config['trim_options'])

			# Save trimmed image
			out_name = config['out_dir'] + 'midproc/' + bias_im.split('/')[-1].replace('.fits','_T.fits')
			data_t.write(out_name)
			processed_bias_images.append(out_name)
			
	# Make zero-frame bias image
	if bias_type == 'zero':
		# Loop through bias images
		for bias_im in raw_bias_images:
			# Get data 
			data = CCDData.read(bias_im, unit=config['data_units'])

			# Subtract overscan 
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['overscan_region'])
			data_o = ccdproc.subtract_overscan(data, overscan=data[xmin:xmax,ymin:ymax], **config['overscan_options'])

			# Trim overscan region
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data_ot = ccdproc.trim_image(data_o[xmin:xmax,ymin:ymax], **config['trim_options'])

			# Save overscan subtracted, trimmed image
			out_name = config['out_dir'] + 'midproc/' + bias_im.split('/')[-1].replace('.fits','_OT.fits')
			data_ot.write(out_name)
			processed_bias_images.append(out_name)

	# Combine images to make master bias
	master_name = config['out_dir'] + 'midproc/masterbias_' + bias_type + '.fits' 
	masterbias = ccdproc.combine(processed_bias_images, output_file=master_name, **config['combine_options'])

	# Get feedback on master bias

	return masterbias, config

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
