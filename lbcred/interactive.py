'''
Functions for allowing user to check output of processing pipeline:
	- Overscan and trim
	- Constructing master biases, flat fields, etc.
	- Stacking images
'''
from . import image
from astropy.stats import mad_std
import numpy as np
import ccdproc, os, sys, time, shutil, warnings, yaml
from astropy.utils.exceptions import AstropyUserWarning

affirmative = ['y','Y','yes','Yes']
negative = ['n','N','no','No']
default_acceptable = affirmative + negative

def get_input(question, acceptable_responses=default_acceptable, anything_acceptable=False, exit_response='quit', is_dir=False):

	# Keep asking for input until acceptable input is entered
	response = None
	# For
	if not is_dir and not anything_acceptable:
		while response not in acceptable_responses and response != exit_response:
			response = input(question)
			if response not in acceptable_responses and response != exit_response:
				warnings.warn(f'Response invalid. Try again or exit lbcreduce by entering \'{exit_response}\'.', AstropyUserWarning)

	elif is_dir:
		dir_okay = False
		while not dir_okay and response != exit_response:
			response = input(question)
			if response != exit_response:
				dir_okay = os.path.isdir(response)
				if not dir_okay:
					warnings.warn(f'Input directory not found. Try again or exit lbcreduce by entering \'{exit_response}\'.', AstropyUserWarning)

	else:
		response = input(question)

	# Stop program if desired
	if response == exit_response:
		sys.exit('lbcreduce stopped.')

	return response

def initialize_directories(config, check_in_dir=True, check_out_dir=True):
	'''
	This function checks to make sure input and output directories are valid based on user input.

	Parameters
	----------
	config : dict
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
	image_dir = config['image_dir']
	out_dir = config['out_dir']
	overwrite = config['overwrite']
	check_output = config['check_output']

	# Check image_dir exists
	image_dir_exists = os.path.isdir(image_dir)
	if not image_dir_exists:
		warnings.warn('Given image_dir doesn\'t exist.', AstropyUserWarning)
		image_dir = get_input('Enter a new image_dir (where lbcreduce looks for images to use in reduction): ', is_dir=True)

	# Fix directory formats
	if out_dir is None:
		datetime = time.strftime('%Y%m%d_%H%M%S',time.gmtime())
		out_dir = os.path.join(image_dir, f'../lbcreduce_{datetime}')

	directory_exists = os.path.isdir(out_dir)
	dir_overwritten = False

	if directory_exists and overwrite:
		warnings.warn('The output directory already exists and will be overwritten.', AstropyUserWarning)
		shutil.rmtree(out_dir)
		os.mkdir(out_dir)
		dir_overwritten = True

	elif directory_exists and not overwrite:
		print('Error: Output directory already exists and lbcreduce was not given permission to overwrite it.')
		print('Either specify a different out_dir or set overwrite to True.')
		sys.exit('lbcreduce stopped.')

	else:
		os.mkdir(out_dir)

	os.mkdir(os.path.join(out_dir,'midproc'))
	os.mkdir(os.path.join(out_dir,'plots'))

	# Update options			################## DO THIS BETTER - IN A DIFFERENT FUNCTION ####################
	config['image_dir'] = image_dir
	config['out_dir'] = out_dir

	return config, dir_overwritten
