'''
Functions for allowing user to check output of processing pipeline:
	- Overscan and trim
	- Constructing master biases, flat fields, etc.
	- Stacking images
'''
import image
from astropy.stats import sigma_clip, mad_std
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

def initialize_directories(options, check_in_dir=True, check_out_dir=True):
	'''
	This function checks to make sure input and output directories are valid based on user input.

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

	# Check image_dir exists
	image_dir_exists = os.path.isdir(image_dir)
	if not image_dir_exists:
		warnings.warn('Given image_dir doesn\'t exist.', AstropyUserWarning)
		image_dir = get_input('Enter a new image_dir (where lbcreduce looks for images to use in reduction): ', is_dir=True)

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
			if response in affirmative:
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

def change_file_selection(options, zeros=False, darks=False, flats=False, objects=False, start_over=False):
	'''
	This function takes the reduction options as an input and asks the user for input
	about which options should be changed to affect the input files for reduction.

	Parameters
	----------
	options : dict
		Dictionary containing information about data reduction options, including how to select files for reduction.
		Necessary items in dictionary are:
			- image_dir : str
				Directory where raw images are saved
			- include_filenames : list
				List of filenames to be included in image processing (should include any necessary object images, flats, darks, zero frames, etc.)
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

	zeros : bool
		A boolean argument letting the function know that only the zero frames need be changed

	darks : bool
		A boolean argument letting the function know that only the dark frames need be changed

	flats : bool
		A boolean argument letting the function know that only the flat fields need be changed

	objects : bool
		A boolean argument letting the function know that only the object files need be changed
	'''
	# Select all new file collection
	if start_over:
		print('The available options that affect which files are selected for processing are: image_dir, include_filenames, glob_include, glob_exclude, and object.')
		# Check for input directory
		change_image_dir = get_input('Would you like to change image_dir (the directory where lbcreduce looks for images to process)? [y/n]: ')
		if change_image_dir:
			options['image_dir'] = get_input('Enter new image_dir: ', is_dir=True)
		# Check for list of files
		change_filelist = get_input('Would you like to specify a list of files to be included in image reduction (applied before glob_include; recommended only for a short list of files)? [y/n]: ')
		if change_filelist:
			options['include_filenames'] = get_input('Enter a list of filenames to be included in reduction: ')
		# Check for glob include/exclude
		change_glob_include = get_input('Would you like to specify a new glob_include (a Unix-style filename segment used to select included files)? [y/n]: ')
		if change_glob_include:
			options['glob_include'] = get_input('Enter new glob_include: ', anything_acceptable=True)
		change_glob_exclude = get_input('Would you like to specify a new glob_exclude (a Unix-style filename segment used to select excluded files)? [y/n]: ')
		if change_glob_exclude:
			options['glob_exclude'] = get_input('Enter new glob_exclude: ', anything_acceptable=True)
		# Check for specific object
		change_object = get_input('Would you like to specify a single object for image reduction? [y/n]: ')
		if change_object:
			options['object'] = get_input('Enter the object name as it appears in the OBJECT keyword of the header: ', anything_acceptable=True)

	else:
		if zeros:
			print('Okay, let\'s look for some zero frames to include in the analysis.')
			######################################################## DO STUFF HERE STILL ########################################################

		if darks:
			print('Okay, let\'s look for some dark frames to include in the analysis.')
			######################################################## DO STUFF HERE STILL ########################################################

		if flats:
			print('Okay, let\'s look for some flat fields to include in the analysis.')
			######################################################## DO STUFF HERE STILL ########################################################

		if objects:
			print('Okay, let\'s look for some objects to include in the analysis.')
			######################################################## DO STUFF HERE STILL ########################################################


	return options
