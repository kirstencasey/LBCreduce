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
from astropy.nddata import CCDData
from astropy.stats import sigma_clip, mad_std
from . import interactive, tools
import sys, glob, warnings, os
from astropy.modeling import models, fitting
from scipy.stats import chisquare


affirmative = ['y','Y','yes','Yes']
negative = ['n','N','no','No']
keys = ['imagetyp', 'object', 'propid']


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
		all_images = glob.glob(os.path.join(image_dir,config['glob_include']))

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
		hdulist = fits.open(os.path.join(image_dir, fi))
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
	save_dir = os.path.join(config['out_dir'],'midproc')

	# Loop through files
	for fi in image_info['filename']:

		# Get current HDU data
		hdul = fits.open(os.path.join(config['image_dir'], fi))
		primary_hdu = hdul[0]
		primary_hdu.header['NEXTEND'] = 2
		chip1_hdu = hdul[1]
		chip2_hdu = hdul[2]
		chip3_hdu = hdul[3]
		chip4_hdu = hdul[4]

		# Save selected chips
		if config['reduce_selected_chips']['chip1']:
			fits.HDUList([primary_hdu,chip1_hdu]).writeto(os.path.join(save_dir, fi.split('/')[-1].replace('.fits','-chip1.fits')))
		if config['reduce_selected_chips']['chip2']:
			fits.HDUList([primary_hdu,chip2_hdu]).writeto(os.path.join(save_dir, fi.split('/')[-1].replace('.fits','-chip2.fits')))
		if config['reduce_selected_chips']['chip3']:
			fits.HDUList([primary_hdu,chip3_hdu]).writeto(os.path.join(save_dir, fi.split('/')[-1].replace('.fits','-chip3.fits')))
		if config['reduce_selected_chips']['chip4']:
			fits.HDUList([primary_hdu,chip4_hdu]).writeto(os.path.join(save_dir, fi.split('/')[-1].replace('.fits','-chip4.fits')))

		hdul.close()

	files_by_chip = get_image_info(config, image_dir=os.path.join(config['out_dir'],'midproc'))

	return files_by_chip

def get_ccd_section(image_header, section_name):
	'''
	This function takes a string in the format, '[xmin:xmax,ymin:ymax]',
	and converts the appropriate str-type FITS indices to int-type python indices
	'''
	section = image_header[section_name]

	xmin = int(section.split('[')[1].split(':')[0]) - 1 # NOTE: x and y here are in the FITS convention
	xmax = int(section.split(':')[1].split(',')[0])
	ymin = int(section.split(',')[1].split(':')[0]) - 1
	ymax = int(section.split(':')[-1].split(']')[0])

	#if section_name == 'BIASSEC':	############ DO THIS BETTER!!!!!! ############
	#	xmin += 10

	return ymin, ymax, xmin, xmax # NOTE: To fit the python convention when this function is called, y and x are swapped

def check_flat_counts(flat_info, config):

	# Get flat counts
	counts = []
	for flat in flat_info['filename']:
		data = CCDData.read(os.path.join(config['out_dir'], 'midproc',  flat), unit=config['data_units'], hdu=config['ext'])
		counts.append(np.sum(data))
	counts = np.asarray(counts)

	# Throw out flats with bad counts
	clipped_counts = stats.sigma_clip(counts, **config['sigma_clip_flats_options'])
	mask = ~clipped_counts.mask

	return flat_info[mask].copy()


def find_best_masterframe(mastertype, sci_info, config, date = None):
	if sci_info['instrument'] == config['lbc_red']:
		instrument = '_R'
	elif sci_info['instrument'] == config['lbc_blue']:
		instrument = '_B'
	chip = '-' + sci_info['filename'].split('-')[-1].split('.fits')[0]

	if mastertype == 'zero':
		master_name = os.path.join(config['out_dir'], 'midproc', 'masterbias_zero' + chip + instrument + '.fits')############ CHANGE THIS TO CONSIDER DIFFERENT MASTER DATES ###########
		master_data = CCDData.read(master_name, unit=config['data_units'])

	elif mastertype == 'bias2D':
		master_name = os.path.join(config['out_dir'], 'midproc', 'masterbias_2Dbias' + chip + instrument + '.fits')############ CHANGE THIS TO CONSIDER DIFFERENT MASTER DATES ###########
		master_data = CCDData.read(master_name, unit=config['data_units'])

	elif mastertype == 'dark':
		master_name = os.path.join(config['out_dir'], 'midproc', 'masterdark' + chip + '_' + sci_info['filter'] + '.fits')############ CHANGE THIS TO CONSIDER DIFFERENT MASTER DATES ###########
		master_data = CCDData.read(master_name, unit=config['data_units'])

	elif mastertype == 'flat':
		master_name = os.path.join(config['out_dir'], 'midproc', 'masterflat' + '_' + date + chip + '_' + sci_info['filter'] + '.fits')
		'''
		if os.path.isfile(master_name):
			master_data = CCDData.read(master_name, unit=config['data_units'])
		else:
			warnings.warn(f'File not found: {master_name.split('/')[-1]}', AstropyUserWarning)
			warnings.warn('Select one of the following masterflats to use instead or type \'STOP\' to stop looking.', AstropyUserWarning)
		'''
		master_data = CCDData.read(master_name, unit=config['data_units'])
	return master_data



def find_best_overscan_legendre_model(overscan, image_info, legendre_options, max_degree_allowed, std_dev_func = mad_std, sigma=4, bin_width=12, bin_func='sigmaclip_mean', sigma_clip_options=None, save_plot=True, out_dir=None):

	# Correct bin_width, prep for binning overscan
	imgname = image_info['filename']
	bin_width = tools.find_nearest_multiple(overscan.shape[0], bin_width)
	num_rows = int(overscan.shape[0]/bin_width - 1)
	binned_overscan = np.zeros(num_rows)
	x = np.arange(num_rows)

	if save_plot: fig = plt.figure()

	# Bin overscan data
	idx = 0
	row_num = 0
	saverow = []

	for row in overscan.data:
		saverow = np.append(saverow,row)
		if row_num % bin_width is 0 and row_num is not 0:

			if bin_func == 'median':
				binned_overscan[idx] = np.median(saverow)
			elif bin_func == 'sigmaclipped_mean':
				binned_overscan[idx] = np.mean(sigma_clip(saverow, **sigma_clip_options))
			elif bin_finc == 'mean':
				binned_overscan[idx] = np.mean(saverow)

			saverow=[]
			idx += 1

		row_num += 1

	# Create overscan models up to and including max_degree_allowed (and with all the other parameters in legendre_options) based on binned data
	fit = fitting.LevMarLSQFitter()
	all_models = []
	chisquare_results = []
	model_num = 1
	while model_num <= max_degree_allowed:
		legendre_options['degree'] = model_num
		all_models.append(fit(models.Legendre1D(**legendre_options), x, binned_overscan))
		chisquare_results.append(chisquare(binned_overscan,all_models[-1](x), ddof=model_num).statistic)
		model_num += 1

	# Identify and deal with problem cases

	# Decide which overscan model is best (use lowest order possible)
	limiting_chi = np.median(chisquare_results) + sigma * std_dev_func(chisquare_results)
	best_model = all_models[np.where(chisquare_results < limiting_chi)[0][0]]

	# Plot overscan and models if requested

	if save_plot:
		plt.ylabel('counts')
		plt.xlabel('image row (binned)')
		plt.plot(binned_overscan,label=f'overscan, {bin_func}',c='black')
		plt.plot(x, best_model(x), label=f'Legendre1D, degree = {best_model.degree}',c='dodgerblue')
		plt.savefig(os.path.join(out_dir, f'overscan_model_{imgname}.png'))
		plt.close()

	return best_model
