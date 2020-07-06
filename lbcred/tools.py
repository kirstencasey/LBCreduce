'''
Functions for reducing LBC images:
	- Overscan and trim
	- Constructing master biases, flats
	- Calibrating flat fields, etc.
	- Stacking images
'''
import numpy as np
from . import image
from astropy.stats import sigma_clip, mad_std
import ccdproc, os, sys, time, shutil, warnings, yaml, random
from astropy.utils.exceptions import AstropyUserWarning
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Column
import logging
from .log import logger
from astropy.modeling import models
from lbcred import project_dir


def setup_logger(level, log_fn=None):
    """
    Setup the pipeline logger.
    Parameters
    ----------
    level : str
        The log level (debug, info, warn, or error).
    log_fn : str (optional)
       Log file name.
    """
    if log_fn is not None:
        fh = logging.FileHandler(log_fn)
        formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s: %(message)s',
            '%Y-%m-%d | %H:%M:%S')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    logger.setLevel(level.upper())

def inv_median(arr):
	return 1/np.median(arr)

def initialize_config(config_filename, input_options = {}):
	# Open and read config file
	with open(config_filename, 'r') as filename:
		config = yaml.load(filename, Loader=yaml.FullLoader)

	# Replace options input via command line into config
	for key in input_options:
		if input_options[key] != None:
			config[key] = input_options[key]

	chips = []
	if config['reduce_selected_chips']['chip1']: chips.append('-chip1')
	if config['reduce_selected_chips']['chip2']: chips.append('-chip2')
	if config['reduce_selected_chips']['chip3']: chips.append('-chip3')
	if config['reduce_selected_chips']['chip4']: chips.append('-chip4')
	config['chips'] = chips

	for typ in ['default','object','flat','bias','dark']:
		window = [ config[f'{typ}_options']['legendre_options']['window_lower'], config[f'{typ}_options']['legendre_options']['window_upper'] ]
		config[f'{typ}_options']['legendre_options'].pop('window_lower')
		config[f'{typ}_options']['legendre_options'].pop('window_upper')
		config[f'{typ}_options']['legendre_options']['window'] = window

	combine_arg_dict = {'median': np.ma.median,
						'mean' : np.ma.mean,
						'std' : np.ma.std,
						'mad_std' : mad_std,
						'inverse_median' : inv_median,
						None : None}

	for typ in ['default','object','flat','bias','dark']:
		config[f'{typ}_options']['combine_options']['sigma_clip_func'] = combine_arg_dict[config[f'{typ}_options']['combine_options']['sigma_clip_func']]
		config[f'{typ}_options']['combine_options']['sigma_clip_dev_func'] = combine_arg_dict[config[f'{typ}_options']['combine_options']['sigma_clip_dev_func']]
		config[f'{typ}_options']['combine_options']['scale'] = combine_arg_dict[config[f'{typ}_options']['combine_options']['scale']]
		config[f'{typ}_options']['variable_legendre_options']['sigma_clip_options']['stdfunc'] = combine_arg_dict[config[f'{typ}_options']['variable_legendre_options']['sigma_clip_options']['stdfunc']]
		config[f'{typ}_options']['variable_legendre_options']['std_dev_func'] = combine_arg_dict[config[f'{typ}_options']['variable_legendre_options']['std_dev_func']]
		if config[f'{typ}_options']['overscan_options']['model'] == 'legendre':
			config[f'{typ}_options']['overscan_options']['model'] = models.Legendre1D(**config[f'{typ}_options']['legendre_options'])

	return config

def textfile(config, end_notes, dir_overwritten):
	'''
	This function creates a textfile in the directory containing the processed images which details all the decisions made during processing
	'''

	# Make text file with notes
	out_dir = config['out_dir']
	image_dir = config['image_dir']
	filename = os.path.join(out_dir,'image-production-details.txt')
	file = open(filename,'w')
	datetime = time.strftime('%Y-%m-%d at %H:%M:%S',time.gmtime())
	lines = [f'This file describes the options used to process the images in {image_dir}. \nThe image processing was performed by lbcreduce on {datetime}.',
			'\n\nNotes input by user at runtime:\n', config['notes'], '\n\nNotes input by user after all image processing steps:\n', end_notes]

	# Make sure inputs are what were actually used if the user changed anything at all while running the program interactively
	file.writelines(lines)

	print('\nImage processing complete! :)\n')
	return

def get_config_options(image_type, config):

	if config[f'use_specified_{image_type}_options']:
		combine_options = config[f'{image_type}_options']['combine_options']
		overscan_options = config[f'{image_type}_options']['overscan_options']
		legendre_options = config[f'{image_type}_options']['legendre_options']
		if overscan_options['model'] is not None:
			if overscan_options['model'].__class__.name == 'Legendre1D':
				use_variable_legendre = config[f'{image_type}_options']['allow_variable_legendre']
				variable_legendre_options = config[f'{image_type}_options']['variable_legendre_options']
		else: use_variable_legendre = False

	else:
		combine_options = config['default_options']['combine_options']
		overscan_options = config['default_options']['overscan_options']
		legendre_options = config['default_options']['legendre_options']
		if overscan_options['model'] is not None:
			if overscan_options['model'].__class__.name == 'Legendre1D':
				use_variable_legendre = config['default_options']['allow_variable_legendre']
				variable_legendre_options = config['default_options']['variable_legendre_options']
		else: use_variable_legendre = False

	return combine_options, overscan_options, legendre_options, use_variable_legendre, variable_legendre_options

def bias(bias_type, config, file_info):

	raw_bias_info = file_info[np.where(file_info['imagetyp']==config['bias_image_keyword'])]
	processed_bias_info = raw_bias_info.copy()
	out_names = []

	combine_options, overscan_options, legendre_options, use_variable_legendre, variable_legendre_options = get_config_options('bias', config)

	# Make 2D bias image
	if bias_type == '2Dbias':
		# Loop through bias images
		for bias_im in raw_bias_info:
			# Get data
			data = CCDData.read(os.path.join(config['out_dir'], 'midproc', bias_im['filename']), unit=config['data_units'], hdu=config['ext'])

			# Trim overscan region
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data_trimmed = ccdproc.trim_image(data[xmin:xmax,ymin:ymax])
			data_trimmed.meta[config['data_region']] = f'[1:{data_trimmed.shape[1]},1:{data_trimmed.shape[0]}]'

			# Save trimmed image
			out_names.append(os.path.join(config['out_dir'], 'midproc', bias_im['filename'].replace('.fits','_T.fits')))
			primary_hdu = fits.open(os.path.join(config['out_dir'], 'midproc',  bias_im['filename']),ignore_blank=True)[0]
			image_hdu = fits.ImageHDU(data=data_trimmed,header=data_trimmed.header)
			fits.HDUList([primary_hdu,image_hdu]).writeto(out_names[-1])


	# Make zero-frame bias image
	if bias_type == 'zero':
		# Loop through bias images
		for bias_im in raw_bias_info:
			# Get data
			data = CCDData.read(os.path.join(config['out_dir'], 'midproc', bias_im['filename']), unit=config['data_units'], hdu=config['ext'])

			# Subtract overscan   ################################ TRYING NEW FUNC HERE ################################
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['overscan_region'])
			if use_variable_legendre:
				new_model = image.find_best_overscan_legendre_model(data[xmin:xmax,ymin:ymax], bias_im, legendre_options, out_dir=os.path.join(config['out_dir'], 'plots'), **variable_legendre_options)
				overscan_options['model'] = new_model
			data_o = ccdproc.subtract_overscan(data, overscan=data[xmin:xmax,ymin:ymax], **overscan_options)

			# Trim overscan region
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data_ot = ccdproc.trim_image(data_o[xmin:xmax,ymin:ymax])
			data_ot.meta[config['data_region']] = f'[1:{data_ot.shape[1]},1:{data_ot.shape[0]}]'

			# Save overscan subtracted, trimmed image
			out_names.append(os.path.join(config['out_dir'], 'midproc', bias_im['filename'].replace('.fits','_OT.fits')))
			primary_hdu = fits.open(os.path.join(config['out_dir'], 'midproc',  bias_im['filename']),ignore_blank=True)[0]
			image_hdu = fits.ImageHDU(data=data_ot,header=data_ot.header)
			fits.HDUList([primary_hdu,image_hdu]).writeto(out_names[-1])

	processed_bias_info['filename'] = Column(data=out_names, name='filename')

	# Combine images to make master bias
	for chip in config['chips']:
		master_name = os.path.join(config['out_dir'], 'midproc', 'masterbias_' + bias_type + chip)

		# Get correct chip
		mask = [idx for idx,fi in enumerate(processed_bias_info['filename']) if chip in fi]
		chip_info = processed_bias_info[mask]
		chip_info_R = chip_info[np.where(chip_info['instrument']==config['lbc_red'])]
		chip_info_B = chip_info[np.where(chip_info['instrument']==config['lbc_blue'])]

		# Sort by instrument
		masterbias_R = ccdproc.combine(chip_info_R['filename'], output_file=master_name+'_R.fits', unit=config['data_units'], **combine_options)
		masterbias_B = ccdproc.combine(chip_info_B['filename'], output_file=master_name+'_B.fits', unit=config['data_units'], **combine_options)

	# Get feedback on master bias

	return


# Flat fielding
def flat(config, file_info):
	'''
	'''
	# Get flats
	raw_flats_info = file_info[np.where(file_info['imagetyp']==config['flat_image_keyword'])]
	out_names = []
	dates = []

	# Examine flatfield counts - throw out flats with bad counts
	if config['examine_flat_counts']:
		processed_flats_info = image.check_flat_counts(raw_flats_info, config)
	else: processed_flats_info = raw_flats_info.copy()

	combine_options, overscan_options, legendre_options, use_variable_legendre, variable_legendre_options = get_config_options('flat', config)

	# Loop through files to calibrate flats (subtract bias, trim overscan):
	for flat_im in processed_flats_info:
		# Get data
		data = CCDData.read(os.path.join(config['out_dir'], 'midproc',  flat_im['filename']), unit=config['data_units'], hdu=config['ext'])
		dates.append(data.meta['DATE_OBS'].split('T')[0])

		# Trim overscan region
		xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['overscan_region'])   ############################# Use bias or overscan??? ################
		if use_variable_legendre:
			new_model = image.find_best_overscan_legendre_model(data[xmin:xmax,ymin:ymax], flat_im, legendre_options, out_dir=os.path.join(config['out_dir'], 'plots'), **variable_legendre_options)
			overscan_options['model'] = new_model
		data_o = ccdproc.subtract_overscan(data, overscan=data[xmin:xmax,ymin:ymax], **overscan_options)

		xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])  ############################# Use bias or overscan??? ################
		data_ot = ccdproc.trim_image(data_o[xmin:xmax,ymin:ymax])
		data_ot.meta[config['data_region']] = f'[1:{data_ot.shape[1]},1:{data_ot.shape[0]}]'


		'''
		xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])  ############################# Use bias or overscan??? ################
		data_t = ccdproc.trim_image(data[xmin:xmax,ymin:ymax])
		data_t.meta[config['data_region']] = f'[1:{data_t.shape[1]},1:{data_t.shape[0]}]'

																				############################# Use bias or overscan??? ################
		# Subtract 2D bias image
		if flat_im['instrument'] == config['lbc_red']:
			inst_color = '_R'
		elif flat_im['instrument'] == config['lbc_blue']:
			inst_color = '_B'
		chip = '-' + flat_im['filename'].split('-')[-1].split('.fits')[0]

		bias_name = os.path.join(config['out_dir'], 'midproc', 'masterbias_2Dbias' + chip + inst_color + '.fits')
		bias2D = CCDData.read(bias_name, unit=config['data_units'])
		data_ot = CCDData.subtract(data_t, bias2D)
		data_ot.meta = data_t.meta
		'''
		# Save calibrated flat
		out_names.append(os.path.join(config['out_dir'], 'midproc', flat_im['filename'].replace('.fits','_OT.fits')))
		primary_hdu = fits.open(os.path.join(config['out_dir'], 'midproc', flat_im['filename']),ignore_blank=True)[0]
		image_hdu = fits.ImageHDU(data=data_ot,header=data_ot.header)
		fits.HDUList([primary_hdu,image_hdu]).writeto(out_names[-1])

	processed_flats_info['filename'] = Column(data=out_names, name='filename')

	# Make master flat for each filter, day
	filts = np.unique(processed_flats_info['filter'])
	dates = np.unique(np.asarray(dates))
	for filt in filts:
		flats_to_combine = processed_flats_info['filename'][np.where(processed_flats_info['filter']==filt)].copy()
		for day in dates:
			date = day.split('-')[0] + day.split('-')[1] + day.split('-')[2]     ################# DO THIS BETTER
			mask = [idx for idx,fi in enumerate(flats_to_combine) if date in fi]
			date_info = flats_to_combine[mask]
			if len(date_info) == 0: continue
			for chip in config['chips']:
				master_name = os.path.join(config['out_dir'], 'midproc', 'masterflat' + '_' + day + chip + '_' + filt + '.fits')
				mask = [idx for idx,fi in enumerate(date_info) if chip in fi]
				chip_info = date_info[mask]
				if len(chip_info) == 0: continue
				if len(chip_info) == 1:
					shutil.copyfile(chip_info[0], master_name)
					warnings.warn(f'Only one flat to \'combine\' for {master_name}.', AstropyUserWarning)
					continue
				masterflat = ccdproc.combine(chip_info, output_file=master_name, unit=config['data_units'], **combine_options)

	# Get feedback on master flats

	return


# Calibrate dark frames
def dark(config, file_info):
	'''
	'''
	return

# Process images
def process(config, file_info):

	# Get images
	sci_ims = file_info[np.where(file_info['imagetyp']==config['object_image_keyword'])]
	out_dir  = config['out_dir']
	if config['stack']:
		out_dir = os.path.join(out_dir, 'midproc')

	combine_options, overscan_options, legendre_options, use_variable_legendre, variable_legendre_options = get_config_options('object', config)

	# Loop through images
	for im in sci_ims:
		# Get image data
		data = CCDData.read(os.path.join(config['out_dir'], 'midproc',  im['filename']), unit=config['data_units'], hdu=config['ext'])
		proc_name = os.path.join(out_dir, im['filename'].split('.fits')[0] + '_' + im['object'] + '_')
		sci_date = data.meta['DATE_OBS'].split('T')[0]
		primary_hdu = fits.open(os.path.join(config['out_dir'], 'midproc', im['filename']),ignore_blank=True)[0]
		image_hdu = fits.ImageHDU(data=data,header=data.header)
		header = image.combine_headers(primary_hdu, image_hdu)
		#fits.HDUList([primary_hdu,image_hdu]).writeto(proc_name + '.fits')
		fits.writeto(proc_name + '.fits', data, header=header)

		# Subtract overscan, trim
		if config['overscan']:
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['overscan_region'])
			if use_variable_legendre:
				new_model = image.find_best_overscan_legendre_model(data[xmin:xmax,ymin:ymax], im, legendre_options, out_dir=os.path.join(config['out_dir'], 'plots'), **variable_legendre_options)
				overscan_options['model'] = new_model
			data = ccdproc.subtract_overscan(data, overscan=data[xmin:xmax,ymin:ymax], **overscan_options)
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data = ccdproc.trim_image(data[xmin:xmax,ymin:ymax])
			data.meta[config['data_region']] = f'[1:{data.shape[1]},1:{data.shape[0]}]'
			proc_name += 'OT'
			image_hdu = fits.ImageHDU(data=data,header=data.header)
			header = image.combine_headers(primary_hdu, image_hdu)
			#fits.HDUList([primary_hdu,image_hdu]).writeto(proc_name + '.fits')
			fits.writeto(proc_name + '.fits', data, header=header)

		# Subtract zero frame
		if config['zero']:
			data_temp = data
			# Find best master zero frame
			zero = image.find_best_masterframe('zero', im, config)
			data = CCDData.subtract(data_temp, zero)
			data.meta = data_temp.meta
			proc_name += 'Z'
			image_hdu = fits.ImageHDU(data=data,header=data.header)
			header = image.combine_headers(primary_hdu, image_hdu)
			#fits.HDUList([primary_hdu,image_hdu]).writeto(proc_name + '.fits')
			fits.writeto(proc_name + '.fits', data, header=header)

		# Subtract dark frame
		if config['dark']:
			data_temp = data
			# Find best master dark frame
			dark = image.find_best_masterframe('dark', im, config)
			data = CCDData.subtract(data_temp, dark)
			data.meta = data_temp.meta
			proc_name += 'D'
			image_hdu = fits.ImageHDU(data=data,header=data.header)
			header = image.combine_headers(primary_hdu, image_hdu)
			#fits.HDUList([primary_hdu,image_hdu]).writeto(proc_name + '.fits')
			fits.writeto(proc_name + '.fits', data, header=header)

		# Find the best master flat frame, divide
		if config['flat']:
			data_temp = data
			flat = image.find_best_masterframe('flat', im, config, date = sci_date)
			flat.data[np.where(flat.data==0)]='Nan'
			data = CCDData.divide(data_temp, flat)
			data.meta = data_temp.meta
			proc_name += 'F'
			image_hdu = fits.ImageHDU(data=data,header=data.header)
			header = image.combine_headers(primary_hdu, image_hdu)
			#fits.HDUList([primary_hdu,image_hdu]).writeto(proc_name + '.fits')
			fits.writeto(proc_name + '.fits', data, header=header)

		# Save image
		final_name = os.path.join(config['out_dir'], im['filename'].split('.fits')[0] + '_' + im['object'] + '.proc.fits')
		image_hdu = fits.ImageHDU(data=data,header=data.header)
		header = image.combine_headers(primary_hdu, image_hdu)
		#fits.HDUList([primary_hdu,image_hdu]).writeto(final_name)
		fits.writeto(final_name, data, header=header)

	return

def astrometry(config, file_info):
	# Get object files
	obj_ims = file_info[np.where(file_info['imagetyp']==config['object_image_keyword'])]
	scamp_ref = config['scamp_ref_cat']
	out_dir = config['out_dir']
	astrometry_out = os.path.join(out_dir,'astrometry')

	# For each image...
	for im in obj_ims:
		base_name = im['filename'].split('.fits')[0] + '_' + im['object']
		proc_name = os.path.join(out_dir, base_name + '.proc.fits')
		scamp_cat = os.path.join(astrometry_out, base_name + '_scamp.cat')
		scamp_config = os.path.join(astrometry_out, base_name + '_scamp.config')
		xml_out = os.path.join(astrometry_out, f'scamp.{base_name}.xml')
		# Run through astrometry.net
		os.system(f'solve-field -i {scamp_cat} -n {scamp_config} -D {astrometry_out} -o {base_name} --scale-units arcsecperpix --scale-low 0.2 --scale-high 0.3 -t 3 {proc_name}') #--ra {} --dec {} --radius {}

		# Run through SCAMP
		os.system(f'scamp {scamp_cat} -c {scamp_config} -AHEADER_SUFFIX .wcs -ASTREF_CATALOG {scamp_ref} -XML_NAME {xml_out} -CHECKPLOT_DEV PDF')

	return

# Stacking
def stack(config, file_info):

	# Get object files
	obj_ims = file_info[np.where(file_info['imagetyp']==config['object_image_keyword'])]
	out_dir = config['out_dir']

	# Find images to combine
	swarp_config = os.path.join(project_dir,config['swarp_config'])
	objects = np.unique(obj_ims['object'])
	filts = np.unique(obj_ims['filter'])

	# Run through SWarp
	for obj in objects:
		get_objs = obj_ims[np.where(obj_ims['object']==obj)]

		for filt in filts:
			get_filts = get_objs['filename'][np.where(get_objs['filter']==filt)]
			imgs_to_combine = []
			for file in get_filts:
				imgs_to_combine.append(os.path.join(out_dir, 'astrometry', file.split('.fits')[0] + f'_{obj}.new'))
			#print(imgs_to_combine)
			img_list_filename = os.path.join(out_dir, f'{obj}_{filt}.dat')
			with open(img_list_filename, 'w') as f:
			    for img_name in imgs_to_combine:
			        f.write(f'{img_name}\n')
			out_name = os.path.join(out_dir,f'{obj}_{filt}.fits')
			out_name_weight = os.path.join(out_dir,f'{obj}_{filt}.weight.fits')
			out_name_xml = os.path.join(out_dir,f'swarp.{obj}_{filt}.xml')
			os.system(f'swarp @{img_list_filename} -c {swarp_config} -IMAGEOUT_NAME {out_name} -WEIGHTOUT_NAME {out_name_weight} -WRITE_FILEINFO Y -XML_NAME {out_name_xml} -HEADER_SUFFIX _scamp.head -VERBOSE_TYPE FULL')

	return


def find_nearest_multiple(num, multiple_guess):

	multiple_up = multiple_guess
	multiple_down = multiple_guess

	while num % multiple_up != 0:
	    multiple_up+=1

	while num%multiple_down !=0:
	    multiple_down-=1

	if abs(multiple_guess-multiple_up) < abs(multiple_guess-multiple_down): nearest_multiple = multiple_up
	else: nearest_multiple = multiple_down

	return nearest_multiple
