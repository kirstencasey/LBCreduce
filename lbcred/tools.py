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
from astropy.io import fits
from astropy.table import Column


def initialize_config(input_options, config_filename):
	# Open and read config file
	with open(config_filename, 'r') as filename:
		config = yaml.load(filename, Loader=yaml.FullLoader)

	# Replace options input via command line into config
	for key in input_options:
		if input_options[key] != None:
			config[key] = input_options[key]

	return config

def textfile(config, end_notes, dir_overwritten):
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
	out_dir = config['out_dir']
	image_dir = config['image_dir']
	file = open(out_dir+'image-production-details.txt','w')
	datetime = time.strftime('%Y-%m-%d at %H:%M:%S',time.gmtime())
	lines = [f'This file describes the options used to process the images in {image_dir}. \nThe image processing was performed by lbcreduce on {datetime}.',
			'\n\nNotes input by user at runtime:\n', config['notes'], '\n\nNotes input by user after all image processing steps:\n', end_notes]

	# Make sure inputs are what were actually used if the user changed anything at all while running the program interactively
	file.writelines(lines)

	print('\nImage processing complete! :)\n')
	return

def bias(bias_type, config, file_info):

	raw_bias_info = file_info[np.where(file_info['imagetyp']==config['bias_image_keyword'])]
	processed_bias_info = raw_bias_info.copy()
	out_names = []

	# Make 2D bias image
	if bias_type == '2Dbias':
		# Loop through bias images
		for bias_im in raw_bias_info:
			# Get data
			data = CCDData.read(config['out_dir'] + 'midproc/' +  bias_im['filename'], unit=config['data_units'], hdu=config['ext'])

			# Trim overscan region
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data_trimmed = ccdproc.trim_image(data[xmin:xmax,ymin:ymax], **config['trim_options'])

			# Save trimmed image
			out_names.append(config['out_dir'] + 'midproc/' + bias_im['filename'].replace('.fits','_T.fits'))
			primary_hdu = fits.open(config['out_dir'] + 'midproc/' +  bias_im['filename'])[0]
			image_hdu = fits.ImageHDU(data=data_trimmed,header=data_trimmed.header)
			fits.HDUList([primary_hdu,image_hdu]).writeto(out_names[-1])


	# Make zero-frame bias image
	if bias_type == 'zero':
		# Loop through bias images
		for bias_im in raw_bias_info:
			# Get data
			data = CCDData.read(config['out_dir'] + 'midproc/' +  bias_im['filename'], unit=config['data_units'], hdu=config['ext'])

			# Subtract overscan
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['overscan_region'])
			data_o = ccdproc.subtract_overscan(data, overscan=data[xmin:xmax,ymin:ymax], **config['overscan_options'])

			# Trim overscan region
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data_ot = ccdproc.trim_image(data_o[xmin:xmax,ymin:ymax], **config['trim_options'])

			# Save overscan subtracted, trimmed image
			out_names.append(config['out_dir'] + 'midproc/' + bias_im['filename'].replace('.fits','_OT.fits'))
			primary_hdu = fits.open(config['out_dir'] + 'midproc/' +   bias_im['filename'])[0]
			image_hdu = fits.ImageHDU(data=data_ot,header=data_ot.header)
			fits.HDUList([primary_hdu,image_hdu]).writeto(out_names[-1])

	processed_bias_info['filename'] = Column(data=out_names, name='filename')

	# Combine images to make master bias
	for chip in ['-chip1','-chip2','-chip3','-chip4']:
		master_name = config['out_dir'] + 'midproc/masterbias_' + bias_type + chip

		# Get correct chip
		mask = [idx for idx,fi in enumerate(processed_bias_info['filename']) if chip in fi]
		chip_info = processed_bias_info[mask]
		chip_info_R = chip_info[np.where(chip_info['instrument']==config['lbc_red'])]
		chip_info_B = chip_info[np.where(chip_info['instrument']==config['lbc_blue'])]

		# Sort by instrument
		masterbias_R = ccdproc.combine(chip_info_R['filename'], output_file=master_name+'_R.fits', unit=config['data_units'], **config['combine_options'])
		masterbias_B = ccdproc.combine(chip_info_B['filename'], output_file=master_name+'_B.fits', unit=config['data_units'], **config['combine_options'])

	# Get feedback on master bias

	return

# Overscan and trim
def model_overscan(options, file_info):
	'''

	'''
	# Get biases (image.get_images)

	# Loop through files:

		# Subtract overscan

		# Trim image

		# Save the result
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

	# Loop through files to calibrate flats (subtract bias, trim overscan):
	for flat_im in processed_flats_info:
		# Get data
		data = CCDData.read(config['out_dir'] + 'midproc/' +  flat_im['filename'], unit=config['data_units'], hdu=config['ext'])
		dates.append(data.meta['DATE_OBS'].split('T')[0])

		# Trim overscan region
		xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
		data_t = ccdproc.trim_image(data[xmin:xmax,ymin:ymax], **config['trim_options'])

		# Subtract 2D bias image
		if flat_im['instrument'] == config['lbc_red']:
			inst_color = '_R'
		elif flat_im['instrument'] == config['lbc_blue']:
			inst_color = '_B'
		chip = '-' + flat_im['filename'].split('-')[-1].split('.fits')[0]

		bias_name = config['out_dir'] + 'midproc/' +  'masterbias_2Dbias' + chip + inst_color + '.fits'
		bias2D = CCDData.read(bias_name, unit=config['data_units'])
		data_ot = CCDData.subtract(data_t, bias2D)
		data_ot.meta = data_t.meta

		# Save calibrated flat
		out_names.append(config['out_dir'] + 'midproc/' + flat_im['filename'].replace('.fits','_OT.fits'))
		primary_hdu = fits.open(config['out_dir'] + 'midproc/' +   flat_im['filename'])[0]
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
			for chip in ['-chip1','-chip2','-chip3','-chip4']:
				master_name = config['out_dir'] + 'midproc/masterflat' + '_' + day + '_' + chip + '_' + filt + '.fits'
				mask = [idx for idx,fi in enumerate(date_info) if chip in fi]
				chip_info = date_info[mask]
				if len(chip_info) == 0: continue
				if len(chip_info) == 1:
					shutil.copyfile(chip_info[0], master_name)
					warnings.warn(f'Only one flat to \'combine\' for {master_name}.', AstropyUserWarning)
					continue
				masterflat = ccdproc.combine(chip_info, output_file=master_name, unit=config['data_units'], **config['combine_options'])

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
		out_dir += 'midproc/'

	# Loop through images
	for im in sci_ims:
		# Get image data
		data = CCDData.read(config['out_dir'] + 'midproc/' +  im['filename'], unit=config['data_units'], hdu=config['ext'])
		proc_name = out_dir + im['filename'].split('.fits')[0] + '_' + im['object'] + '_'

		# Subtract overscan, trim
		if config['overscan']:
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['overscan_region'])
			data = ccdproc.subtract_overscan(data, overscan=data[xmin:xmax,ymin:ymax], **config['overscan_options'])
			xmin, xmax, ymin, ymax = image.get_ccd_section(data.meta, config['science_region'])
			data = ccdproc.trim_image(data[xmin:xmax,ymin:ymax], **config['trim_options'])
			proc_name += 'OT'

		# Subtract zero frame
		if config['zero']:
			data_temp = data
			# Find best master zero frame
			zero = image.find_best_masterframe('zero', im, config)
			data = CCDData.subtract(data_temp, zero)
			data.meta = data_temp.meta
			proc_name += 'Z'

		# Subtract dark frame
		if config['dark']:
			data_temp = data
			# Find best master dark frame
			dark = image.find_best_masterframe('dark', im, config)
			data = CCDData.subtract(data_temp, dark)
			data.meta = data_temp.meta
			proc_name += 'D'

		# Find the best master flat frame, divide
		if config['flat']:
			data_temp = data
			flat = image.find_best_masterframe('flat', im, config)
			flat.data[np.where(flat.data==0)]='Nan'
			data = CCDData.divide(data_temp, flat)
			data.meta = data_temp.meta
			proc_name += 'F'

		# Save image
		fits.writeto(proc_name + '.fits', data.data, fits.Header(data.meta))

	return

# Stacking
def stack(config, file_info):
	'''
	Note: talk to Chris and Johnny about this
	'''

	# Run through astrometry.net

	# Run through SCAMP (first SDSS/Pan-STARRS catalog then Gaia)

	# Re-combine chips

	# Run through SWarp

	return
