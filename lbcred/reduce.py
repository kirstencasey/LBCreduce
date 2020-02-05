'''
Functions for reducing LBC images:
	- Overscan and trim
	- Constructing master biases, flat fields, etc.
	- Stacking images
'''
import numpy as np
from lbcred import image
from astropy.stats import sigma_clip, mad_std
import ccdproc


# Overscan and trim
def overscan():
	'''
	
	'''
	# Get biases (image.get_images)

	# Loop through files: 

		# Subtract overscan

		# Trim image

		# Save the result

# Flat fielding
def flat():
	'''
	'''
	# Get flats (image.get_images)

	# Loop through files: 
	
		# 

# Combine calibrated images
def combine():
	'''
	'''
	# Combine calibrated bias images

	# Or, combine calibrated flats

		# Rescale individual flats

		# Combine flats

# Calibrate dark frames
def dark():
	'''
	'''

# Stacking
def stack():
	'''
	Note: talk to Chris and Johnny about this 
	'''
