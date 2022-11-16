from lbcred import tools, interactive, image
import numpy as np
import os, time, shutil, warnings
from lbcred.log import logger
from shutil import copyfile


affirmative = ['y','Y','yes','Yes']
negative = ['n','N','no','No']
keys = ['imagetyp']
#warnings.simplefilter('ignore', FITSFixedWarning)

def reduce(config_filename, options = {}):

    # Read configuration file - NOTE THAT THIS FUNC ASSUMES DEFAULT OPTIONS FROM COMMAND LINE ARE None
    logger.info('Getting ready...')
    initial_config = tools.initialize_config(config_filename, options) # Overwrite anything in the config file that was supplied via the command line
    tools.setup_logger(initial_config['logger_level'], log_fn=initial_config['log_to_file'])
    # Save copy of config to output directory
    copyfile(config_filename,os.path.join(initial_config['out_dir'],config_filename.split('/')[-1]))

    # Do directory stuff
    config, dir_overwritten = interactive.initialize_directories(initial_config)

    # Get raw images
    logger.info('Gathering image information...')
    image_info = image.get_image_info(config)
    image.check_files(image_info, config)
    image_info = image.separate_chips(image_info, config)

    # Create master bias images (2D bias and zero frame)
    if config['zero'] or config['flat']:
        logger.info('Calibrating bias images...')
        tools.bias('2Dbias', config, image_info) # Needed for flat fields
        tools.bias('zero', config, image_info)   #### MAKE THESE MORE EFFICIENT (COMBINE THEM WHERE POSSIBLE) #####

    # Calibrate dark frames
    if config['dark']:
        logger.info('Calibrating dark frames...')
        tools.dark(config, image_info)

    # Check counts, calibrate flat fields
    if config['flat']:
        logger.info('Calibrating flats...')
        tools.flat(config, image_info)

    # Process images
    if config['reduce_objects']:
        logger.info('Processing science images...')
        tools.process(config, image_info)

    if config['astrometry']:
        logger.info('Finding astrometric solutions...')
        tools.astrometry(config, image_info)

    # Stack images
    if config['stack']:
        logger.info('Stacking science images...')
        tools.stack(config, image_info)

    # Ask for final notes (if check_output is not 0)
    end_notes = ''
    if config['check_output']:
        # TEMP:
        '''
        response = input('Are there any other notes you would like to add to the output textfile? [y/n]: ')
        if response in affirmative:
            end_notes = input('Notes to add: ')
        '''
    # Make textfile
    tools.textfile(config, end_notes, dir_overwritten)

    # End
    return
