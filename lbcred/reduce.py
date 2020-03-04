import tools, interactive, image
import numpy as np
import os, time, shutil, warnings

affirmative = ['y','Y','yes','Yes']
negative = ['n','N','no','No']
keys = ['imagetyp']

def reduce(options, config_filename):

    # Read configuration file - NOTE THAT THIS FUNC ASSUMES DEFAULT OPTIONS FROM COMMAND LINE ARE None
    initial_config = tools.initialize_config(options, config_filename) # Overwrite anything in the config file that was supplied via the command line

    # Do directory stuff
    config, dir_overwritten = interactive.initialize_directories(initial_config)

    # Get raw images
    raw_files, config = image.get_images(config)
    #raw_files, config= image.check_files(raw_files, config)		###### Doesn't work when nothing in config is changed

    # Create master bias images (2D bias and zero frame)
    if config['zero'] or config['flat']:
        masterbias_2D, config = tools.bias('2Dbias', config, raw_files) # Needed for flat fields
        masterbias_zero, config = tools.bias('zero', config, raw_files)
    '''
    # Calibrate dark frames
    if config['dark']:
        options_now = tools.dark(options_now, raw_files)

    # Check counts, calibrate flat fields
    if config['flat']:
        options_now = tools.flat(options_now, raw_files)

    # Process images
    if config['reduce_objects']:
    options_now = tools.process(options_now, raw_files)

    # Stack images
    if config['stack']:
        options_now = tools.stack(options_now, raw_files)

    # Ask for final notes (if check_output is not 0)
    if config['check_output']:
        response = input('Are there any other notes you would like to add to the output textfile? [y/n]: ')
        if response in affirmative:
            end_notes = input('Notes to add: ')
        else:
            end_notes = ''
    '''
    # Make textfile specifying options used
    if config != initial_config:
        final_config = config
    else:
        final_config = None

    # Make textfile
    tools.textfile(initial_config, end_notes, dir_overwritten, final_config)

    # End
    return


# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Process raw LBC images.')
    parser.add_argument('-i','--image_dir', type=str, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
    parser.add_argument('--config', type=str, default='./lbcreduce-config.yml', help='path of the .yml config file ')
    parser.add_argument('-ext','--extension', type=int, help='extension of FITS file data to be used')
    parser.add_argument('-o','--overscan', type=bool, help='model, subtract, and trim overscan')
    parser.add_argument('-z','--zero', type=bool, help='include zero frame subtraction in reduction')
    parser.add_argument('-d','--dark', type=bool, help='include dark subtraction in reduction')
    parser.add_argument('-f', '--flat', type=bool, help='include flat-fielding step in reduction')
    parser.add_argument('-s', '--stack', type=bool, help='stack images for each OBJECT after reduction steps')
    parser.add_argument('--reduce_objects', type=bool, help='include reduction of object images in image processing; only set to False if you\'re only interested in producing master flats/bias images, etc. or stacking previously processed images')
    parser.add_argument('--filenames', type=list, help='list of files to include in image reduction; should include any object images, flats, darks, bias images, etc. necessary for reduction')
    parser.add_argument('--object', type=str, help='OBJECT keyword in header; used when image reduction of only one object is desired')
    parser.add_argument('--glob_include', type=str, help='Unix-style filename segment used to include files; only necessary if you don\'t want all files in image_dir used (eg. \'lbcb*\'); applied after include_filenames')
    parser.add_argument('--glob_exclude', type=str, help='Unix-style filename segment used to exclude files; (eg. \'*20191220*\')')
    parser.add_argument('--overwrite', type=bool, help='overwrite existing directory for output images')
    parser.add_argument('--out_dir', type=str, help='specify output directory name/path; default is located in the same directory as image_dir with name \'lbcreduce_<date>_<time>\'')
    parser.add_argument('-n','--notes', type=str, help='notes added to text file included in output dir along with other options used')
    parser.add_argument('-c','--check_output', type=bool, help='check output of each reduction step before continuing')
    args = parser.parse_args()

    options = {
        'image_dir' : args.image_dir,
        'ext' : args.extension,
        'overscan' : args.overscan,
        'zero' : args.zero,
        'dark' : args.dark,
        'flat' : args.flat,
        'stack' : args.stack,
        'reduce_objects' : args.reduce_objects,
        'filenames' : args.filenames,
        'object' : args.object,
        'glob_include' : args.glob_include,
        'glob_exclude' : args.glob_exclude,
        'overwrite' : args.overwrite,
        'out_dir' : args.out_dir,
        'notes' : args.notes,
        'check_output' : args.check_output
    }


    reduce(options, config_filename=args.config)
