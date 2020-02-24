import tools, interactive, image
import numpy as np
import os, time, shutil, warnings

def reduce(image_dir, do_overscan, overscan_model, combine_func, std_type, std, do_zero, do_dark, do_flat, do_stack, filenames, glob_include, glob_exclude, include_attributes, exclude_attributes, overwrite, out_dir, notes, check_output):

    initial_options = {
        'image_dir' : image_dir,
        'do_overscan' : do_overscan,
        'overscan_model' : overscan_model,
        'combine_func' : combine_func,
        'std_type' : std_type,
        'std' : std,
        'do_zero' : do_zero,
        'do_dark' : do_dark,
        'do_flat' : do_flat,
        'do_stack' : do_stack,
        'filenames' : filenames,
        'glob_include' : glob_include,
        'glob_exclude' : glob_exclude,
        'include_attributes' : include_attributes,
        'exclude_attributes' : exclude_attributes,
        'overwrite' : overwrite,
        'out_dir' : out_dir,
        'check_output' : check_output
    }

    options_now = initial_options

    # Do directory stuff
    options_now, dir_overwritten = tools.initialize_directories(options_now)

    # Get raw images
    all_files, options_now = image.get_images(options_now)

    # Create 2D bias images
    if do_flat:
        masterbias, options_now = tools.bias('2Dbias', options_now)

    # Create zero frame
    if do_zero:
        zeroframe, options_now = tools.bias('zero', options_now)

    # Calibrate dark frames
    if do_dark:
        options_now = tools.dark(options_now)

    # Check counts, calibrate flat fields
    if do_flat:
        options_now = tools.flat(options_now)

    # Process images
    options_now = tools.process(options_now)

    # Stack images
    if do_stack:
        options_now = tools.stack(options_now)

    # Ask for final notes (if check_output is not 0)
    if options_now['check_output'] != 0:
        response = input('Are there any other notes you would like to add to the output textfile? [y/n]: ')
        if response == 'y' or response == 'yes':
            end_notes = input('Notes to add: ')
        else:
            end_notes = ''

    # Make textfile specifying options used
    if initial_options != options_now:
        final_options = options_now
    else:
        final_options = None

    # Make textfile
    tools.textfile(initial_options, notes, end_notes, dir_overwritten, final_options)

    # End
    return


# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Process raw LBC images.')
    parser.add_argument('-i','--image_dir', type=str, required=True, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
    parser.add_argument('-o','--overscan', type=bool, default=True, help='model, subtract, and trim overscan')
    parser.add_argument('-m','--overscan_model', type=str, default='median', choices=['mean', 'median'], help='model fit to overscan before subtraction; mean includes sigma clipping') # Needs work
    #parser.add_argument('-l','--legendre_order', type=int, default=3, choices=[0,1,2,3,4,5,6], help='if Legendre function is requested, this is the order of the function to be used') # Needs work
    parser.add_argument('--combine_func', type=str, default='mean', choices=['mean','median'], help='function used to combine bias images, etc.; mean includes sigma clipping')
    parser.add_argument('--std_type', type=str, default='mad_std', choices=['mad_std','std'], help='function used to identify outliers for modeling overscan, choosing useful flats, etc.')
    parser.add_argument('--std', type=int, default=3, help='number of standard deviations from mean/median used to determine outliers for modeling overscan, choosing useful flats, etc.')
    parser.add_argument('-z','--zero', type=bool, default=True, help='include zero frame subtraction in reduction')
    parser.add_argument('-d','--dark', type=bool, default=False, help='include dark subtraction in reduction')
    parser.add_argument('-f', '--flat', type=bool, default=True, help='include flat-fielding step in reduction')
    parser.add_argument('-s', '--stack', type=str, default=True,help='stack images for each OBJECT after reduction steps')
    parser.add_argument('--include_filenames', type=str or list of str, default=None,help='list of files to include in image reduction; should include any object images, flats, darks, bias images, etc. necessary for reduction')
    parser.add_argument('--glob_include', type=list, default=None, help='Unix-style filename segment used to include files; only necessary if you don\'t want all files in image_dir used (eg. \'lbcb*\')')
    parser.add_argument('--glob_exclude', type=list, default=None, help='Unix-style filename segment used to exclude files; (eg. \'*20191220*\')')
    parser.add_argument('--include_attributes', type=dict, default=None, help='dictionary with FITS header keyword(s) to check for and a list of value(s) of that keyword to include in processing (ex. {\'propid\' : [\'OSU_dwarfsimg\']}); note that this does not affect the inclusion of flats, bias images, or dark frames')
    parser.add_argument('--exclude_attributes', type=dict, default=None, help='dictionary with FITS header keyword(s) to check for and a list of value(s) of that keyword to exclude in processing (ex. {\'propid\' : [\'OSU_dwarfsimg\']}); note that this does not work for flats, bias images, or dark frames and that this does nothing if \'include_attribute\' is specified')
    parser.add_argument('--overwrite', type=bool, default=False, help='overwrite existing directory for output images')
    parser.add_argument('--out_dir', type=str, default=None, help='specify output directory name/path; default is located in the same directory as image_dir with name \'lbcreduce_<date>_<time>\'')
    parser.add_argument('-n','--notes', type=str, default='', help='notes added to text file included in output dir along with other options used')
    parser.add_argument('-c','--check_output', type=int, default=1, choices=[0,1,2], help='check output of each reduction step interactively; 0 is no checking, 1 is some checking, 2 is all checking')
    args = parser.parse_args()

    reduce(image_dir=args.image_dir, do_overscan=args.overscan, overscan_model=args.overscan_model, combine_func=args.combine_func, std_type=args.std_type, std=args.std, do_zero=args.zero, do_dark=args.dark, do_flat=args.flat, do_stack=args.stack, filenames=args.include_filenames, glob_include=args.glob_include, glob_exclude=args.glob_exclude, include_attributes=args.include_attributes, exclude_attributes=args.exclude_attributes, overwrite=args.overwrite, out_dir=args.out_dir, notes=args.notes, check_output=args.check_output)
