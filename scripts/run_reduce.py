import os
import lbcred

# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Process raw LBC images.')
    default_config = os.path.join(lbcred.project_dir, 'lbcreduce-config.yml')
    parser.add_argument('-i','--image_dir', type=str, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
    parser.add_argument('--config', type=str, default=default_config, help='path of the .yml config file ')
    parser.add_argument('-ext','--extension', type=int, help='extension of FITS file data to be used')
    parser.add_argument('-o','--overscan', type=bool, help='model, subtract, and trim overscan')
    parser.add_argument('-z','--zero', type=bool, help='include zero frame subtraction in reduction')
    parser.add_argument('-d','--dark', type=bool, help='include dark subtraction in reduction')
    parser.add_argument('-f', '--flat', type=bool, help='include flat-fielding step in reduction')
    parser.add_argument('-s', '--stack', type=bool, help='stack images for each OBJECT after reduction steps')
    parser.add_argument('--reduce_objects', type=bool, help='include reduction of object images in image processing; only set to False if you\'re only interested in producing master flats/bias images, etc. or stacking previously processed images')
    parser.add_argument('--filenames', type=list, help='list of files to include in image reduction; should include any object images, flats, darks, bias images, etc. necessary for reduction')
    parser.add_argument('--object', type=str, help='OBJECT keyword in header; used when image reduction of only one object is desired')
    parser.add_argument('--propid', type=str, help='PROPID keyword in header; used when image reduction of only one propid is desired')
    parser.add_argument('--glob_include', type=str, help='Unix-style filename segment used to include files; only necessary if you don\'t want all files in image_dir used (eg. \'lbcb*\'); applied after include_filenames')
    parser.add_argument('--exclude', type=str, help='Unix-style filename segment used to exclude files; (eg. \'*20191220*\')')
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
        'propid' : args.propid,
        'glob_include' : args.glob_include,
        'exclude' : args.exclude,
        'overwrite' : args.overwrite,
        'out_dir' : args.out_dir,
        'notes' : args.notes,
        'check_output' : args.check_output
    }


    lbcred.reduce(config_filename=args.config, options=options)
