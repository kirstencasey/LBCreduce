import lbcred, os

# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Model galaxy using imfit. Measure SBF distance.')
    default_config = os.path.join(lbcred.project_dir, 'modeling-config_sbf.yml')
    parser.add_argument('-i','--image_dir', type=str, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
    parser.add_argument('--model_fn', type=str, help='filename of the galaxy model located in image_dir')
    parser.add_argument('--resid_fn', type=str, help='filename of the galaxy residual located in image_dir')
    parser.add_argument('--model_summary_fn', type=str, help='filename of sersic fitting result from imfit - should contain info on radius, ellipticity, position, and position angle')
    parser.add_argument('--radius', type=float, help='major axis of model in pixels')
    parser.add_argument('--xpos', type=float, help='x-position (FITS standard) of model')
    parser.add_argument('--ypos', type=float, help='x-position (FITS standard) of model')
    parser.add_argument('--ellip', type=float, help='ellipticity of model')
    parser.add_argument('--pa', type=float, help='position angle of model in degrees')
    parser.add_argument('-n','--num_iters', type=int, help='number of iterations for sbf analysis')
    parser.add_argument('--text_fn', type=str, help='name of text file to write SBF measurement results (saved in out_dir)')
    parser.add_argument('--config', type=str, default=default_config, help='path of the .yml config file ')
    parser.add_argument('--overwrite', type=bool, help='overwrite existing directory for output images')
    parser.add_argument('--out_dir', type=str, help='specify output directory name/path; default is located in the same directory as image_dir with name \'lbcreduce_modeling_<date>_<time>\'')
    parser.add_argument('-ext','--extension', type=int, help='extension of FITS file data to be used')
    args = parser.parse_args()

    options = {
        'image_dir' : args.image_dir,
        'model_fn' : args.model_fn,
        'resid_fn' : args.resid_fn,
        'model_summary_fn' : args.model_summary_fn,
        'radius' : args.radius,
        'xpos' : args.xpos,
        'ypos' : args.ypos,
        'ellip' : args.ellip,
        'pa' : args.pa,
        'num_iters' : args.num_iters,
        'text_fn' : args.text_fn,
        'ext' : args.extension,
        'overwrite' : args.overwrite,
        'out_dir' : args.out_dir
    }

    lbcred.modeling_sbf(config_filename=args.config, options=options)
