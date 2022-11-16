import lbcred, os

# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Model galaxy using imfit. Measure SBF distance.')
    default_config = os.path.join(lbcred.project_dir, 'modeling-config_imfit.yml')
    parser.add_argument('-i','--image_dir', type=str, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
    parser.add_argument('--config', type=str, default=default_config, help='path of the .yml config file ')
    parser.add_argument('--run_imfit', type=bool, help='do imfit modelling')
    parser.add_argument('--run_artpop', type=bool, help='do artpop modelling')
    parser.add_argument('--overwrite', type=bool, help='overwrite existing directory for output images')
    parser.add_argument('--out_dir', type=str, help='specify output directory name/path; default is located in the same directory as image_dir with name \'lbcreduce_modeling_<date>_<time>\'')
    parser.add_argument('-ext','--extension', type=int, help='extension of FITS file data to be used')
    args = parser.parse_args()
    '''
    options = {
        'image_dir' : args.image_dir,
        'ext' : args.extension,
        'run_imfit' : args.run_imfit,
        'run_artpop' : args.run_artpop,
        'overwrite' : args.overwrite,
        'out_dir' : args.out_dir
    }
    '''
    options = {}
    lbcred.modeling_imfit(config_filename=args.config, options=options)
