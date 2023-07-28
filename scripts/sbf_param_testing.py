import lbcred, os, yaml, random, glob
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def sbf_test(config_filename, options):
    
    # Create relevant directories, move relevant files if necessary
    
    
    # Run SBF calculation
    sbf_mag, dist_a, dist_b = lbcred.modeling_sbf(config_filename=sb_config, options=option, imfit_functions=functions,run_id=run_id)
    

    return


# Run 'sbf_test' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Model galaxy using imfit. Measure SBF distance.')
    default_config = os.path.join(lbcred.project_dir, 'modeling-config_sbf.yml')
    parser.add_argument('-i','--init_dir', type=str, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
   parser.add_argument('--radius', type=float, help='major axis of model in pixels')
    parser.add_argument('--ellip', type=float, help='ellipticity of model')
    parser.add_argument('--pa', type=float, help='position angle of model in degrees')
    parser.add_argument('--text_fn', type=str, help='name of text file to write SBF measurement results (saved in out_dir)')
    parser.add_argument('--config', type=str, default=default_config, help='path of the .yml config file ')
    parser.add_argument('--out_dir', type=str, help='specify output directory name/path; default is located in the same directory as image_dir with name \'lbcreduce_modeling_<date>_<time>\'')
    args = parser.parse_args()

    options = {
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

    sbf_test(config_filename=args.config, options=options)
