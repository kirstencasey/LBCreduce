import lbcred, os

# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Model galaxy using imfit. Measure SBF distance.')
    default_imfit_config = os.path.join(lbcred.project_dir, 'modeling-config_imfit.yml')
    default_sbf_config = os.path.join(lbcred.project_dir, 'modeling-config_sbf.yml')
    parser.add_argument('-i','--image_dir', type=str, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
    parser.add_argument('--out_dir', type=str, help='specify output directory name/path; default is located in the same directory as image_dir with name \'lbcreduce_modeling_<date>_<time>\'')
    parser.add_argument('--imfit_config', type=str, default=default_imfit_config, help='path of the .yml config file ')
    parser.add_argument('--sbf_config', type=str, default=default_sbf_config, help='path of the .yml config file ')
    parser.add_argument('--run_imfit', type=bool, default=True, help='do imfit modelling')
    parser.add_argument('--run_artpop', type=bool, help='do artpop modelling before running imfit')
    parser.add_argument('--run_sbf', type=bool, default=True, help='do sbf measurement')
    args = parser.parse_args()

    options = {
        'image_dir' : args.image_dir,
        'out_dir' : args.out_dir,
        'run_imfit' : args.run_imfit,
        'run_artpop' : args.run_artpop
    }

    if args.run_imfit:
        bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, functions, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true = lbcred.modeling_imfit(config_filename=args.imfit_config, options=options)

    options.pop('run_imfit')
    options.pop('run_artpop')
    options['color'] = color
    options['model_fn'] = model1_fn
    options['resid_fn'] = resid1_fn
    options['model_summary_fn'] = bestfit1_fn

    if args.run_sbf:
        sbf_mag, dist_a, dist_b = lbcred.modeling_sbf(config_filename=args.sbf_config, options=options, imfit_functions=functions)
