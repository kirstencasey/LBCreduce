import os, lbcred


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    default_config = os.path.join(lbcred.project_dir, 'register_calibration-config.yml')
    parser.add_argument('-t', '--tmp-path', default='/tmp', help='temporary path for intermediate steps')
    parser.add_argument('--register-images', default=True)
    parser.add_argument('-c','--calibrate-images', default=True)
    parser.add_argument('-r', '--ref-cat', help='reference catalog file name')
    parser.add_argument('--index-path', default=None, help='path to astrometry.net index files')
    parser.add_argument('--make-plots', action='store_true')
    parser.add_argument('--log-level', default='info', help='log level (debug, info, warn, error, critical)')
    parser.add_argument('--log-fn', default=None, help='log file name')
    parser.add_argument('--config_fn', type=str, default=default_config, help='path of the .yml config file ')
    args = parser.parse_args()
    lbcred.utils.setup_logger(args.log_level, args.log_fn)

    if args.register_images:
        # Register r-band images
        config, sky_pos, src = lbcred.register_images(tmp_path=args.tmp_path, bandpass='R', ref_cat=args.ref_cat, index_path=args.index_path, make_plots=args.make_plots, config_fn=args.config_fn)

        # Register b-band images
        config, sky_pos, src = lbcred.register_images(tmp_path=args.tmp_path, bandpass='B', ref_cat=args.ref_cat, index_path=args.index_path, make_plots=args.make_plots, config_fn=args.config_fn)

        config['image_dir'] = config['out_dir']
        if config['subtract_background']: config['glob_select'] = 'backsub_' + config['glob_select'].replace('.fits','_reg.fits')
        else: config['glob_select'] = config['glob_select'].replace('.fits','_reg.fits')

    # Calibrate images
    if args.calibrate_images:
        if args.register_images: config = config
        else : config = args.config_fn

        lbcred.calibrate_images(config=config)
