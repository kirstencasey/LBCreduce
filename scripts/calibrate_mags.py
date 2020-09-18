import glob, os, random
from astropy.table import vstack
from lbcred import detection, utils



def calibrate_mags(data_path,ref_cat,run_label,txt_filename,flux_apers,make_plots):

    # Check for necessary files, etc.
    files = glob.glob(os.path.join(data_path,'lbc*.fits'))
    num_apers = len(utils.list_of_strings(flux_apers))
    flux_params = f'FLUX_AUTO,FLUX_APER({num_apers}),ALPHA_J2000,DELTA_J2000'
    catalogs = []
    ##########################################
    rand_file = random.randint(0,len(files)-2)
    files = files[rand_file:rand_file+2]
    ##########################################
    # Run SE on each image
    for file in files:
        fn_base = file.split('.fits')[0]
        cat_fn = f'{fn_base}.cat'
        cat = detection.sextractor.run(file,DETECT_MINAREA=3,DETECT_THRESH=10,PIXEL_SCALE=0.225,PHOT_APERTURES=flux_apers,
                             catalog_path=cat_fn,extra_params=flux_params)

        '''
        star_query = 'FLAGS==0 and ISOAREA_IMAGE > 5 and \
                      FWHM_IMAGE > 1 and FWHM_IMAGE < 26'
        cat = cat[cat.to_pandas().query(star_query).index.values]
        '''
        utils.sextractor_cat_to_ds9reg(cat, outfile=f'{fn_base}.reg')

        cat['filename'] = file

        if len(catalogs) == 0:
            catalogs = cat
        else:
            catalogs = vstack([catalogs,cat])

    # For each object, plot flux given different aperature types: ex. FLUXAPER(5), FLUXAPER(10), etc.
    print(catalogs)
    ### Get random sample of objects in catalog (save N objects)
    #for obj in catalogs:

        # Plot the flux measured given the size of the aperature (aperature size on the x-axis, flux on the y-axis)

        # Plot a horiz. line where the "auto" aperature is


    # Decide on a good flux aperature (also decide how to exclude non-stellar objects) - Look at the plots to see how to do this

    #------------------------------------------------------------------

    # Calculate instrumental magnitude of object using measured flux: m_inst,0 = -2.5 log10(flux/t_exp)

    # Calculate correction to m_inst due to extinction,

    # Calculate color term and color index for each star

    # Cross-match with Pan-STARRS catalog - find m_calib for given star

    # Calculate zero-point for given star - make density map of zero points for each filter

    # Decide on best zero-point for each filter (consider average/median/etc. zero-point from many images/starts)

    return



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-d', '--data-path', required=True,
                        help='path to data, which have had basic reductions')
    parser.add_argument('-r', '--ref-cat', help='reference catalog file name')
    parser.add_argument('-l', '--run-label', help='a unique label for sextractor temporary files')
    parser.add_argument('-f', '--sum-file', default='calibration_summary.txt',
                        help='filename for file summarizing magnitude calibration')
    parser.add_argument('-a', '--flux-aper', default='2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40',
                        help='SExtractor FLUXAPER(#)s to compare (in pixels)')
    parser.add_argument('--pixscale', default=0.226, type=float)
    parser.add_argument('--make-plots', action='store_true')
    parser.add_argument('--log-level', default='info',
                        help='log level (debug, info, warn, error, critical)')
    parser.add_argument('--log-fn', default=None, help='log file name')
    args = parser.parse_args()
    utils.setup_logger(args.log_level, args.log_fn)

    calibrate_mags(
        data_path=args.data_path,
        ref_cat=args.ref_cat,
        run_label=args.run_label,
        txt_filename = args.sum_file,
        flux_apers = args.flux_aper,
        make_plots=args.make_plots
    )
