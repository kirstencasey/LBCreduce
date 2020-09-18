import os
from glob import glob
import numpy as np
from shutil import copyfile
from astropy.io import fits
from astropy.table import Table
from astropy.modeling import models, fitting
from lbcred.detection import extract_bright_stars, sextractor_sky_model
from lbcred.astrometry import solve_field, check_astrometric_solution, TweakWCS
from lbcred import improc, utils, logger

# astropy fits warnings are so annoying
import warnings
from tqdm import tqdm
from astropy.wcs import FITSFixedWarning
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=FITSFixedWarning)
warnings.simplefilter('ignore', category=VerifyWarning)


def fetch_files(path, band='R', exptime_range=[0, 50],
                name_must_contain='M81blob.proc'):
    files = glob(os.path.join(path, f'lbc{band.lower()}*{name_must_contain}*'))
    files = [f for f in files if\
             float(fits.getheader(f)['EXPTIME']) > exptime_range[0] if\
             float(fits.getheader(f)['EXPTIME']) < exptime_range[1]]
    # sort files by time
    # NOTE: This assumes the file
    # names contain the exposure number
    files.sort()
    return files


@utils.func_timer
def register_images(data_path, out_path, center, tmp_path,
                    bandpass, subtract_back, back_model, stack, stack_type,
                    index_path=None, ref_cat=None, output_dim=[4000, 6000],
                    pixscale=0.226, make_plots=False):

    # fetch reference catalog if necessary
    if ref_cat is None:
        fn =  f'panstarrs-{center[0]:.1f}-{center[1]:.1f}.dat'
        fn = os.path.join(out_path, fn)
        # if has been fetched previously, load local copy
        if os.path.isfile(fn):
            logger.info('reading reference catalog from ' + fn)
            ref_cat = Table.read(fn)
        else:
            # otherwise fetch using pyvo
            from lbcred.tap import PanstarrsTAP
            panstarrs_tap = PanstarrsTAP()
            logger.warning('no reference catalog. fetching PanSTARRS catalog')
            query_template = panstarrs_tap.\
                             default_query_template.\
                             replace('10', '16').\
                             replace('23', '22')
            ref_cat = panstarrs_tap.quey_region(*center, radius=0.25,
                                                query_template=query_template)
            ref_cat.rename_column('raMean', 'ra')
            ref_cat.rename_column('decMean', 'dec')
            ref_cat.write(fn, overwrite=True)
    else:
        # if you passed a file name, load that
        ref_cat = Table.read(ref_cat)

    if subtract_back:
        back_out = os.path.join(out_path, 'back_subtracted')
        utils.mkdir_if_needed(back_out)
        all_files = glob(os.path.join(data_path, f'lbc{bandpass.lower()}*M81blob.proc.fits'))
        for fi in all_files:
            fi_base = fi.split('/')[-1].replace('M81blob.proc.fits','backsub_M81blob.proc.fits')
            copyfile(fi, os.path.join(back_out,fi_base))
        data_path = back_out

    # short exposures for finding the astrometric solution
    files_cali = fetch_files(data_path, bandpass, [0, 50])
    num_cali = len(files_cali)

    # long exposures for making science images
    files_sci = fetch_files(data_path, bandpass, [200, 500])
    num_sci = len(files_sci)

    logger.info(f'Solving astrometry for {bandpass}-band frames')

    # loop over file types (calibration and science)
    for ftype, files, num in zip(['CALI', 'SCI'], [files_cali, files_sci], [num_cali, num_sci]):
        logger.info(f'Solving astrometry for {bandpass}-band {ftype} frames')

        frame_out = os.path.join(out_path, ftype.lower())
        utils.mkdir_if_needed(frame_out)

        astrom = []
        logger.start_tqdm()

        # run solve-field on each image
        for fn in tqdm(files):

            if subtract_back:
                logger.info('Subtracting background for ' + fn)
                sky = sextractor_sky_model(fn)

                if back_model == 'plane':
                    p_init = models.Polynomial2D(degree=1)
                    fit_p = fitting.LinearLSQFitter()
                    y, x = np.mgrid[:sky.shape[0], :sky.shape[1]]
                    p = fit_p(p_init,x,y,sky)

                    with fits.open(fn, mode='update') as hdul:
                        hdul[0].data -= p(x,y)
                        hdul[0].header['BACKSUB_TYPE'] = 'Polynomial2D, degree=1'
                        for i in range(len(p.param_names)):
                            hdul[0].header[f'BACKSUB_PARAMS_{i}'] = p.param_names[i]
                            hdul[0].header[f'BACKSUB_PARAMVALS_{i}'] = p.parameters[i]
                        hdul.close()

                if back_model == 'median':
                    median = np.median(sky)
                    with fits.open(fn, mode='update') as hdul:
                        hdul[0].data -= median
                        hdul[0].header['BACKSUB_TYPE'] = 'median'
                        hdul[0].header['BACKSUB_VAL'] = median
                        hdul.close()

            logger.info('Solving field for ' + fn)
            solution = solve_field(fn, index_path=index_path,
                                   tmp_path=tmp_path,
                                   target_radec=center,
                                   search_radius=0.5,
                                   identifier='OBJECT')
            fn_base = fn.split('/')[-1].split('.proc.fits')[0]
            utils.mkdir_if_needed(os.path.join(frame_out,'catalogs'))
            cat_fn = os.path.join(frame_out,'catalogs',f'{fn_base}.cat')
            cat = extract_bright_stars(fn,catalog_path=cat_fn)
            check = check_astrometric_solution(ref_cat,
                                               header=solution.header,
                                               cat=cat, max_sep=1)
            tweak = TweakWCS(solution.header, check.cat_match, check.ref_match)
            tweak.optimize()
            tweak.update_header(solution.header)
            tweak.update_header(solution.fitsio_header)
            astrom.append(solution)
            if make_plots:
                fig_dir = os.path.join(frame_out, 'diagnostic-plots')
                utils.mkdir_if_needed(fig_dir)
                check = check_astrometric_solution(ref_cat,
                                                   header=solution.header,
                                                   cat=cat, max_sep=1,
                                                   make_plot=make_plots,
                                                   xlim=[-0.3, 0.3],
                                                   ylim=[-0.3, 0.3])
                fig_fn = os.path.basename(fn).replace('.fits',
                                                      '_check_wcs.png')
                check.fig.savefig(os.path.join(fig_dir, fig_fn), dpi=250)


        logger.end_tqdm()

        # resample images
        logger.info(f'Registering and writing {ftype} frames')
        reg_files = []
        for fn, sol in zip(files, astrom):
            resamp = improc.resample_image(
                fn, center[0], center[1], pixscale, args.output_dim[0],
                args.output_dim[1], sol.fitsio_header)
            header = fits.getheader(fn)
            for k, v in resamp.wcs.to_header().items():
                header[k] = v
            out_fn = os.path.basename(fn).replace('.fits', '_reg.fits')
            out_fn = os.path.join(frame_out, out_fn)
            fits.writeto(out_fn, resamp.pixels, header, overwrite=True)
            reg_files.append(out_fn)

        # stack images
        if stack:
            logger.info(f'Stacking and writing {ftype} frames')

            stack_data = np.ndarray((num,output_dim[1],output_dim[0]))
            exposure_map = np.zeros((output_dim[1],output_dim[0]))
            im_num = 0
            logger.start_tqdm()
            for fn in tqdm(reg_files):
                hdul = fits.open(fn)
                exptime = hdul[0].header['EXPTIME']
                stack_data[im_num] = hdul[0].data
                ind = np.nonzero(hdul[0].data)
                exposure_map[ind] += exptime
                hdul.close()
                im_num+=1
            logger.end_tqdm()
            # find median image
            if stack_type == 'median':
                stacked = np.median(stack_data,axis=0)

            hdu = fits.PrimaryHDU(stacked)
            hdu.writeto(os.path.join(frame_out,f'lbc{bandpass.lower()}_M81blob_{stack_type}stack.fits'),overwrite=True)
            hdu_exp = fits.PrimaryHDU(exposure_map)
            hdu_exp.writeto(os.path.join(frame_out,f'lbc{bandpass.lower()}_M81blob_{stack_type}stack_exposuremap.fits'),overwrite=True)



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument('-d', '--data-path', required=True,
                        help='path to data, which have had basic reductions')
    parser.add_argument('-o', '--out-path', required=True,
                        help='send all output here')
    parser.add_argument('-c', '--center', type=float, nargs=2, required=True,
                        help='central ra and dec')
    parser.add_argument('-t', '--tmp-path', default='/tmp',
                        help='temporary path for intermediate steps')
    parser.add_argument('-b', '--bandpass', default='R')
    parser.add_argument('--subtract-back', default=True)
    parser.add_argument('--back-model', default='plane',
                        help='background model type (median, plane)')
    parser.add_argument('-s','--stack-images', default=True)
    parser.add_argument('--stack-type', default='median',
                        help='stacking type (median, sigclipmean)')
    parser.add_argument('-r', '--ref-cat', help='reference catalog file name')
    parser.add_argument('--output-dim', default=[3000, 5500], nargs=2,
                        type=float, help='output image dimensions')
    parser.add_argument('--index-path', default=None,
                        help='path to astrometry.net index files')
    parser.add_argument('--pixscale', default=0.226, type=float)
    parser.add_argument('--make-plots', action='store_true')
    parser.add_argument('--log-level', default='info',
                        help='log level (debug, info, warn, error, critical)')
    parser.add_argument('--log-fn', default=None, help='log file name')
    args = parser.parse_args()
    utils.setup_logger(args.log_level, args.log_fn)

    register_images(
        data_path=args.data_path,
        out_path=args.out_path,
        center=args.center,
        tmp_path=args.tmp_path,
        bandpass=args.bandpass,
        subtract_back=args.subtract_back,
        back_model=args.back_model,
        stack=args.stack_images,
        stack_type=args.stack_type,
        ref_cat=args.ref_cat,
        index_path=args.index_path,
        output_dim=args.output_dim,
        pixscale=args.pixscale,
        make_plots=args.make_plots
    )
