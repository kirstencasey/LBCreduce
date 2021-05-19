import numpy as np
import os, yaml
from lbcred.log import logger
from lbcred import tools, sbf_functions
from lbcred.utils import misc
from lbcred.model import imfit
from astropy.io import fits


def modeling(config_filename, options = {}):

    logger.info('Getting ready...')
    # Open and read config file
    with open(config_filename, 'r') as filename:
        config = yaml.load(filename, Loader=yaml.FullLoader)

	# Replace options input via command line into config
    for key in options:
        if options[key] != None: config[key] = options[key]
    tools.setup_logger(config['logger_level'], log_fn=config['log_to_file'])

    # inject model into image if necessary
    color1 = config['color1']
    color2 = config['color2']
    if config['inject_model1'] is not None:
        image1 = misc.inject_model(os.path.join(config['image_dir'],config['image1']), config['inject_model1'], config['xpos'], config['ypos'])
        config['image1'] = config['image1'].replace(f'_{color1}.fits', f'_mock_injected_{color1}.fits')
        io.write_pixels(os.path.join(config['image_dir'],config['image1']), image1)

    if config['inject_model2'] is not None:
        image2 = misc.inject_model(os.path.join(config['image_dir'],config['image2']), config['inject_model2'], config['xpos'], config['ypos'])
        config['image2'] = config['image2'].replace(f'_{color2}.fits', f'_mock_injected_{color2}.fits')
        io.write_pixels(os.path.join(config['image_dir'],config['image2']), image2)

    # Get imfit mask
    mask1, mask1_fn = imfit.create_imfit_mask(config, os.path.join(config['image_dir'],config['image1']), color1, config['color1_masking_imfit'])
    mask2, mask2_fn = imfit.create_imfit_mask(config, os.path.join(config['image_dir'],config['image2']), color2, config['color2_masking_imfit'])

    # run imfit, b-band sersic only
    rdnoise = config['readnoise']
    gain = config['gain']
    ncomb = config['ncombined']
    sky1 = config['sky1']
    sky2 = config['sky2']
    options = f'--readnoise {rdnoise} --gain {gain } --ncombined {ncomb}'
    options1 = options + f'--sky {sky1} '
    options2 = options + f'--sky {sky2} '
    sersic_results2, _, _ = imfit.run_imfit(os.path.join(config['image_dir'],config['image2']), mask2_fn, os.path.join(config['image_dir'],config['psf2']), config['color2'], config, options2, sersic=True, tiltedplane=False, fixedsersic=None, viz=True)

    # run imfit, r-band fixed sersic + tilted plane
    sersic_tp_results1, model1_fn, resid1_fn = imfit.run_imfit(os.path.join(config['image_dir'],config['image1']), mask1_fn, os.path.join(config['image_dir'],config['psf1']), config['color1'], config, options1, sersic=True, tiltedplane=True, fixedsersic=sersic_results2, viz=True)

    # run imfit, b-band fixed sersic + tilted plane
    sersic_tp_results2, model2_fn, resid2_fn = imfit.run_imfit(os.path.join(config['image_dir'],config['image2']), mask2_fn, os.path.join(config['image_dir'],config['psf2']), config['color2'], config, options2, sersic=True, tiltedplane=True, fixedsersic=sersic_results2, viz=True)

    # Summarize major findings (mags, color, etc.)
    mag1, mag2, color = imfit.summarize_results(model1_fn, config['zpt1'], model2_fn, config['zpt2'], config, sersic_tp_results1, sersic_tp_results2)

    logger.info(f'Imfit result : \n{color1}-band magnitude: {mag1}\n{color2}-band magnitude: {mag2}\n{color2}-{color1} color: {color}\n')

    if config['sbf']:

        sbf_resid1, sbf_mask1 = sbf_functions.get_sbf_mask_resid(model1_fn, resid1_fn, sersic_results2.results['comp_1'], color1, config['color1_masking_sbf'], config)

        # run sbf measurement
        psf1 = fits.open(os.path.join(config['image_dir'],config['psf1']))[0].data
        psf1 = psf1/np.sum(psf1)
        psf2 = fits.open(os.path.join(config['image_dir'],config['psf2']))[0].data
        psf2 = psf2/np.sum(psf2)

        results1 = sbf_functions.measure_sbf(sbf_resid1[0].data, psf1, mask=sbf_mask1, k_range=[config['k_min'], config['k_max']],
                        fit_param_guess=[100, 50], num_radial_bins=config['num_radial_bins'],
                        use_sigma=config['use_sigma'])
        sbf_resid1[0].data[sbf_mask1.astype(bool)] = 0.0
        sbf_functions.sbf_results(results1, sbf_resid1[0].data, subplots=None, xlabel=r'Spacial Frequency (pixel$^{-1}$)',
                        ylabel=f'Power ({color1}-band)', xscale='linear', percentiles=[config['plot_percentiles_min'], config['plot_percentiles_max']],
                        yscale='log', plot_errors=True, ylim_factors=[0.5, 1.1],
                        cmap='gray_r', save_fn=model1_fn.replace(f'.fits',f'_sbf_results.png'))

        sbf_mag, d_a, d_b = sbf_functions.get_sbf_distance(results1, config['zpt1'], color, config['gain'], config['exposure_time'], colorterm = config['color_term1'])

        logger.info(f'SBF result:\nSBF magnitude: {round(sbf_mag,3)}\nSBF distance (using Eqn. 2 from Jerjen 2000): {round(d_a/1e6,3)} Mpc\nSBF distance (using Eqn. 3 from Jerjen 2000): {round(d_b/1e6,3)} Mpc')

    return
