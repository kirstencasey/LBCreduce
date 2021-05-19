import numpy as np
import os, yaml
from lbcred.log import logger
from lbcred import tools, sbf_functions
from lbcred.utils import misc, io
from lbcred.model import imfit, artpop_functions
from astropy.io import fits


def modeling_imfit(config_filename, options = {}, iter=None, fn_stub=None):

    logger.info('Getting ready...')

    # Open and read config file
    with open(config_filename, 'r') as filename:
        config = yaml.load(filename, Loader=yaml.FullLoader)

	# Replace options input via command line into config
    for key in options:
        if options[key] != None: config[key] = options[key]
    tools.setup_logger(config['logger_level'], log_fn=config['log_to_file'])

    # Make sure images have odd number of pixels on a side
    if config['check_image_size']:
        misc.save_image_odd_shape(os.path.join(config['image_dir'],config['color1']['image_fn']))
        misc.save_image_odd_shape(os.path.join(config['image_dir'],config['color2']['image_fn']))

    # Get Legacy Survey images if necessary
    if config['use_legacy_survey']:
        logger.info('Getting Legacy Survey Images...')
        # Get psf
        psf_fns = misc.fetch_psf(config['legacy_survey']['ra'], config['legacy_survey']['dec'], save_files=True, fn_root=config['legacy_survey']['fn_root'], out_dir=config['out_dir'])

        # Get images, inverse variance images
        image_fns, invvar_fns = misc.fetch_cutout(config['legacy_survey']['ra'], config['legacy_survey']['dec'], save_files=True, fn_root=config['legacy_survey']['fn_root'], out_dir=config['out_dir'])

    # Run ArtPop if necessary
    if config['run_artpop']:
        logger.info('Running ArtPop...')

        if config['use_legacy_survey']:
            model1, model2, src = artpop_functions.run_idealimager(config, image_fns, psf_fns) #################### DO THIS STILL ####################

        else:
            model1, model2, src = artpop_functions.run_artimager(config)

        if not config['run_imfit']:
            color1 = config['color1']['name']
            color2 = config['color2']['name']
            src_color1 = src.sp.total_mag('Bessell_R')
            src_color2 = src.sp.total_mag('Bessell_B')
            src_sbfmag1 = src.sp.sbf_mag('Bessell_R')
            logger.info(f'ArtPop model : \n{color1}-band magnitude: {src_color1}\n{color2}-band magnitude: {src_color2}\n{color2}-{color1} color: {src_color2-src_color1}\n{color1}-band SBF magnitude: {src_sbfmag1}\n')

    # Inject model into image(s) if necessary
    color1 = config['color1']['name']
    color2 = config['color2']['name']
    if config['inject_artpop_model']:
        logger.info('Injecting model galaxy into image...')
        if config['use_legacy_survey']:
            image1 = misc.inject_model(os.path.join(config['out_dir'],image_fns['r']), os.path.join(config['out_dir'],config['color1']['artpop_model_fn']), config['ypos_inject'],config['xpos_inject'])
            config['color1']['image_fn'] = image_fns['r'].replace('_r.fits', '_mock_injected_r.fits')
            io.write_pixels(os.path.join(config['out_dir'],config['color1']['image_fn']), image1)

            image2 = misc.inject_model(os.path.join(config['out_dir'],image_fns['g']), os.path.join(config['out_dir'],config['color2']['artpop_model_fn']), config['ypos_inject'],config['xpos_inject'])
            config['color2']['image_fn'] = image_fns['g'].replace('_g.fits', '_mock_injected_g.fits')
            io.write_pixels(os.path.join(config['out_dir'],config['color2']['image_fn']), image2)

        else:
            image1 = misc.inject_model(os.path.join(config['out_dir'],config['color1']['image_fn']), os.path.join(config['out_dir'],config['color1']['artpop_model_fn']), config['ypos_inject'],config['xpos_inject'])
            config['color1']['original_fn'] = config['color1']['image_fn']
            config['color1']['image_fn'] = config['color1']['image_fn'].replace(f'_{color1}.fits', f'_mock_injected_{color1}.fits')
            io.write_pixels(os.path.join(config['out_dir'],config['color1']['image_fn']), image1)


            image2 = misc.inject_model(os.path.join(config['out_dir'],config['color2']['image_fn']), os.path.join(config['out_dir'],config['color2']['artpop_model_fn']), config['ypos_inject'],config['xpos_inject'])
            config['color2']['original_fn'] = config['color2']['image_fn']
            config['color2']['image_fn'] = config['color2']['image_fn'].replace(f'_{color2}.fits', f'_mock_injected_{color2}.fits')
            io.write_pixels(os.path.join(config['out_dir'],config['color2']['image_fn']), image2)

        config['image_dir'] = config['out_dir']

    # Run imfit if necessary
    if config['run_imfit']:

        logger.info('Creating Imfit Masks...')
        # Get imfit mask
        mask1, mask1_fn = imfit.create_imfit_mask(config, config['color1'])
        mask2, mask2_fn = imfit.create_imfit_mask(config, config['color2'])

        if config['combine_imfit_masks']:
            final_mask = mask1.astype(bool) | mask2.astype(bool)
            final_mask = final_mask.astype(float)
            mask_hdu = fits.PrimaryHDU(final_mask)
            mask_hdu.writeto(mask1_fn,overwrite=True)
            mask_hdu.writeto(mask2_fn,overwrite=True)
        '''
        # Smooth images if necessary
        if config['smooth_images']:
            im1_fn = imfit.smooth_image(config['color1']['image_fn'], config)
            im2_fn = imfit.smooth_image(config['color2']['image_fn'], config)


        else:
            im1_fn = config['color1']['image_fn']
            im2_fn = config['color2']['image_fn']
        '''
        rdnoise = config['readnoise']
        gain = config['gain']
        ncomb = config['ncombined']
        '''
        sky1 = config['color1']['sky']
        sky2 = config['color2']['sky']
        options = f'--readnoise {rdnoise} --gain {gain} --ncombined {ncomb}'
        options1 = options + f'--sky {sky1} '
        options2 = options + f'--sky {sky2} '

        # run imfit, b-band sersic only
        logger.info('Running Imfit : Sersic-only')
        #results2, model2_fn, resid2_fn = imfit.run_imfit(im2_fn, mask2_fn, config['color2'], config, options2, sersic=True, tiltedplane=False, fixedsersic=None, viz=True, iter=iter, fn_stub=fn_stub)

        # run imfit, r-band fixed sersic + tilted plane
        logger.info('Running Imfit : Sersic + Tilted Plane')
        #results1, model1_fn, resid1_fn = imfit.run_imfit(im1_fn, mask1_fn, config['color1'], config, options1, sersic=True, tiltedplane=True, fixedsersic=results2, viz=True, iter=iter, fn_stub=fn_stub)
        results1, model1_fn, resid1_fn = imfit.run_imfit(im1_fn, mask1_fn, config['color1'], config, options1, sersic=True, tiltedplane=False, fixedsersic=None, viz=True, iter=iter, fn_stub=fn_stub)

        # run imfit, b-band fixed sersic + tilted plane
        #results2, model2_fn, resid2_fn = imfit.run_imfit(im2_fn, mask2_fn, config['color2'], config, options2, sersic=True, tiltedplane=True, fixedsersic=results2, viz=True, iter=iter, fn_stub=fn_stub)
        results2, model2_fn, resid2_fn = imfit.run_imfit(im2_fn, mask2_fn, config['color2'], config, options2, sersic=True, tiltedplane=False, fixedsersic=results1, viz=True, iter=iter, fn_stub=fn_stub)

        # Summarize major findings (mags, color, etc.)
        mag1, mag2, color = imfit.summarize_results(config, results1.results['comp_1'], results2.results['comp_1'])


        '''
        #### New run_imfit steps

        imfit_results = {}
        step_num = 1
        for key,step in config['imfit_steps'].items():

            # Initialize results dict, get correct color info, filenames, etc.
            imfit_results[key] = {}
            options = f'--readnoise {rdnoise} --gain {gain} --ncombined {ncomb}'
            if step['color'] == config['color1']['name']:
                color_options = config['color1']
                mask_fn = mask1_fn
                final_step1 = key

            else:
                color_options = config['color2']
                mask_fn = mask2_fn
                final_step2 = key

            im_fn = color_options['image_fn']

            if config['smooth_images']:
                im_fn = imfit.smooth_image(im_fn, config)

            sky = color_options['sky']
            options += f'--sky {sky} '
            imfit_results[key]['functions'] = misc.list_of_strings(step['functions'])

            if step_num == 1 or step['fix_struc_step'] is None:
                fixed_struc = None
            else:
                fixed_struc = imfit_results[step['fix_struc_step']]['results']

            # Run imfit
            model_sum = ''
            for func in imfit_results[key]['functions']: model_sum += f'{func} '
            logger.info(f'Running Imfit : {model_sum}')
            results, model_fn, resid_fn = imfit.run_imfit_updated(im_fn, mask_fn, color_options, config, model_funcs=imfit_results[key]['functions'], options=options, fixedsersic=fixed_struc, viz=True, iter=iter, fn_stub=fn_stub)

            # Update results dict
            imfit_results[key]['results'] = results
            imfit_results[key]['model_fn'] = model_fn
            imfit_results[key]['resid_fn'] = resid_fn
            step_num += 1

        # Summarize major findings
        comp1 = imfit_results[final_step1]['functions'].index('Sersic')+1
        comp2 = imfit_results[final_step2]['functions'].index('Sersic')+1
        results1 = imfit_results[final_step1]['results'].results[f'comp_{comp1}']
        results2 = imfit_results[final_step2]['results'].results[f'comp_{comp2}']
        mag1, mag2, color = imfit.summarize_results(config, results1, results2)

        ######### END

        logger.info(f'Imfit result : \n{color1}-band magnitude: {mag1}\n{color2}-band magnitude: {mag2}\n{color2}-{color1} color: {color}\n')

        if config['run_artpop']:
            src_color1 = src.sp.total_mag('Bessell_R')
            src_color2 = src.sp.total_mag('Bessell_B')
            src_sbfmag1 = src.sp.sbf_mag('Bessell_R')
            logger.info(f'ArtPop model : \n{color1}-band magnitude: {src_color1}\n{color2}-band magnitude: {src_color2}\n{color2}-{color1} color: {src_color2-src_color1}\n{color1}-band SBF magnitude: {src_sbfmag1}\n')
            logger.info(f'Errors : \n{color1}-band magnitude: {mag1-src_color1}\n{color2}-band magnitude: {mag2-src_color2}\n{color2}-{color1} color: {color - (src_color2-src_color1)}\n')
        elif config['inject_artpop_model'] and config['color1']['artpop_mag'] is not None and config['color2']['artpop_mag'] is not None:
            src_color1 = config['color1']['artpop_mag']
            src_color2 = config['color2']['artpop_mag']
            src_sbfmag1 = config['color1']['artpop_sbfmag']
            logger.info(f'ArtPop model : \n{color1}-band magnitude: {src_color1}\n{color2}-band magnitude: {src_color2}\n{color2}-{color1} color: {src_color2-src_color1}\n{color1}-band SBF magnitude: {src_sbfmag1}\n')
            logger.info(f'Errors : \n{color1}-band magnitude: {mag1-src_color1}\n{color2}-band magnitude: {mag2-src_color2}\n{color2}-{color1} color: {color - (src_color2-src_color1)}\n')
        else:
            src_sbfmag1 = None
            src_color1 = None
            src_color2 = None

        ''' OLD IMFIT
        if iter is None:
            bf1_fn = os.path.join(config['out_dir'],model1_fn.split('_model')[0]+f'_bestfit-params_{color1}.txt')
            bf2_fn = os.path.join(config['out_dir'],model2_fn.split('_model')[0]+f'_bestfit-params_{color2}.txt')
        else:
            bf1_fn = os.path.join(config['out_dir'],model1_fn.split('_model')[0]+f'_bestfit-params_{fn_stub}_iter{iter}_{color1}.txt')
            bf2_fn = os.path.join(config['out_dir'],model2_fn.split('_model')[0]+f'_bestfit-params_{fn_stub}_iter{iter}_{color2}.txt')
        #bestfit1 = io.read_results(bf1_fn, ['TiltedSkyPlane', 'Sersic'])
        bestfit1 = io.read_results(bf1_fn, ['Sersic'])
        #bestfit2 = io.read_results(bf2_fn, ['TiltedSkyPlane', 'Sersic'])
        bestfit2 = io.read_results(bf2_fn, ['Sersic'])
        #bestfit1 = io.read_results(bf1_fn, ['FlatSky', 'Sersic'])
        #bestfit2 = io.read_results(bf2_fn, ['FlatSky', 'Sersic'])

        return bestfit1, bestfit2, bf1_fn.split('/')[-1], bf2_fn.split('/')[-1], mag1, mag2, color, model1_fn, resid1_fn, model2_fn, resid2_fn, src_sbfmag1, src_color1, src_color2

        '''
        ######### NEW
        if iter is None:
            bf1_fn = os.path.join(config['out_dir'],imfit_results[final_step1]['model_fn'].split('_model')[0]+f'_bestfit-params_{color1}.txt')
            bf2_fn = os.path.join(config['out_dir'],imfit_results[final_step2]['model_fn'].split('_model')[0]+f'_bestfit-params_{color2}.txt')
        else:
            bf1_fn = os.path.join(config['out_dir'],imfit_results[final_step1]['model_fn'].split('_model')[0]+f'_bestfit-params_{fn_stub}_iter{iter}_{color1}.txt')
            bf2_fn = os.path.join(config['out_dir'],imfit_results[final_step2]['model_fn'].split('_model')[0]+f'_bestfit-params_{fn_stub}_iter{iter}_{color2}.txt')

        bestfit1 = io.read_results(bf1_fn, imfit_results[final_step1]['functions'])
        bestfit2 = io.read_results(bf2_fn, imfit_results[final_step2]['functions'])

        return bestfit1, bestfit2, bf1_fn.split('/')[-1], bf2_fn.split('/')[-1], mag1, mag2, color, imfit_results[final_step1]['model_fn'], imfit_results[final_step1]['resid_fn'], imfit_results[final_step1]['functions'], imfit_results[final_step2]['model_fn'], imfit_results[final_step2]['resid_fn'], src_sbfmag1, src_color1, src_color2
        ######### END New

    return
