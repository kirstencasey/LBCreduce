import numpy as np
import os, yaml, pandas
from shutil import copyfile
from lbcred.log import logger
from lbcred import tools, sbf_functions
from lbcred.utils import misc, io
from lbcred.model import imfit
from astropy.io import fits
from astropy.table import Table


def modeling_sbf(config_filename, options = {}, run_iter=None, imfit_functions=None, run_id=None, plot_blank_fields=False):

    logger.info('Getting ready...')
    if type(config_filename) is not dict:
        # Open and read config file
        with open(config_filename, 'r') as filename:
            config = yaml.load(filename, Loader=yaml.FullLoader)
        # Save copy of config to output directory
        copyfile(config_filename,os.path.join(config['out_dir'],config_filename.split('/')[-1]))
    else:
        config = config_filename
        # Save copy of config to output directory
        with open(os.path.join(config['out_dir'],'modeling-config_sbf.yml'), 'w') as outfile:
            yaml.dump(config, outfile, default_flow_style=False)

	# Replace options input via command line into config
    for key in options:
        if options[key] != None: config[key] = options[key]
    tools.setup_logger(config['logger_level'], log_fn=config['log_to_file'])

    # Create residual if necessary
    if config['create_resid_from_model']:
        model = fits.getdata(os.path.join(config['image_dir'],config['model_fn']), ext=config['ext'])
        orig_im = fits.open(os.path.join(config['image_dir'],config['orig_image_fn']), ext=config['ext'])
        resid = orig_im.copy()
        resid[0].data = orig_im[0].data - model
        resid.writeto(os.path.join(config['image_dir'],config['resid_fn']),overwrite=config['overwrite'])

    # Make sure images have odd number of pixels on a side
    shape = misc.save_image_odd_shape(os.path.join(config['image_dir'],config['model_fn']), ext=config['ext'])
    misc.save_image_odd_shape(os.path.join(config['image_dir'],config['resid_fn']),ext=config['ext'])

    psf = fits.open(os.path.join(config['image_dir'],config['psf']))[config['ext']].data
    psf = psf/np.sum(psf)

    if config['model_summary_fn'] != None:
        bestfit_params = io.read_results(os.path.join(config['image_dir'], config['model_summary_fn']), imfit_functions)
        sersic_idxs = np.where(np.asarray(imfit_functions)=='Sersic')[0]
        
        if len(sersic_idxs) > 1: # Choose Sersic component centered closest to the center of the image 
            pix_dist = 10000
            for sersic_idx in sersic_idxs:
                sersic = bestfit_params[f'comp_{sersic_idx+1}']
                sersic_dist = np.sqrt((sersic['X0']-shape[0]/2)**2 + (sersic['Y0']-shape[1]/2)**2)
                if sersic_dist < pix_dist:
                    pix_dist = sersic_dist
                    comp = sersic_idx+1
        
        else: comp = imfit_functions.index('Sersic')+1
        sersic_results = bestfit_params[f'comp_{comp}']
    
    else:
        sersic_results = {'function' : 'Sersic', 'X0' : config['xpos'], 'Y0' : config['ypos'],
                           'PA' : config['pa'], 'ell' : config['ellip'], 'n' : config['n'],
                           'I_e' : config['I_e'], 'r_e' : config['radius']}

    model_cutout_fn = config['model_fn'].replace('.fits','_cutout.fits')
    resid_cutout_fn = config['resid_fn'].replace('.fits','_cutout.fits')

    if sersic_results['ell'] < 0: smajor_axis = sersic_results['r_e']*(1-sersic_results['ell'])
    else: smajor_axis = sersic_results['r_e']

    if config['masking_sbf']['randomly_vary_mask_radius']:
        cutout_size = 2*(smajor_axis*config['masking_sbf']['max_random_frac_of_radius'])//1
        if cutout_size < psf.shape[0]:
            cutout_size = psf.shape[0]
        if cutout_size%2==0: cutout_size+=1
    else:
        cutout_size = 2*(smajor_axis)//1
        if cutout_size < psf.shape[0]:
            cutout_size = psf.shape[0]
        if cutout_size%2==0: cutout_size+=1
    cutout = misc.make_cutout(os.path.join(config['image_dir'],config['model_fn']), (sersic_results['X0'],sersic_results['Y0']), (cutout_size,cutout_size), ext=config['ext'], cutout_fn=os.path.join(config['image_dir'],model_cutout_fn))
    misc.make_cutout(os.path.join(config['image_dir'],config['resid_fn']), (sersic_results['X0'],sersic_results['Y0']), (cutout_size,cutout_size), ext=config['ext'], cutout_fn=os.path.join(config['image_dir'],resid_cutout_fn))
    config['model_fn'] = model_cutout_fn
    config['resid_fn'] = resid_cutout_fn
    sersic_results['X0'] = cutout.center_cutout[0]
    sersic_results['Y0'] = cutout.center_cutout[1]
    
    # Get blank field info if necessary
    if config['include_blank_field_correction']:
         blank_field_positions = Table.from_pandas(pandas.read_csv(config['pre_selected_fields']))

    # run sbf measurement
    color = config['color_name'].upper()
    chi_squares = []
    sbf_mags = []
    dists_a = []
    dists_b = []
    k_mins = []
    k_maxs = []
    all_results = []
    blank_results = None
    used_blank_nums = []
    blank_inputs = {}
    signals = []
    bg_signals = []
    
    if config['masking_sbf']['randomly_vary_grow_obj']: grow_objs =[]
    if config['masking_sbf']['randomly_vary_mask_radius']: scales = []
    iter=0
    logger.info('Iterating through SBF measurement parameters...')
    while iter < config['num_iters']:
        if config['masking_sbf']['randomly_vary_mask_radius']:
            scale = np.random.uniform(config['masking_sbf']['min_random_frac_of_radius'],config['masking_sbf']['max_random_frac_of_radius'])
            scales.append(scale)
        else: scale=config['masking_sbf']['fixed_frac_of_radius']

        grow_obj = config['masking_sbf']['grow_obj']
        if config['masking_sbf']['randomly_vary_grow_obj']:
            grow_obj = np.random.uniform(config['masking_sbf']['grow_obj']-config['masking_sbf']['random_growth_max_deviation'],config['masking_sbf']['grow_obj']+config['masking_sbf']['random_growth_max_deviation'])
            grow_objs.append(grow_obj)
        else: grow_obj=config['masking_sbf']['grow_obj']
        
        # Get normalized residual and mask for sbf measurement
        sbf_resid, sbf_mask = sbf_functions.get_sbf_mask_resid(os.path.join(config['image_dir'],config['model_fn']), os.path.join(config['image_dir'],config['resid_fn']), sersic_results, grow_obj, scale, config)

        # Get normalized residual and mask for sbf measurement on blank field
        if config['include_blank_field_correction']:
            # Select random blank field from csv, prep cutouts
            blank_num = np.random.randint(0,len(blank_field_positions))
            pos = blank_field_positions[blank_num]
            misc.make_cutout(config['stack_fn'], (pos['xpos'],pos['ypos']), (cutout_size,cutout_size), ext=config['ext'], cutout_fn=os.path.join(config['image_dir'],f'blank_field_{blank_num}_r.fits'))
            
            # Get resid, mask
            blank_resid, blank_mask = sbf_functions.get_sbf_mask_resid(os.path.join(config['image_dir'],config['model_fn']), os.path.join(config['image_dir'],f'blank_field_{blank_num}_r.fits'), sersic_results, grow_obj, scale, config, blank_field=True)
            if blank_num not in used_blank_nums: 
                blank_inputs[blank_num] = {'resid' : blank_resid,
                                           'mask' : blank_mask}

        if config['randomly_vary_k_limits']:
            kmins = misc.list_of_floats(config['k_min'])
            kmaxs = misc.list_of_floats(config['k_max'])
            k_min = np.random.uniform(kmins[0],kmins[1])
            k_max = np.random.uniform(kmaxs[0],kmaxs[1])
        else: 
            k_min = float(config['k_min'])
            k_max = float(config['k_max'])
        k_mins.append(k_min)
        k_maxs.append(k_max)

        results = sbf_functions.measure_sbf(sbf_resid[config['ext']].data, psf, mask=sbf_mask, k_range=[k_min, k_max],
                        fit_param_guess=[100, 50], num_radial_bins=config['num_radial_bins'],
                        use_sigma=config['use_sigma'])
        signal = (results.p[0])*config['gain']/config['exposure_time']/results.npix
        
        bg_signal = None
        if config['include_blank_field_correction']:
            blank_results = sbf_functions.measure_sbf(blank_resid[config['ext']].data, psf, mask=blank_mask, k_range=[k_min, k_max],
                        fit_param_guess=[100, 50], num_radial_bins=config['num_radial_bins'],
                        use_sigma=config['use_sigma'])
            signal = (results.p[0]-blank_results.p[0])*config['gain']/config['exposure_time']/results.npix
            bg_signal = (blank_results.p[0])*config['gain']/config['exposure_time']/results.npix
            used_blank_nums.append(blank_num)

        # Should probably modify this in the case where blank_results!=None
        chisq = misc.get_chisquare(f_obs=(results.ps_image / results.npix), f_exp= (results.fit_func(results.k, *results.p) / results.npix))
        chi_squares.append(chisq)

        sbf_mag, d_a, d_b = sbf_functions.get_sbf_distance(results, config['zpt'], config['color'], config['gain'], config['exposure_time'], colorterm = config['color_term'], extinction_correction=config['extinction'], blank_field_results=blank_results)

        sbf_mags.append(sbf_mag)
        dists_a.append(d_a)
        dists_b.append(d_b)
        all_results.append(results)
        signals.append(signal)
        bg_signals.append(bg_signal)


        iter+=1


    chi_squares=np.asarray(chi_squares)
    if config['include_blank_field_correction']:
        wavg_mag = np.median(sbf_mags)
        wavg_dist_a = np.median(dists_a)
        wavg_dist_b = np.median(dists_b)
        uncert_mag = np.std(sbf_mags)
        print('\nMedian signal: ',np.median(np.asarray(signals)))
        print('Std of background signal: ',np.std(np.asarray(bg_signals)))
        print('S/N: ',np.median(np.asarray(signals))/np.std(np.asarray(bg_signals)),'\n')
        s_n = np.median(np.asarray(signals))/np.std(np.asarray(bg_signals))
    else:
        wavg_mag = np.average(sbf_mags,weights=1/chi_squares)
        wavg_dist_a = np.average(dists_a,weights=1/chi_squares)
        wavg_dist_b = np.average(dists_b,weights=1/chi_squares)
        uncert_mag = np.std(sbf_mags)
        s_n = 0.
    logger.info(f'SBF result:\nSBF magnitude: {round(wavg_mag,3)}\nSBF distance (using Eqn. 2 from Jerjen 2000): {round(wavg_dist_a/1e6,3)} Mpc\nSBF distance (using Eqn. 3 from Jerjen 2000): {round(wavg_dist_b/1e6,3)} Mpc')

    k_mins=np.asarray(k_mins)
    k_maxs=np.asarray(k_maxs)

    # Save results
    for name, arr in zip(['sbf_mags','chi_squares','k_min','k_max','dist_a','dist_b','s_n'],[sbf_mags,chi_squares,k_mins,k_maxs,dists_a,dists_b,s_n]):
        out_name = f'sbf_calculation_{name}'
        if run_id is not None: out_name += f'_{run_id}'
        if run_iter is not None: out_name += f'_iter{run_iter}'

        file = open(os.path.join(config['out_dir'],out_name), "wb")
        np.save(file, arr)
        file.close

    idx = np.where(chi_squares==chi_squares.min())[0][0]
    logger.info('Calculating best-fit SBF...')
    k_min = k_mins[idx]
    k_max = k_maxs[idx]
    if config['masking_sbf']['randomly_vary_grow_obj']:
        grow_objs = np.asarray(grow_objs)
        out_name = 'sbf_calculation_grow_obj'
        if run_id is not None: out_name += f'_{run_id}'
        if run_iter is not None: out_name += f'_iter{run_iter}'
        file = open(os.path.join(config['out_dir'],out_name), "wb")
        np.save(file, grow_objs)
        file.close
        grow_obj = grow_objs[idx]
        print('Grow obj: ', grow_obj)
    if config['masking_sbf']['randomly_vary_mask_radius']:
        scales = np.asarray(scales)
        out_name = 'sbf_calculation_mask_radius_scale'
        if run_id is not None: out_name += f'_{run_id}'
        if run_iter is not None: out_name += f'_iter{run_iter}'
        file = open(os.path.join(config['out_dir'],out_name), "wb")
        np.save(file, scales)
        file.close
        scale = scales[idx]
        print('Scale : ', scale)
    #print('kmin, kmax: ', k_min, k_max)

    model_fn = os.path.join(config['image_dir'],config['model_fn'])
    resid_fn = os.path.join(config['image_dir'],config['resid_fn'])
    sbf_resid, sbf_mask = sbf_functions.get_sbf_mask_resid(model_fn, resid_fn, sersic_results, grow_obj, scale, config)
    results = sbf_functions.measure_sbf(sbf_resid[config['ext']].data, psf, mask=sbf_mask, k_range=[k_min, k_max],
                    fit_param_guess=[100, 50], num_radial_bins=config['num_radial_bins'],
                    use_sigma=config['use_sigma'])

    sbf_resid[config['ext']].data[sbf_mask.astype(bool)] = 0.0
    if run_iter == None:
        save_fn = os.path.join(config['out_dir'],config['model_fn'].replace('.fits','_sbf-results.png'))
    else:
        save_fn = os.path.join(config['out_dir'],config['model_fn'].replace('.fits',f'_sbf-results_iter{run_iter}.png'))

    sbf_functions.sbf_results(results, sbf_resid[config['ext']].data, subplots=None, xlabel=r'Spacial Frequency (pixel$^{-1}$)',
                    ylabel=f'Power ({color}-band)', xscale='linear', percentiles=[config['plot_percentiles_min'], config['plot_percentiles_max']],
                    yscale='log', plot_errors=False, ylim_factors=[0.5, 1.1],
                    cmap='gray_r',save_fn=save_fn)

    # Plot SBF mag
    if plot_blank_fields:
        for key in blank_inputs:
            
            blank_result = sbf_functions.measure_sbf(blank_inputs[key]['resid'][config['ext']].data, psf, mask=blank_inputs[key]['mask'], k_range=[k_min, k_max],
                    fit_param_guess=[100, 50], num_radial_bins=config['num_radial_bins'],
                    use_sigma=config['use_sigma'])
            blank_inputs[key]['resid'][config['ext']].data[blank_inputs[key]['mask'].astype(bool)] = 0.0
            blank_inputs[key]['results'] = blank_result
            
        # Make it so that sbf_results function plots all the blank fields.
        sbf_functions.sbf_results(results, sbf_resid[config['ext']].data, subplots=None, xlabel=r'Spacial Frequency (pixel$^{-1}$)',
                ylabel=f'Power (${color}$-band)', xscale='linear', percentiles=[config['plot_percentiles_min'], config['plot_percentiles_max']],
                yscale='log', plot_errors=False, ylim_factors=[0.5, 1.1],
                cmap='gray_r',save_fn=save_fn.replace('.png','_withblankfields.png'), normalize_ps = True, plot_blank_fields = True, blank_results = blank_inputs)
        
    return wavg_mag, wavg_dist_a, wavg_dist_b, uncert_mag
