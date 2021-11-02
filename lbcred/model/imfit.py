import pymfit, subprocess, os, yaml, pandas
import numpy as np
from lbcred.utils import io, misc, mkdir_if_needed
from astropy.io import fits
from astropy.table import Table
from collections import OrderedDict



def create_imfit_mask(config, color_specific_info):

    '''
    # Taken from sbf_functions.get_sbf_mask_resid
    mask_fn = model_fn.replace(f'.fits','_sbf_mask.fits')
    dist_mod = 5*np.log10(config['assumed_distance']*10**6)-5
    app_mag = dist_mod+config['given_abs_mag_to_mask']
    resid_no_norm = np.ascontiguousarray(resid_no_norm)
    obj, seg, bkg, img = pymfit.masking.detect_sources(resid_no_norm, config['masking_sbf']['thresh'], config['masking_sbf']['backsize'], mask=aper_mask, return_all=True, kern_sig=config['masking_sbf']['kern_sig'])
    flux=obj['flux']

    obj = obj[(-2.5*np.log10(flux*config['gain']/config['exposure_time'])+config['zpt']) < app_mag]
    seg_mask = pymfit.masking.make_seg_mask(seg, grow_obj, config['masking_sbf']['thresh'])
    obj_mask = pymfit.masking.make_obj_mask(obj, resid[0].data.shape, grow_obj)
    mask = (seg_mask | obj_mask).astype(int)
    fits.writeto(mask_fn, mask, overwrite=config['overwrite'])
    '''

    color = color_specific_info['name']
    if color_specific_info['mask_fn'] is not None :
        mask_fn = os.path.join(config['image_dir'],color_specific_info['mask_fn'])
        mask = fits.getdata(mask_fn)
        return mask, mask_fn

    else :
        mask_fn = os.path.join(config['out_dir'],color_specific_info['image_fn'].replace(f'_{color}.fits',f'_mask_{color}.fits'))

    if config['inject_artpop_model'] :
        x = config['xpos_inject']
        y = config['ypos_inject']
    else:
        x = config['sersic_params']['xpos_guess']
        y = config['sersic_params']['ypos_guess']
    mask_kws = dict(out_fn=mask_fn, thresh=color_specific_info['masking_imfit']['thresh'], kern_sig=color_specific_info['masking_imfit']['kern_sig'], backsize=color_specific_info['masking_imfit']['backsize'],
                        obj_rmin=color_specific_info['masking_imfit']['obj_rmin'], grow_obj=color_specific_info['masking_imfit']['grow_obj'], use_hsc_mask=False, gal_pos=(x,y), seg_rmin=color_specific_info['masking_imfit']['seg_rmin']) #, sep_extract_kws=sep_extract_kws)
    mask = pymfit.make_mask(os.path.join(config['image_dir'],color_specific_info['image_fn']), **mask_kws)

    if config['mask_bright_star']:

        r_circ = config['circle_mask']['radius']
        i_c, j_c = (config['circle_mask']['ypos'], config['circle_mask']['xpos'])
        ii, jj = np.mgrid[:mask.shape[0], :mask.shape[1]]
        circ_mask = ((ii - i_c)**2 + (jj - j_c)**2) - r_circ**2 < 0

        mask = mask.astype(bool) | circ_mask
        mask = mask.astype(float)
        mask_hdu = fits.PrimaryHDU(mask)
        mask_hdu.writeto(mask_fn,overwrite=True)

    if config['random_inject_position']:
        im = fits.open(os.path.join(config['image_dir'],color_specific_info['original_fn']))
        star_mask = fits.open('/Users/kirstencasey/m81blob_out/sci/lbc-reduce_testing/star_mask.fits')
        edge_mask = im[0].data == 0
        mask = mask.astype(bool) | edge_mask.astype(bool)
        mask = mask.astype(bool) | star_mask[0].data.astype(bool)
        mask = mask.astype(float)
        mask_hdu = fits.PrimaryHDU(mask)
        mask_hdu.writeto(mask_fn,overwrite=True)

    return mask, mask_fn

def run_makeimage(bestfit_fn, psf_fn=None, ref_fn=None, output_root=None, out_fn=None, del_temp=False, options=''):

    cmd = f'makeimage {bestfit_fn} '

    if psf_fn is not None:
        cmd += f'--psf {psf_fn} '

    if ref_fn is not None:
        cmd += f'--refimage {ref_fn} '

    if output_root is not None:
        cmd += f'--output-functions {output_root}'

    elif out_fn is not None:
        cmd += f'--output {out_fn}'

    cmd += options
    subprocess.call(cmd, shell=True)

    if del_temp: os.remove(temp_fn)

    return

def organize_initial_params(config, model, fixedsersic=None, color=None):


    if model == 'Sersic' and fixedsersic == None:
        if config['sersic_params']['fix_PA']: PA = [config['sersic_params']['PA_guess'], 'fixed']
        else: PA = [config['sersic_params']['PA_guess'], config['sersic_params']['PA_min'], config['sersic_params']['PA_max']]
        if config['sersic_params']['fix_n']: n = [config['sersic_params']['n_guess'], 'fixed']
        else: n = [config['sersic_params']['n_guess'], config['sersic_params']['n_min'], config['sersic_params']['n_max']]
        if config['sersic_params']['fix_ell']: ell = [config['sersic_params']['ell_guess'], 'fixed']
        else: ell = [config['sersic_params']['ell_guess'], config['sersic_params']['ell_min'], config['sersic_params']['ell_max']]
        if config['sersic_params']['fix_r_e']: r_e = [config['sersic_params']['r_e_guess'], 'fixed']
        else: r_e = [config['sersic_params']['r_e_guess'], config['sersic_params']['r_e_min'], config['sersic_params']['r_e_max']]
        if config['sersic_params']['fix_I_e']:
            if color==config['color1']['name']:
                I_e = [45.7, 'fixed']
            elif color==config['color2']['name']:
                I_e = [16.85, 'fixed']

            #I_e = [config['sersic_params']['I_e_guess'], 'fixed']
        else: I_e = [config['sersic_params']['I_e_guess'], config['sersic_params']['I_e_min'], config['sersic_params']['I_e_max']]

        if config['inject_artpop_model']:
            center = [config['xpos_inject'],config['ypos_inject']]
        else:
            center = [config['sersic_params']['xpos_guess'],config['sersic_params']['ypos_guess']]

        init_params = dict(PA=PA, n=n, ell=ell,r_e=r_e,I_e=I_e)
        dcent = config['sersic_params']['pos_err']

    elif fixedsersic is not None:
        fixedparams = fixedsersic.results
        for i in range(len(fixedparams)-1):

            if fixedparams[f'comp_{i+1}']['function'] == 'Sersic':
                PA = [fixedparams[f'comp_{i+1}']['PA'], 'fixed']
                n = [fixedparams[f'comp_{i+1}']['n'], 'fixed']
                ell = [fixedparams[f'comp_{i+1}']['ell'], 'fixed']
                r_e = [fixedparams[f'comp_{i+1}']['r_e'], 'fixed']
                I_e = [config['sersic_params']['I_e_guess'], config['sersic_params']['I_e_min'], config['sersic_params']['I_e_max']]

                init_params = dict(PA=PA, n=n, ell=ell,r_e=r_e,I_e=I_e)
                center = [fixedparams[f'comp_{i+1}']['X0'],fixedparams[f'comp_{i+1}']['Y0']]
                dcent = 0
                continue

    elif model == 'TiltedSkyPlane':
        I_0 = [config['tilted_plane_params']['I_0_guess'],config['tilted_plane_params']['I_0_min'],config['tilted_plane_params']['I_0_max']]
        m_x = [config['tilted_plane_params']['m_x_guess'],config['tilted_plane_params']['m_x_min'],config['tilted_plane_params']['m_x_max']]
        m_y = [config['tilted_plane_params']['m_y_guess'],config['tilted_plane_params']['m_y_min'],config['tilted_plane_params']['m_y_max']]
        init_params = dict(I_0=I_0, m_x=m_x, m_y=m_y)
        '''
        center = None
        dcent = None
        '''
        if config['inject_artpop_model']:
            center = [config['xpos_inject'],config['ypos_inject']]
            dcent=1000
        else:
            center = [config['sersic_params']['xpos_guess'],config['sersic_params']['ypos_guess']]
            dcent=1000

    elif model == 'FlatSky':
        I_sky = [config['flat_sky_params']['I_sky_guess'],config['flat_sky_params']['I_sky_min'],config['flat_sky_params']['I_sky_max']]
        init_params = dict(I_sky=I_sky)
        '''
        center = None
        dcent = None
        '''
        if config['inject_artpop_model']:
            center = [config['xpos_inject'],config['ypos_inject']]
            dcent=1000
        else:
            center = [config['sersic_params']['xpos_guess'],config['sersic_params']['ypos_guess']]
            dcent=1000

    elif model == 'Gaussian':
        PA = [config['gauss_params']['PA_guess'],config['gauss_params']['PA_min'],config['gauss_params']['PA_max']]
        ell = [config['gauss_params']['ell_guess'],config['gauss_params']['ell_min'],config['gauss_params']['ell_max']]
        I_0 = [config['gauss_params']['I_0_guess'],config['gauss_params']['I_0_min'],config['gauss_params']['I_0_max']]
        sigma = [config['gauss_params']['sigma_guess'],config['gauss_params']['sigma_min'],config['gauss_params']['sigma_max']]
        center = [config['gauss_params']['xpos_guess'],config['gauss_params']['ypos_guess']]
        init_params = OrderedDict([('PA', PA),('ell', ell),('I_0', I_0),('sigma', sigma)])
        dcent = config['gauss_params']['pos_err']

    return init_params, center, dcent

def run_imfit(img_fn, mask_fn, color_specific_info, config, model_funcs, options='', fixedsersic=None, viz=False, iter=None, fn_stub=None, alt_out_dir=None, alt_image_dir=None, glob_select=None):

    psf_fn = os.path.join(config['image_dir'],color_specific_info['psf'])
    color = color_specific_info['name']
    if alt_out_dir is None: out_dir = config['out_dir']
    else: out_dir = alt_out_dir
    if alt_image_dir is None: image_dir = config['image_dir']
    else: image_dir = alt_image_dir

    # Get model/resid filenames
    model_sum = ''
    for func in model_funcs: model_sum+=f'{func}_'
    if f'_{color}.fits' in img_fn:
        model_fn = os.path.join(out_dir,img_fn.replace(f'_{color}.fits',f'_{model_sum}model_{color}.fits'))
        resid_fn = os.path.join(out_dir,img_fn.replace(f'_{color}.fits',f'_{model_sum}resid_{color}.fits'))
    else:
        model_fn = os.path.join(out_dir,img_fn.replace(glob_select,f'{model_sum}model.fits'))
        resid_fn = os.path.join(out_dir,img_fn.replace(glob_select,f'{model_sum}resid.fits'))
    if iter is None:
        if f'_{color}.fits' in img_fn: bestfit_fn = os.path.join(out_dir,img_fn.replace(f'_{color}.fits',f'_{model_sum}bestfit-params_{color}.txt'))
        else: bestfit_fn = os.path.join(out_dir,img_fn.replace(glob_select,f'{model_sum}bestfit-params.txt'))
    else: bestfit_fn = os.path.join(out_dir,img_fn.replace(f'_{color}.fits',f'_{model_sum}bestfit-params_{fn_stub}_iter{iter}_{color}.txt'))
    config_fn = f'config_{model_sum}{color}.txt'

    params = []
    centers = []
    dcents = []
    for func in model_funcs:
        init_params, center, dcent = organize_initial_params(config, func, fixedsersic, color)
        params.append(init_params)
        centers.append(center)
        dcents.append(dcent)
        if dcent is not None: dcent_final = dcent

    model = pymfit.Model(funcs = model_funcs, params = params, centers = centers, dcent=dcent_final)

    fitter = pymfit.PymFitter(model,save_files=True)
    print(fitter.model.funcs)
    print(fitter.model)
    if config['variance_image'] is not None:
        print(mask_fn,psf_fn,bestfit_fn)
        fitter.run(os.path.join(image_dir,img_fn), var_fn=os.path.join(image_dir,config['variance_image']), mask_fn=mask_fn, psf_fn=psf_fn, out_fn=bestfit_fn, outdir=out_dir, config_fn=config_fn, save_model=True, save_residual=True, will_viz=True)
    else:
        fitter.run(os.path.join(image_dir,img_fn), mask_fn=mask_fn, psf_fn=psf_fn, out_fn=bestfit_fn, outdir=out_dir, config_fn=config_fn, save_model=True, save_residual=True, will_viz=True, options=options)

    io.rename_fits_file(os.path.join(image_dir,img_fn.replace('.fits','_model.fits')), model_fn, delete_old=True, overwrite=config['overwrite'])
    io.rename_fits_file(os.path.join(image_dir,img_fn.replace('.fits','_res.fits')), resid_fn, delete_old=True, overwrite=config['overwrite'])

    if len(model_funcs) > 1 and 'Sersic' in model_funcs:
        run_makeimage(bestfit_fn, psf_fn=psf_fn, ref_fn=os.path.join(image_dir,img_fn), output_root=model_fn.replace('.fits','_'))
        sersic_comp = model_funcs.index('Sersic')+1
        model_fn = model_fn.replace('.fits',f'_{sersic_comp}_Sersic.fits')

    if viz:
        if iter is None: fn = model_fn.replace('.fits','.png')
        else: fn = model_fn.replace('.fits',f'_iter{iter}.png')
        fitter.viz_results(show=False, save_fn=fn, dpi=200)

    return fitter, model_fn, resid_fn, bestfit_fn #ADDED BESTFIT_FN TO LIST OF RETURNS!!


def summarize_results(config, sersic_params1, sersic_params2=None):

    zpt1 = config['color1']['zpt'] - 2.5*np.log10(config['gain']/config['exposure_time'])
    params1 = {'I_e': sersic_params1['I_e'],
          'r_e': sersic_params1['r_e'],
          'n': sersic_params1['n'],
          'X0': sersic_params1['X0'],
          'Y0': sersic_params1['Y0'],
          'ell': sersic_params1['ell'],
          'PA': sersic_params1['PA']}

    sersic1 = pymfit.sersic.Sersic(params1, zpt=zpt1, pixscale=config['pixscale'])

    if sersic_params2 is not None:
        zpt2 = config['color2']['zpt'] - 2.5*np.log10(config['gain']/config['exposure_time'])
        params2 = {'I_e': sersic_params2['I_e'],
              'r_e': sersic_params2['r_e'],
              'n': sersic_params2['n'],
              'X0': sersic_params2['X0'],
              'Y0': sersic_params2['Y0'],
              'ell': sersic_params2['ell'],
              'PA': sersic_params2['PA']}

        sersic2 = pymfit.sersic.Sersic(params2, zpt=zpt2, pixscale=config['pixscale'])

        color = sersic2.m_tot - sersic1.m_tot

        # Re-calculate mags using color terms and extinction
        mag1_corrected = sersic1.m_tot - config['color1']['extinction'] - config['color1']['color_term']*color
        mag2_corrected = sersic2.m_tot - config['color2']['extinction'] - config['color2']['color_term']*color

        color_corrected = mag2_corrected - mag1_corrected

        return mag1_corrected , mag2_corrected , color_corrected

    return sersic1.m_tot - config['color1']['extinction']

def determine_imfit_comps(imfit_config):

    step_num = 1
    for step in imfit_config['imfit_steps']:
        if imfit_config['imfit_steps'][step]['color'] == imfit_config['color1']['name']:
            step1 = step_num
        elif imfit_config['imfit_steps'][step]['color'] == imfit_config['color2']['name']:
            step2 = step_num
        step_num+=1
    comp2 = misc.list_of_strings(imfit_config['imfit_steps'][f'step{step2}']['functions']).index('Sersic')+1
    comp1 = misc.list_of_strings(imfit_config['imfit_steps'][f'step{step1}']['functions']).index('Sersic')+1
    comp2 = f'comp_{comp2}'
    comp1 = f'comp_{comp1}'

    return comp1,comp2

def subtract_bright_star(filename, config, color, glob_select, use_cutout=False, cutout_size=(2001,2001)):

    # Read in imfit_config and csv
    with open(config['imfit_config'], 'r') as config_imfit:
        imfit_config = yaml.load(config_imfit, Loader=yaml.FullLoader)

    if config['custom_mask_csv_fn'] is not None:
        custom_mask_details = Table.from_pandas(pandas.read_csv(config['custom_mask_csv_fn']))
        custom_mask_details = custom_mask_details[np.where(custom_mask_details['filename'] == filename.split('/')[-1].replace(glob_select,config['glob_select']))]

    # Deal with directories, files, etc. - assumes flat variance image for now
    mkdir_if_needed(os.path.join(config['out_dir'],'star_subtracted'))
    file = fits.open(filename)[0]
    if use_cutout:
        var = fits.PrimaryHDU(np.zeros((cutout_size))+1)
    else:
        var = fits.PrimaryHDU(np.zeros((file.data.shape))+1)
    var.writeto(os.path.join(config['out_dir'],'star_subtracted',imfit_config['variance_image']),overwrite=True)

    # Create mask
    mask_fn = os.path.join(config['out_dir'],'star_subtracted', filename.replace(config['glob_select'],'mask_' + config['glob_select']))
    mask_kws = dict(out_fn=mask_fn, gal_pos=(custom_mask_details['x_star'][0],custom_mask_details['y_star'][0]), **imfit_config['masking_imfit_star'])
    mask = pymfit.make_mask(filename, **mask_kws)
    mask[np.where(file.data==0.)]=1

    # Create custom mask if necessary
    if config['custom_mask_csv_fn'] is not None:
        xpos = [] # xpos and ypos are in the FITS standard
        ypos = []
        radii = []
        for col in custom_mask_details.colnames:
            if col == 'filename' or 'star' in col: continue
            elif 'x_' in col: xpos.append(custom_mask_details[col][0])
            elif 'y_' in col: ypos.append(custom_mask_details[col][0])
            elif 'radius_' in col: radii.append(custom_mask_details[col][0])
            else: logger.warning('Unrecognized column name in ' + config['custom_mask_csv_fn'])

        # Add custom masks
        mask_orig = mask
        ii, jj = np.mgrid[:mask.shape[0], :mask.shape[1]]
        for idx in range(len(radii)):
            circ_mask = ((ii - ypos[idx])**2 + (jj - xpos[idx])**2) - radii[idx]**2 < 0
            mask = mask.astype(bool) | circ_mask
            mask = mask.astype(float)

        # Include star anti-mask
        star_mask = ((ii - custom_mask_details['y_star'][0])**2 + (jj - custom_mask_details['x_star'][0])**2) - custom_mask_details['radius_star'][0]**2 < 0
        mask[star_mask] = 0
        # Include diffraction spike masks, misc. masks
        #mask[0:file.data.shape[1],custom_mask_details['x_star'][0]-15:custom_mask_details['x_star'][0]+15]=1
        mask[np.where(file.data>config['saturation_limit'])]=1

        # Save mask
        mask_hdu = fits.PrimaryHDU(mask)
        mask_hdu.writeto(mask_fn, overwrite=True)

    # Create cutouts for imfit
    imfit_file = filename
    if use_cutout:
        imfit_file = misc.make_cutout(filename,(custom_mask_details['x_star'][0],custom_mask_details['y_star'][0]),cutout_size,cutout_fn=filename.replace(config['glob_select'],'cutout_'+config['glob_select']),force_shape=True)
        mask_file = misc.make_cutout(mask_fn,(custom_mask_details['x_star'][0],custom_mask_details['y_star'][0]),cutout_size,cutout_fn=mask_fn.replace(config['glob_select'],'cutout_'+config['glob_select']),force_shape=True)
        # Modify mask if necessary
        mask_file = fits.open(mask_fn.replace(config['glob_select'],'cutout_'+config['glob_select']))
        mask_file[0].data[np.where(imfit_file.data==0)] = 1
        mask_file.writeto(mask_fn.replace(config['glob_select'],'cutout_'+config['glob_select']), overwrite=True)
        # Save cutout filenames
        imfit_file = filename.replace(config['glob_select'],'cutout_'+config['glob_select'])
        mask_fn = mask_fn.replace(config['glob_select'],'cutout_'+config['glob_select'])


    # Do imfit stuff for star cutouts
    rdnoise = config['readnoise']
    gain = config['gain']
    options = f'--readnoise {rdnoise} --gain {gain}'
    color_options = {'name' : color, 'psf' : config[f'psf_fn_{color}']}

    for step in imfit_config['star_imfit']:
        band = imfit_config['star_imfit'][step]['color']
        funcs = misc.list_of_strings(imfit_config['star_imfit'][step]['functions'])
        if band != color : continue
        if 'Gaussian' in funcs:
            if use_cutout:
                imfit_config['gauss_params']['xpos_guess'] = int(cutout_size[0]/2)
                imfit_config['gauss_params']['ypos_guess'] = int(cutout_size[1]/2)
            else:
                imfit_config['gauss_params']['xpos_guess'] = custom_mask_details['x_star'][0]
                imfit_config['gauss_params']['ypos_guess'] = custom_mask_details['y_star'][0]
        # Run imfit
        results, model_fn, resid_fn, bf_fn = run_imfit(imfit_file, mask_fn, color_options, imfit_config, model_funcs=funcs, options=options, viz=True, alt_out_dir=os.path.join(config['out_dir'],'star_subtracted'), alt_image_dir=os.path.join(config['out_dir'],'star_subtracted'),glob_select=config['glob_select'])

    # Reconstruct model for whole image
    if use_cutout:
        # Change positions in bestfit file
        new_lines=[]
        funcs = []
        with open(os.path.join(config['out_dir'],'star_subtracted',bf_fn), 'r') as bf:

            for line in bf.readlines():
                if line[0] == 'X':
                    # Get cutout position
                    cutout_x = line.split('\t')[2].split(' ')[0]
                    # Get new whole image position
                    whole_im_x = custom_mask_details['x_star'][0] - (int(cutout_size[0]/2)-float(cutout_x))
                    new_line = line.replace(line.split('\t')[2].split(' ')[0],str(whole_im_x))
                    new_lines.append(new_line)

                elif line[0] == 'Y':
                    # Get cutout position
                    cutout_y = line.split('\t')[2].split(' ')[0]
                    # Get new whole image position
                    whole_im_y = custom_mask_details['y_star'][0] - (int(cutout_size[1]/2)-float(cutout_y))
                    new_line = line.replace(line.split('\t')[2].split(' ')[0],str(whole_im_y))
                    new_lines.append(new_line)

                elif line[0] == 'F':
                    funcs.append(line.split(' ')[-1].split('\n')[0])
                    if 'FlatSky' not in line: new_lines.append(line)
                    else: new_lines = new_lines[0:-2]

                elif 'sky' not in line: new_lines.append(line)

        os.rename(os.path.join(config['out_dir'],'star_subtracted',bf_fn),os.path.join(config['out_dir'],'star_subtracted',bf_fn.replace('.txt','_cutoutpositions.txt')))

        with open(os.path.join(config['out_dir'],'star_subtracted',bf_fn), 'w') as f:
            f.writelines(new_lines)

        # Run makeimage to recreate model for whole image
        run_makeimage(os.path.join(config['out_dir'],'star_subtracted',bf_fn) , psf_fn=os.path.join(config['out_dir'],config[f'psf_fn_{color}']), ref_fn=filename, output_root=None, out_fn=filename.replace('.fits','_model.fits'))

        # Subtract model from original image, save files
        file.data = file.data - fits.getdata(filename.replace('.fits','_model.fits'))
        file.writeto(filename,overwrite=True)

    return results, model_fn, resid_fn, bf_fn

def smooth_image():

    return
