import lbcred, os, yaml, random, glob, warnings, shutil, pandas
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from lbcred.log import logger
from lbcred.utils import misc, io
from lbcred.model import imfit
from astropy.wcs import WCS
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning
from pathlib import Path
from astropy.table import Table
from lbcred import sbf_functions


def make_cutout(original_img_fn, position, shape, ext = 0, cutout_fn=None):
    img = fits.open(original_img_fn)
    cutout = Cutout2D(img[ext].data, position, shape)
    img[ext].data = cutout.data
    if cutout_fn is not None:
        img.writeto(cutout_fn, overwrite=True)
    img.close()
    return cutout
    
def check_full_glob_select(config,ftype='sci'):
    
    glob_select = config['glob_select']
    chip_num = config['chip_num']
    
    # Get files
    fn_example = glob.glob(os.path.join(config['image_dir'],ftype,f'lbc*{glob_select}*'))[0]
    
    fn_glob = fn_example.split(f'-chip{chip_num}_')[-1]
    
    if fn_glob != glob_select: logger.warning(f'Changing glob_select from {glob_select} to {fn_glob}.')
    
    return fn_glob


def run_registration_calibration_steps(reg_config, register_config, im_config, iter, back_model, first_run, xpos, ypos, cutout_size, orig_image_dir, artpop_ids=None, inject_artpop_in_registration_step=False, num_artpop_models_to_test=0, use_premade_masks=False, construct_artpop=False,overwrite_out_dir=None,psf_num=None,artpop_dist=None):

    with open(register_config, 'r') as filename:
        reg_config = yaml.load(filename, Loader=yaml.FullLoader)
        if overwrite_out_dir != None: reg_config['out_dir'] = overwrite_out_dir

    if first_run and inject_artpop_in_registration_step:
        reg_config['construct_artpop'] = construct_artpop
    else: reg_config['construct_artpop'] = False

    if not inject_artpop_in_registration_step:
        reg_config['construct_artpop'] = False

    if num_artpop_models_to_test > 1 and inject_artpop_in_registration_step:
        reg_config['artpop_model_fn_r'] = artpop_ids[-1]+'_r.fits'
        reg_config['artpop_model_fn_b'] = artpop_ids[-1]+'_b.fits'

    # Change config based on desired background model, etc.
    reg_config['background_model'] = back_model
    reg_config['xpos_inject'] = xpos[iter]
    reg_config['ypos_inject'] = ypos[iter]
    if psf_num != None: 
        reg_config['psf_stars_idx'] = psf_num
        reg_config['choose_from_csv'] = True
        reg_config['update_psf_fns'] = True

    print('Registering images')

    back_id = f'{back_model}{iter+1}'
    if psf_num != None: psf_id = f'_psf{psf_num}'
    else: psf_id = ''
    if artpop_dist != None: art_id = f'_artpop{artpop_dist}'
    else: art_id = ''
    if num_artpop_models_to_test > 1 and inject_artpop_in_registration_step: back_id+=artpop_ids[-1]
    print('HERE WE GOOOOOOOO!!!!!!!')
    # Check for temp paths
    if not os.path.isdir(f'/fs/scratch/PCON0003/osu10713/temp/temp_{back_id}{psf_id}{art_id}'):
        os.mkdir(f'/fs/scratch/PCON0003/osu10713/temp/temp_{back_id}{psf_id}{art_id}')
        
    # Register r-band images (sky_pos=[] unless inject_artpop_in_registration_step = True)
    orig_glob_select = reg_config['glob_select']
    _, _, src = lbcred.register_images(tmp_path=f'/fs/scratch/PCON0003/osu10713/temp/temp_{back_id}{psf_id}{art_id}', bandpass='R', make_plots=True, config_fn=reg_config, back_fn_id=back_id)

    # Register b-band images (sky_pos=[] unless inject_artpop_in_registration_step = True)
    reg_config['glob_select'] = orig_glob_select
    reg_config, sky_pos, _ = lbcred.register_images(tmp_path=f'/fs/scratch/PCON0003/osu10713/temp/temp_{back_id}{psf_id}{art_id}', bandpass='B', make_plots=True, config_fn=reg_config, back_fn_id=back_id)

    reg_config['image_dir'] = reg_config['out_dir']
    reg_config['glob_select'] = reg_config['glob_select'].replace('.fits','_reg.fits')
    
    reg_config['glob_select'] = check_full_glob_select(reg_config)

    print('Calibrating images')
    psf_fn_r, psf_fn_b = lbcred.calibrate_images(config=reg_config)
    reg_config['psf_fn_r'] = psf_fn_r
    reg_config['psf_fn_b'] = psf_fn_b
    im_config['color1']['psf'] = psf_fn_r
    im_config['color2']['psf'] = psf_fn_b

    if reg_config['stack_images']:
        out_fn_r = reg_config['final_stack_fn_r'].replace('_r.fits',f'_cutout{iter+1}_r.fits')
        out_fn_b = reg_config['final_stack_fn_b'].replace('_b.fits',f'_cutout{iter+1}_b.fits')
        in_fn_r = reg_config['final_stack_fn_r']
        in_fn_b = reg_config['final_stack_fn_b']
    else:
        out_fn_r = im_config['color1']['image_fn'].replace('_r.fits',f'_cutout{iter+1}_r.fits')
        out_fn_b = im_config['color2']['image_fn'].replace('_b.fits',f'_cutout{iter+1}_b.fits')
        in_fn_r = im_config['color1']['image_fn']
        in_fn_b = im_config['color2']['image_fn']

    if num_artpop_models_to_test > 1 and inject_artpop_in_registration_step:
        out_fn_b = out_fn_b.replace('b.fits',artpop_ids[-1]+'_b.fits')
        out_fn_r = out_fn_r.replace('r.fits',artpop_ids[-1]+'_r.fits')

    im_config['color1']['image_fn'] = out_fn_r
    im_config['color2']['image_fn'] = out_fn_b

    if inject_artpop_in_registration_step:
        w = WCS(fits.open(os.path.join(reg_config['image_dir'],'sci',in_fn_b))[0].header)
        x_cutout, y_cutout = w.world_to_pixel(sky_pos[0])
    else: x_cutout, y_cutout = xpos[iter], ypos[iter]

    cutout_b = make_cutout(os.path.join(reg_config['image_dir'],'sci',in_fn_b), (x_cutout,y_cutout), cutout_size, ext = 0, cutout_fn=os.path.join(reg_config['image_dir'],'sci',out_fn_b))
    cutout_r = make_cutout(os.path.join(reg_config['image_dir'],'sci',in_fn_r), (x_cutout,y_cutout), cutout_size, ext = 0, cutout_fn=os.path.join(reg_config['image_dir'],'sci',out_fn_r))

    # Put everything in the directory where imfit will be run if necessary
    if im_config['run_imfit']:
        imfit_dir_exists = os.path.isdir(im_config['out_dir'])
        if imfit_dir_exists and first_run:
            warnings.warn('The imfit output directory already exists and will be overwritten.', AstropyUserWarning)
            shutil.rmtree(im_config['out_dir'])
            os.mkdir(im_config['out_dir'])
        elif not imfit_dir_exists:
            os.mkdir(im_config['out_dir'])
        # Copy stuff into imfit dir
        shutil.copyfile(os.path.join(reg_config['image_dir'],'sci',out_fn_r),os.path.join(orig_image_dir,out_fn_r))
        shutil.copyfile(os.path.join(reg_config['image_dir'],'sci',out_fn_b),os.path.join(orig_image_dir,out_fn_b))
        shutil.copyfile(os.path.join(reg_config['image_dir'],'sci',reg_config['psf_fn_r']),os.path.join(im_config['out_dir'],im_config['color1']['psf']))
        shutil.copyfile(os.path.join(reg_config['image_dir'],'sci',reg_config['psf_fn_b']),os.path.join(im_config['out_dir'],im_config['color2']['psf']))
        if first_run and im_config['variance_image'] is not None:
            warnings.warn('Assuming a flat variance image.', AstropyUserWarning)
            var_im = fits.PrimaryHDU(np.zeros(cutout_size)+1)
            var_im.writeto(os.path.join(im_config['out_dir'],im_config['variance_image']),overwrite=True)
        if use_premade_masks and first_run:
            all_masks = glob.glob(os.path.join(im_config['image_dir'],'*mask*'))
            for mask in all_masks:
                shutil.copyfile(mask,os.path.join(im_config['out_dir'],mask.split('/')[-1]))

    return reg_config, im_config, out_fn_r, out_fn_b, os.path.join(reg_config['image_dir'],'sci',in_fn_r)


#####################################################################################

def run_background_subtraction_tests(back_models,overwrite_out_dir=None,psf_num=None):
    
    imfit_config = '/users/PCON0003/osu10713/NGC672dwC/modeling-config_imfit.yml'
    sbf_config = '/users/PCON0003/osu10713/NGC672dwC/modeling-config_sbf.yml'
    register_config = '/users/PCON0003/osu10713/NGC672dwC/register_calibration-config.yml'
    use_premade_masks =  True
    construct_artpop = False
    select_premade_artpops = False
    inject_artpop_in_registration_step = False # If this is False and select_premade_artpops is True, injects models in imfit step
    
    use_real_zpt_colors = True
    register_calibrate_ims = True
    run_imfit = True
    run_sbf = True
    xpos = [1503] 
    ypos = [2749]
    mask_bright_star_pos = [False]
    size_pix = 1051 
    cutout_size = (size_pix,size_pix)  
    flat_variance = True

    with open(register_config, 'r') as filename:
        reg_config = yaml.load(filename, Loader=yaml.FullLoader)
        if overwrite_out_dir != None: reg_config['out_dir'] = overwrite_out_dir
        orig_reg_out_dir = reg_config['out_dir']
    
    with open(imfit_config, 'r') as filename: 
        im_config = yaml.load(filename, Loader=yaml.FullLoader)
        if overwrite_out_dir != None: 
            im_config['image_dir'] = overwrite_out_dir
            im_config['out_dir'] = overwrite_out_dir
        orig_image_dir = im_config['image_dir']
        if run_imfit: 
            im_config['run_imfit'] = True
            im_config['variance_image'] = f'flat_variance_image_{size_pix}.fits'
            if flat_variance and not os.path.exists(os.path.join(im_config['image_dir'],f'flat_variance_image_{size_pix}.fits')):
                flat_var = fits.PrimaryHDU(data=np.zeros(cutout_size) + 1.)
                flat_var.writeto(os.path.join(im_config['image_dir'],f'flat_variance_image_{size_pix}.fits'))
        else: im_config['run_imfit'] = False

    im_config['use_src_counts'] = True
    
    # Things for plots
    measured_mags_r = []
    measured_mags_b = []
    measured_sbfmags = []
    measured_dists = []
    measured_radii = []
    measured_ellip = []
    measured_n = []
    measured_pa = []
    measured_xpos = []
    measured_ypos = []
    measured_Ier = []
    measured_Ieb = []
    bckgnd_models = []
    uncertainty_sbf_mags = []
    
    measured_radii_b = []
    measured_ellip_b = []
    measured_n_b = []
    measured_pa_b = []
    iter=0
    
    
    for back_model in back_models:

        logger.info(f'Working on background model: {back_model}')
        bckgnd_models.append(back_model)
        im_config['out_dir'] = os.path.join(orig_image_dir,back_model)
    
        if register_calibrate_ims:

            reg_config['out_dir'] = os.path.join(orig_reg_out_dir,back_model)
            im_config['image_dir'] = os.path.join(orig_image_dir,back_model)
           
         
            reg_config, im_config, out_fn_r, out_fn_b, stack_fn_r   = run_registration_calibration_steps(reg_config, register_config, im_config, iter, back_model, True, xpos, ypos, cutout_size, orig_image_dir, use_premade_masks=use_premade_masks,overwrite_out_dir=overwrite_out_dir,psf_num=psf_num)
            mag_info1 = fits.getheader(os.path.join(reg_config['out_dir'],'sci',reg_config['final_stack_fn_r']))
            mag_info2 = fits.getheader(os.path.join(reg_config['out_dir'],'sci',reg_config['final_stack_fn_b']))
            im_config['color1']['zpt'] = mag_info1['STACK_ZPT']
            im_config['color2']['zpt'] = mag_info2['STACK_ZPT']
            im_config['color1']['color_term'] = mag_info1['STACK_COLORTERM_R']
            im_config['color2']['color_term'] = mag_info1['STACK_COLORTERM_B']
            im_config['image_dir'] = orig_image_dir
        else: stack_fn_r = os.path.join(im_config['image_dir'],'sci',reg_config['final_stack_fn_r'])

        # Run imfit
        if run_imfit:

            im_config['mask_bright_star'] = mask_bright_star_pos[iter]

            if use_premade_masks :

                # Make sure the masks are in the correct directory
                if not os.path.exists(os.path.join(im_config['image_dir'],back_model,im_config['color1']['mask_fn'])):
                    shutil.copyfile(os.path.join(orig_image_dir,im_config['color1']['mask_fn']),os.path.join(orig_image_dir,back_model,im_config['color1']['mask_fn']))
                    shutil.copyfile(os.path.join(orig_image_dir,im_config['color2']['mask_fn']),os.path.join(orig_image_dir,back_model,im_config['color2']['mask_fn']))
                    
            print('Running imfit')
            bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, functions, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true, orig_color = lbcred.modeling_imfit(config_filename=im_config)
            print('\n\nmodel1_fn: ',model1_fn)
            if len(np.where(np.asarray(functions)=='Sersic')[0]) > 1: imfit_bestfit = bestfit1
            else: imfit_bestfit = None
            r_comp, b_comp = imfit.determine_imfit_comps(im_config,imfit_bestfit,cutout_size)
            print('\n\nBESTFIT R_COMP')
            print(bestfit1[r_comp])

            measured_mags_r.append(mag1)
            measured_mags_b.append(mag2)
            measured_radii.append(bestfit1[r_comp]['r_e'])
            measured_ellip.append(bestfit1[r_comp]['ell'])
            measured_n.append(bestfit1[r_comp]['n'])
            measured_pa.append(bestfit1[r_comp]['PA'])
            measured_xpos.append(bestfit1[r_comp]['X0'])
            measured_ypos.append(bestfit1[r_comp]['Y0'])
            measured_Ier.append(bestfit1[r_comp]['I_e'])
            measured_Ieb.append(bestfit2[b_comp]['I_e'])
            
            measured_radii_b.append(bestfit2[b_comp]['r_e'])
            measured_ellip_b.append(bestfit2[b_comp]['ell'])
            measured_n_b.append(bestfit2[b_comp]['n'])
            measured_pa_b.append(bestfit2[b_comp]['PA'])
        
        elif not os.path.isdir(im_config['out_dir']): os.mkdir(im_config['out_dir'])

        # Run sbf
        if run_sbf:

            with open(sbf_config, 'r') as filename:
                sb_config = yaml.load(filename, Loader=yaml.FullLoader)
                if overwrite_out_dir != None: 
                    sb_config['image_dir'] = overwrite_out_dir
                    sb_config['out_dir'] = overwrite_out_dir

            sb_config['model_summary_fn'] = bestfit1_fn 
            
            sb_config['image_dir'] = im_config['out_dir']
            sb_config['out_dir'] = im_config['out_dir']
            sb_config['model_fn'] = model1_fn
            sb_config['color'] = orig_color
            sb_config['resid_fn'] = resid1_fn
            sb_config['stack_fn'] = stack_fn_r
            sb_config['psf'] = im_config['color1']['psf']
            sb_config['zpt'] = im_config['color1']['zpt']
            sb_config['color_term'] = im_config['color1']['color_term']

            print('Running sbf')
            run_id = f'{back_model}{iter+1}_'
            sbf_mag, dist_a, dist_b, uncert_mag = lbcred.modeling_sbf(config_filename=sb_config, imfit_functions=functions,run_id=run_id, plot_blank_fields=True)

            measured_sbfmags.append(sbf_mag)
            measured_dists.append(dist_a)
            uncertainty_sbf_mags.append(uncert_mag)
            
    arrs = [measured_mags_r,measured_mags_b,measured_sbfmags,measured_dists,measured_radii,measured_ellip,measured_n,measured_pa,measured_xpos,measured_ypos,measured_Ier,measured_Ieb,bckgnd_models,uncertainty_sbf_mags,measured_radii_b,measured_ellip_b,measured_n_b,measured_pa_b]
    params = ['measured_mags_r','measured_mags_b','measured_sbfmags','measured_dists','measured_radii','measured_ellip','measured_n','measured_pa','measured_xpos','measured_ypos','measured_Ier','measured_Ieb','background_models','sbf_mag_uncertainty','measured_radii_b','measured_ellip_b','measured_n_b','measured_pa_b']

    for param, arr in zip(params, arrs):
        nparr = np.asarray(arr)
        file = open(os.path.join(im_config['out_dir'],f'background_subtraction_test_results_{param}'), "wb")
        np.save(file, nparr)
        file.close
  
#####################################################################################

# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Model galaxy using imfit. Measure SBF distance.')
    parser.add_argument('--back_models', type=str)
    parser.add_argument('--out_dir', type=str)
    parser.add_argument('--psf_num', type=int)
    args = parser.parse_args()
    
    back_models = misc.list_of_strings(args.back_models)
    #out_dir = '/fs/scratch/PCON0003/osu10713/NGC672dwC/prelim_testing'
    run_background_subtraction_tests(back_models=back_models,overwrite_out_dir=args.out_dir,psf_num=args.psf_num)
    
    
    
    
    
