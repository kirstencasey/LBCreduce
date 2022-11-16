import lbcred, os, yaml, random, glob, warnings, shutil
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from lbcred.log import logger
from lbcred.utils import misc
from lbcred.model import imfit
from astropy.wcs import WCS
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning


def make_cutout(original_img_fn, position, shape, ext = 0, cutout_fn=None):
    img = fits.open(original_img_fn)
    cutout = Cutout2D(img[ext].data, position, shape)
    img[ext].data = cutout.data
    if cutout_fn is not None:
        img.writeto(cutout_fn, overwrite=True)
    img.close()
    return cutout

imfit_config = '/users/PCON0003/osu10713/projects/LBCreduce/modeling-config_imfit.yml'
sbf_config = '/users/PCON0003/osu10713/projects/LBCreduce/modeling-config_sbf.yml'
register_config = '/users/PCON0003/osu10713/projects/LBCreduce/register_calibration-config.yml'
use_premade_masks =  True
num_artpop_models_to_test = 5 #0 # If 0 uses all models in directory; if this is 1 then just grabs the artpop model that's specified in the relevant config file
construct_artpop = False
select_premade_artpops = True
inject_artpop_in_registration_step = False # If this is False and select_premade_artpops is True, injects models in imfit step

back_models = ['SEsky']#,'polynomial','median']
register_calibrate_ims = True
run_imfit = True
run_sbf = True
xpos = [1815,1477] #[818,1041,1427,1516]
ypos = [3658,1586] #[1716,1189,1939,3364]
mask_bright_star_pos = [False,False]#,False,False,False]
cutout_size = (1051,1051)
num_positions = len(xpos)

with open(register_config, 'r') as filename:
    reg_config = yaml.load(filename, Loader=yaml.FullLoader)

with open(imfit_config, 'r') as filename:
    im_config = yaml.load(filename, Loader=yaml.FullLoader)
    orig_image_dir = im_config['image_dir']

if inject_artpop_in_registration_step:
    artpop_dir = reg_config['artpop_model_dir']
    reg_config['inject_artpop'] = True
    im_config['inject_artpop_model'] = False
else:
    artpop_dir = im_config['pregen_artpop_dir']
    reg_config['inject_artpop'] = False
    im_config['inject_artpop_model'] = True

artpop_models = glob.glob(os.path.join(artpop_dir,'*_r.fits'))
total_num_models = len(artpop_models)
im_config['use_src_counts'] = True

if num_artpop_models_to_test > total_num_models:
    print('WARNING: num_artpop_models_to_test >= total_num_models\n         Using all available artpop models.\n')
elif num_artpop_models_to_test == 0:
    num_artpop_models_to_test = total_num_models

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
true_mags_r = []
true_mags_b = []
true_sbfmags = []
true_dists = []
true_radii = []
true_ellip = []
true_n = []
true_pa = []
position_nums = []
artpop_ids = []
bckgnd_models = []

def run_registration_calibration_steps(reg_config, im_config, iter, back_model, first_run, artpop_ids=None, inject_artpop_in_registration_step=False):
    with open(register_config, 'r') as filename:
        reg_config = yaml.load(filename, Loader=yaml.FullLoader)

    if first_run and inject_artpop_in_registration_step:
        reg_config['construct_artpop'] = construct_artpop
    else: reg_config['construct_artpop'] = False

    if not inject_artpop_in_registration_step:
        reg_config['construct_artpop'] = False

    if num_artpop_models_to_test > 1 and inject_artpop_in_registration_step:
        reg_config['artpop_model_fn_r'] = artpop_ids[-1]+'_r.fits'
        reg_config['artpop_model_fn_b'] = artpop_ids[-1]+'_b.fits'

    # Change config based on desired background model
    reg_config['background_model'] = back_model
    reg_config['xpos_inject'] = xpos[iter]
    reg_config['ypos_inject'] = ypos[iter]

    print('\nRegistering images')

    back_id = f'{back_model}{iter+1}'
    if num_artpop_models_to_test > 1 and inject_artpop_in_registration_step: back_id+=artpop_ids[-1]

    # Register r-band images (sky_pos=[] unless inject_artpop_in_registration_step = True)
    orig_glob_select = reg_config['glob_select']
    _, _, src = lbcred.register_images(tmp_path='/tmp', bandpass='R', make_plots=True, config_fn=reg_config, back_fn_id=back_id)

    # Register b-band images (sky_pos=[] unless inject_artpop_in_registration_step = True)
    reg_config['glob_select'] = orig_glob_select
    reg_config, sky_pos, _ = lbcred.register_images(tmp_path='/tmp', bandpass='B', make_plots=True, config_fn=reg_config, back_fn_id=back_id)

    reg_config['image_dir'] = reg_config['out_dir']
    reg_config['glob_select'] = reg_config['glob_select'].replace('.fits','_reg.fits')

    print('\nCalibrating images')
    lbcred.calibrate_images(config=reg_config)

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
    if run_imfit:
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

first_run = True
for iter in range(num_positions):

    for back_model in back_models:

        logger.info(f'Working on background model: {back_model}')
        used_artpops = []

        if not inject_artpop_in_registration_step:

            reg_config, im_config, out_fn_r, out_fn_b, stack_fn_r  = run_registration_calibration_steps(reg_config, im_config, iter, back_model, first_run)
            first_run = False

        for model_num in range(num_artpop_models_to_test):
            position_nums.append(iter+1)
            bckgnd_models.append(back_model)

            if num_artpop_models_to_test == total_num_models:
                artpop_ids.append(artpop_models[model_num].split('/')[-1].split('_r.fits')[0])
            elif num_artpop_models_to_test > 1 :
                artpop_idx = random.randrange(0, total_num_models)
                while artpop_idx in used_artpops:
                    artpop_idx = random.randrange(0, total_num_models)
                used_artpops.append(artpop_idx)
                artpop_ids.append(artpop_models[artpop_idx].split('/')[-1].split('_r.fits')[0])
            else:
                artpop_ids.append(im_config['color1']['artpop_model_fn'].split('_r.fits')[0])

            if register_calibrate_ims and inject_artpop_in_registration_step:

                reg_config, im_config, _, _, stack_fn_r = run_registration_calibration_steps(reg_config, im_config, iter, back_model, first_run, artpop_ids, inject_artpop_in_registration_step=inject_artpop_in_registration_step)
                first_run = False

            elif register_calibrate_ims:
                im_config['color1']['artpop_model_fn'] = artpop_ids[-1]+'_r.fits'
                im_config['color2']['artpop_model_fn'] = artpop_ids[-1]+'_b.fits'
                im_config['inject_artpop_model'] = True
                shutil.copyfile(os.path.join(orig_image_dir,out_fn_r),os.path.join(orig_image_dir,out_fn_r.replace('r.fits',artpop_ids[-1]+'_r.fits')))
                shutil.copyfile(os.path.join(orig_image_dir,out_fn_b),os.path.join(orig_image_dir,out_fn_b.replace('b.fits',artpop_ids[-1]+'_b.fits')))
                im_config['color1']['image_fn'] = out_fn_r.replace('r.fits',artpop_ids[-1]+'_r.fits')
                im_config['color2']['image_fn'] = out_fn_b.replace('b.fits',artpop_ids[-1]+'_b.fits')
                im_config['image_dir'] = orig_image_dir

            # Run imfit
            if run_imfit:

                im_config['mask_bright_star'] = mask_bright_star_pos[iter]

                if select_premade_artpops :
                    artpop_hdr = fits.getheader(os.path.join(artpop_dir,artpop_ids[-1]+'_r.fits'),hdu_number=reg_config['ext'])
                    im_config['sersic_params']['PA_guess'] = float(artpop_hdr['THETA']) + 90.0
                    im_config['sersic_params']['n_guess'] = float(artpop_hdr['N'])
                    im_config['sersic_params']['ell_guess'] = float(artpop_hdr['ELLIP'])
                    im_config['sersic_params']['r_e_guess'] = misc.parsecs_to_pixels(float(artpop_hdr['R_EFF']), float(artpop_hdr['DISTANCE']) * 1e6, im_config['pixscale'])

                if use_premade_masks :
                    im_config['color1']['mask_fn'] = f'cutout{iter+1}_mask_r.fits'
                    im_config['color2']['mask_fn'] = f'cutout{iter+1}_mask_b.fits'

                print('Running imfit')
                bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, functions, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true = lbcred.modeling_imfit(config_filename=im_config)

                r_comp, b_comp = imfit.determine_imfit_comps(im_config)
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

            # Run sbf
            if run_sbf:

                with open(sbf_config, 'r') as filename:
                    sb_config = yaml.load(filename, Loader=yaml.FullLoader)

                sb_config['image_dir'] = im_config['out_dir']
                sb_config['out_dir'] = im_config['out_dir']
                sb_config['model_fn'] = model1_fn
                sb_config['color'] = color
                sb_config['resid_fn'] = resid1_fn
                sb_config['model_summary_fn'] = bestfit1_fn
                sb_config['stack_fn'] = stack_fn_r

                print('Running sbf')
                run_id = f'{back_model}{iter+1}_'
                if num_artpop_models_to_test > 1 : run_id+=artpop_ids[-1]
                sbf_mag, dist_a, dist_b = lbcred.modeling_sbf(config_filename=sb_config, imfit_functions=functions,run_id=run_id)

                measured_sbfmags.append(sbf_mag)
                measured_dists.append(dist_a)

            if inject_artpop_in_registration_step:
                artpop_fn = os.path.join(reg_config['artpop_model_dir'],reg_config['artpop_model_fn_r'])
            else:
                artpop_fn = os.path.join(im_config['pregen_artpop_dir'],artpop_ids[-1]+'_r.fits')

            artpop_hdr = fits.getheader(artpop_fn,hdu_number=0)
            true_dists.append(float(artpop_hdr['distance']))
            true_radii.append(float(artpop_hdr['r_eff']))
            true_ellip.append(float(artpop_hdr['ellip']))
            true_n.append(float(artpop_hdr['n']))
            true_pa.append(float(artpop_hdr['theta'])+90.)
            true_mags_r.append(float(artpop_hdr['Bessell_R_mag']))
            true_mags_b.append(float(artpop_hdr['Bessell_B_mag']))
            true_sbfmags.append(float(artpop_hdr['Bessell_R_sbfmag']))

if run_imfit: 
    arrs = [measured_mags_r,measured_mags_b,measured_sbfmags,measured_dists,measured_radii,measured_ellip,measured_n,measured_pa,measured_xpos,measured_ypos,measured_Ier,measured_Ieb,true_mags_r,true_mags_b,true_sbfmags,true_dists,true_radii,true_ellip,true_n,true_pa,position_nums,artpop_ids,bckgnd_models]
    params = ['measured_mags_r','measured_mags_b','measured_sbfmags','measured_dists','measured_radii','measured_ellip','measured_n','measured_pa','measured_xpos','measured_ypos','measured_Ier','measured_Ieb','true_mags_r','true_mags_b','true_sbfmags','true_dists','true_radii','true_ellip','true_n','true_pa','position_ids','artpop_ids','background_models']
    
    for param, arr in zip(params, arrs):
        nparr = np.asarray(arr)
        file = open(os.path.join(im_config['out_dir'],f'background_subtraction_test_results_{param}'), "wb")
        np.save(file, nparr)
        file.close
