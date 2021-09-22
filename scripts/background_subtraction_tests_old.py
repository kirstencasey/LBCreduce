import lbcred, os, yaml, random
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from lbcred.log import logger
from lbcred.utils import misc
from astropy.wcs import WCS
from astropy import units as u


def make_cutout(original_img_fn, position, shape, ext = 0, cutout_fn=None):
    img = fits.open(original_img_fn)
    cutout = Cutout2D(img[ext].data, position, shape)
    img[ext].data = cutout.data
    if cutout_fn is not None:
        img.writeto(cutout_fn, overwrite=True)
    img.close()
    return cutout

imfit_config = '/Users/kirstencasey/projects/LBCreduce/modeling-config_imfit.yml'
sbf_config = '/Users/kirstencasey/projects/LBCreduce/modeling-config_sbf.yml'
register_config = '/Users/kirstencasey/projects/LBCreduce/register_calibration-config.yml'
in_fn_b = 'lbcb.20191220.080715-chip2_backsub_M81blob.proc_reg.fits'
in_fn_r = 'lbcr.20191220.080711-chip2_backsub_M81blob.proc_reg.fits'
use_premade_masks = True
num_positions = 1
num_artpop_models_to_test = 50 # If this is 1 then it will just grab the artpop model that's specified in the register_calibration config file
total_num_models = 100
construct_artpop = False
select_premade_artpops = True
r_comp = 'comp_1'
b_comp = 'comp_1'

back_models = ['median']#,'polynomial','SEsky']
register_calibrate_ims = True
run_imfit = True
run_sbf = True
xpos = [818,1041,1427,1516]
ypos = [1716,1189,1939,3364]
mask_bright_star_pos = [True,False,False,False]

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
artpop_nums = []
models = []


for iter in range(num_positions):

    mask_bright_star = mask_bright_star_pos[iter]

    for model in back_models:

        logger.info(f'Working on background model: {model}')
        used_artpops = []

        for model_num in range(num_artpop_models_to_test):
            position_nums.append(iter+1)
            models.append(model)
            if num_artpop_models_to_test > 1 :
                artpop_num = random.randrange(0, total_num_models)
                while artpop_num in used_artpops:
                    artpop_num = random.randrange(0, total_num_models)
                used_artpops.append(artpop_num)
                artpop_nums.append(artpop_num)

            if register_calibrate_ims:

                with open(register_config, 'r') as filename:
                    reg_config = yaml.load(filename, Loader=yaml.FullLoader)

                if iter==0:
                    reg_config['construct_artpop'] = construct_artpop
                else:  reg_config['construct_artpop'] = False
                if num_artpop_models_to_test > 1 :
                    reg_config['artpop_model_fn_r'] = f'artpop_model_{artpop_num}_r.fits'
                    reg_config['artpop_model_fn_b'] = f'artpop_model_{artpop_num}_b.fits'

                # Change config based on desired background model
                reg_config['background_model'] = model
                reg_config['xpos_inject'] = xpos[iter]
                reg_config['ypos_inject'] = ypos[iter]

                print('Registering images')

                back_id = f'{model}{iter+1}'
                if num_artpop_models_to_test > 1 : back_id+=f'_artpop{artpop_num}'

                # Register r-band images
                reg_config, sky_pos, src = lbcred.register_images(tmp_path='/tmp', bandpass='R', make_plots=True, config_fn=reg_config, back_fn_id=back_id)

                # Register b-band images
                reg_config, sky_pos, _ = lbcred.register_images(tmp_path='/tmp', bandpass='B', make_plots=True, config_fn=reg_config, back_fn_id=back_id)

                reg_config['image_dir'] = reg_config['out_dir']
                if reg_config['subtract_background']: reg_config['glob_select'] = '_backsub_' + reg_config['glob_select'].replace('.fits','_reg.fits')
                else: reg_config['glob_select'] = reg_config['glob_select'].replace('.fits','_reg.fits')

                print('Calibrating images')
                lbcred.calibrate_images(config=reg_config, stack=False)

                out_fn_b = f'lbcb.20191220.080715-chip2_backsub_M81blob.proc_reg_cutout_{model}{iter+1}_b.fits'
                out_fn_r = f'lbcr.20191220.080711-chip2_backsub_M81blob.proc_reg_cutout_{model}{iter+1}_r.fits'
                if num_artpop_models_to_test > 1 :
                    out_fn_b = out_fn_b.replace('_b.fits',f'_artpop{artpop_num}_b.fits')
                    out_fn_r = out_fn_r.replace('_r.fits',f'_artpop{artpop_num}_r.fits')

                w = WCS(fits.open(os.path.join(reg_config['image_dir'],'sci',in_fn_b))[0].header)
                x_cutout, y_cutout = w.world_to_pixel(sky_pos[0])

                cutout_b = make_cutout(os.path.join(reg_config['image_dir'],'sci',in_fn_b), (x_cutout,y_cutout), (1051,1051), ext = 0, cutout_fn=os.path.join(reg_config['image_dir'],'sci','imfit_sbf',out_fn_b))
                cutout_r = make_cutout(os.path.join(reg_config['image_dir'],'sci',in_fn_r), (x_cutout,y_cutout), (1051,1051), ext = 0, cutout_fn=os.path.join(reg_config['image_dir'],'sci','imfit_sbf',out_fn_r))

                if src is not None: artpop_src = src

            # Run imfit
            if run_imfit:

                with open(imfit_config, 'r') as filename:
                    im_config = yaml.load(filename, Loader=yaml.FullLoader)

                im_config['image_dir'] = os.path.join(reg_config['out_dir'],'sci','imfit_sbf')
                im_config['out_dir'] = os.path.join(reg_config['out_dir'],'sci','imfit_sbf')
                im_config['color1']['image_fn'] = out_fn_r
                im_config['color2']['image_fn'] = out_fn_b
                im_config['mask_bright_star'] = mask_bright_star

                if select_premade_artpops :
                    artpop_hdr = fits.getheader(os.path.join(reg_config['artpop_model_dir'],reg_config['artpop_model_fn_r']),hdu_number=reg_config['ext'])
                    im_config['sersic_params']['PA_guess'] = float(artpop_hdr['THETA'])
                    im_config['sersic_params']['n_guess'] = float(artpop_hdr['N'])
                    im_config['sersic_params']['ell_guess'] = float(artpop_hdr['ELLIP'])
                    im_config['sersic_params']['r_e_guess'] = misc.parsecs_to_pixels(float(artpop_hdr['R_EFF']), float(artpop_hdr['DISTANCE']) * 1e6, im_config['pixscale'])

                if use_premade_masks :
                    im_config['color1']['mask_fn'] = f'mask{iter+1}_r.fits'
                    im_config['color2']['mask_fn'] = f'mask{iter+1}_b.fits'

                print('Running imfit')
                bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, functions, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true = lbcred.modeling_imfit(config_filename=im_config)

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

                print('Running sbf')
                run_id = f'{model}{iter+1}'
                if num_artpop_models_to_test > 1 : run_id+=f'_artpop{artpop_num}'
                sbf_mag, dist_a, dist_b = lbcred.modeling_sbf(config_filename=sb_config, imfit_functions=functions,run_id=run_id)

                measured_sbfmags.append(sbf_mag)
                measured_dists.append(dist_a)

            artpop_hdr = fits.getheader(os.path.join(reg_config['artpop_model_dir'],reg_config['artpop_model_fn_r']),hdu_number=reg_config['ext'])
            true_dists.append(float(artpop_hdr['distance']))
            true_radii.append(float(artpop_hdr['r_eff']))
            true_ellip.append(float(artpop_hdr['ellip']))
            true_n.append(float(artpop_hdr['n']))
            true_pa.append(float(artpop_hdr['theta']))
            true_mags_r.append(float(artpop_hdr['Bessell_R_mag']))
            true_mags_b.append(float(artpop_hdr['Bessell_B_mag']))
            true_sbfmags.append(float(artpop_hdr['Bessell_R_sbfmag']))

arrs = [measured_mags_r,measured_mags_b,measured_sbfmags,measured_dists,measured_radii,measured_ellip,measured_n,measured_pa,measured_xpos,measured_ypos,measured_Ier,measured_Ieb,true_mags_r,true_mags_b,true_sbfmags,true_dists,true_radii,true_ellip,true_n,true_pa,position_nums,artpop_nums,models]
params = ['measured_mags_r','measured_mags_b','measured_sbfmags','measured_dists','measured_radii','measured_ellip','measured_n','measured_pa','measured_xpos','measured_ypos','measured_Ier','measured_Ieb','true_mags_r','true_mags_b','true_sbfmags','true_dists','true_radii','true_ellip','true_n','true_pa','position_ids','artpop_ids','background_models']

for param, arr in zip(params, arrs):
    nparr = np.asarray(arr)
    file = open(os.path.join(reg_config['image_dir'],'sci','imfit_sbf',f'background_subtraction_test_results_{param}'), "wb")
    np.save(file, nparr)
    file.close
