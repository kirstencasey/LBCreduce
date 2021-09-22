import lbcred, os, yaml, pymfit, glob
from astropy.io import fits
from astropy.nddata import Cutout2D
import numpy as np
import matplotlib.pyplot as plt
import random
from lbcred.model import artpop_functions
from lbcred.utils import misc
from regions import read_ds9


def make_cutout(original_img_fn, position, shape, ext = 0, cutout_fn=None):
    img = fits.open(original_img_fn)
    cutout = Cutout2D(img[ext].data, position, shape)
    img[ext].data = cutout.data
    if cutout_fn is not None:
        img.writeto(cutout_fn, overwrite=True)
    img.close()
    return cutout

# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Model galaxy using imfit. Measure SBF distance.')
    default_imfit_config = os.path.join(lbcred.project_dir, 'modeling-config_imfit.yml')
    default_sbf_config = os.path.join(lbcred.project_dir, 'modeling-config_sbf.yml')
    parser.add_argument('--imfit_config', type=str, default=default_imfit_config, help='path of the .yml config file ')
    parser.add_argument('--sbf_config', type=str, default=default_sbf_config, help='path of the .yml config file ')
    parser.add_argument('--num-iters', type=int, default=10, help='number of positions to inject models')
    parser.add_argument('--random_artpops', type=bool, default=True, help='inject randomly generated artpop models')
    parser.add_argument('--random_sersics', type=bool, default=False, help='inject randonmly generated smooth sersic models')
    parser.add_argument('--use_only_one_artpop', type=bool, default=False, help='use the same artpop model for each iteration')
    parser.add_argument('--run_sbf', type=bool, default=True, help='do sbf measurement')
    args = parser.parse_args()

    if args.random_artpops and args.random_sersics:
        print('Choose either artpop models or sersic models, not both. Assuming sersic models.\n')
        args.random_artpops = False

    fn_stub = 'real_image_tests_fixed_artpop_fix_sersic_on_b_autogenmask'
    use_skymodel = False

    ## The following is ignored if use_pregen_artpops is True
    # For artpop models:
    vary_age = 9.9 #(9.9,10.15)
    vary_feh = -1.0 #(-1.2,-1.0)
    vary_mass = 14900000.0  #noise:8500000
    vary_dist = 3.7 # Mpc
    vary_radius = 900 #(800,900) #350 # orig:894 noise:600 # This should be in parsecs!!
    # For both artpop and sersic models:
    vary_n = 0.5 #(0.3,0.7) #1.0 #0.52078
    vary_theta = 0. #(0.,180) #0.
    vary_ellip = 0. #(0.,0.4) #0.
    # For sersic models:
    vary_Ie_b = (100,300)
    vary_Ie_r = (600,800)
    if args.random_sersics:
        vary_radius = (200,300) # This should be in pixels!!

    # Open and read imfit config file
    with open(args.imfit_config, 'r') as filename:
        options = yaml.load(filename, Loader=yaml.FullLoader)

    options['inject_artpop_model'] = True

    # Get imfit r- and b-band components given imfit steps in config file
    r_comp, b_comp = imfit.determine_imfit_comps(options)

    # Open and read sbf config file
    with open(args.sbf_config, 'r') as filename:
        options_sbf = yaml.load(filename, Loader=yaml.FullLoader)

    full_star_mask_fn = os.path.join(options['out_dir'],'full_star_mask.fits')
    star_mask_fn = os.path.join(options['out_dir'],'star_mask.fits')
    skymodel_b_fn = os.path.join(options['out_dir'],'coadd_skymodel_b.fits')
    skymodel_r_fn = os.path.join(options['out_dir'],'coadd_skymodel_r.fits')
    skymodel_cutout_r_fn = os.path.join(options['out_dir'],'coadd_skymodel_cutout_r.fits')
    skymodel_cutout_b_fn = os.path.join(options['out_dir'],'coadd_skymodel_cutout_b.fits')
    regions_fn = os.path.join(options['out_dir'],'inject_positions_careful.reg')

    orig_fn_r = options['color1']['image_fn']
    orig_fn_b = options['color2']['image_fn']

    radii_recovered = []
    pa_recovered = []
    ell_recovered = []
    n_recovered = []
    I_e_r_recovered = []
    mags_r_recovered = []
    xpos_recovered = []
    ypos_recovered = []
    I_e_b_recovered = []
    mags_b_recovered = []
    sbf_mags_recovered = []
    sbf_mags_true = []
    mags_r_true = []
    mags_b_true = []
    sbf_dist_a_recovered = []
    sbf_dist_b_recovered = []
    ages = []
    fehs = []
    artpop_radii = []
    sersic_radii = []
    n_true = []
    thetas = []
    ells_true = []
    masses = []
    I_e_r_true = []
    I_e_b_true = []
    artpop_dists = []
    used_artpops = []
    artpop_nums = []


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    iter = 0
    skymodel_cutouts = None

    # Run imfit on galaxy at various positions
    while iter < args.num_iters:

        if options['use_pregen_artpops']:

            if options['run_artpop']:
                print('\nChoose either \'use_pregen_artpops\' or \'run_artpop\', not both. Assuming \'use_pregen_artpops\'.\n')
                options['run_artpop'] = False
            model_type = 'artpop'
            total_num_models = len(glob.glob(os.path.join(options['pregen_artpop_dir'],'*')))/2

            artpop_num = random.randrange(0, total_num_models)

            while artpop_num in used_artpops:
                artpop_num = random.randrange(0, total_num_models)
            used_artpops.append(artpop_num)
            artpop_nums.append(artpop_num)

            options['color1']['artpop_model_fn'] = f'artpop_model_{artpop_num}_r.fits'
            options['color2']['artpop_model_fn'] = f'artpop_model_{artpop_num}_b.fits'

            artpop_hdr = fits.getheader(os.path.join(options['pregen_artpop_dir'],options['color1']['artpop_model_fn']))
            options['color1']['artpop_mag'] = artpop_hdr['Bessell_R_mag']
            options['color2']['artpop_mag'] = artpop_hdr['Bessell_B_mag']
            options['color1']['artpop_sbfmag'] = artpop_hdr['Bessell_R_sbfmag']
            ages.append(artpop_hdr['LOG_AGE'])
            fehs.append(artpop_hdr['FEH'])
            artpop_radii.append(artpop_hdr['R_EFF'])
            n_true.append(artpop_hdr['N'])
            thetas.append(artpop_hdr['THETA'])
            ells_true.append(artpop_hdr['ELLIP'])
            masses.append(artpop_hdr['total_mass'])
            artpop_dists.append(artpop_hdr['DISTANCE'])

        elif args.random_artpops:
            # Create random artpop model
            name_stub = 'artpop_model'
            model_type = 'artpop'
            if type(vary_age) is tuple:
                options['artpop_model_params']['log_age'] = random.uniform(vary_age[0],vary_age[1])
            else: options['artpop_model_params']['log_age'] = vary_age
            age = options['artpop_model_params']['log_age']
            print(f'ArtPop Age: {round(age,4)}')
            ages.append(age)
            name_stub += f'_age{round(age,4)}'

            if type(vary_feh) is tuple:
                options['artpop_model_params']['feh'] = random.uniform(vary_feh[0],vary_feh[1])
            else: options['artpop_model_params']['feh'] = vary_feh
            feh = options['artpop_model_params']['feh']
            print(f'ArtPop FeH: {round(feh,4)}')
            fehs.append(feh)
            name_stub += f'_feh{round(feh,4)}'

            if type(vary_radius) is tuple:
                options['artpop_model_params']['r_eff'] = random.uniform(vary_radius[0],vary_radius[1])
            else: options['artpop_model_params']['r_eff'] = vary_radius
            radius = options['artpop_model_params']['r_eff']
            print(f'ArtPop Radius: {round(radius,2)} parsecs')
            artpop_radii.append(radius)
            name_stub += f'_radius{round(radius,2)}'

            if type(vary_n) is tuple:
                options['artpop_model_params']['n'] = random.uniform(vary_n[0],vary_n[1])
            else: options['artpop_model_params']['n'] = vary_n
            n = options['artpop_model_params']['n']
            print(f'ArtPop Sersic Index: {round(n,4)}')
            n_true.append(n)
            name_stub += f'_n{round(n,4)}'

            if type(vary_theta) is tuple:
                options['artpop_model_params']['theta'] = random.uniform(vary_theta[0],vary_theta[1])
            else: options['artpop_model_params']['theta'] = vary_theta
            theta =  options['artpop_model_params']['theta']
            print(f'ArtPop PA: {round(theta,4)}')
            thetas.append(theta)
            name_stub += f'_theta{round(theta,4)}'

            if type(vary_ellip) is tuple:
                options['artpop_model_params']['ellip'] = random.uniform(vary_ellip[0],vary_ellip[1])
            else: options['artpop_model_params']['ellip'] = vary_ellip
            ell = options['artpop_model_params']['ellip']
            print(f'ArtPop Ellipticity: {round(ell,4)}')
            ells_true.append(ell)
            name_stub += f'_ellip{round(ell,4)}'

            if type(vary_mass) is tuple:
                options['artpop_model_params']['total_mass'] = random.uniform(vary_mass[0],vary_mass[1])
            else: options['artpop_model_params']['total_mass'] = vary_mass
            mass = options['artpop_model_params']['total_mass']
            print(f'ArtPop Mass: {round(mass,0)}')
            masses.append(mass)
            name_stub += f'_mass{round(mass,0)}'

            if type(vary_dist) is tuple:
                options['artpop_model_params']['distance'] = random.uniform(vary_dist[0],vary_dist[1])
            else: options['artpop_model_params']['distance'] = vary_dist
            dist = options['artpop_model_params']['distance']
            print(f'ArtPop Distance: {round(dist,4)} Mpc')
            artpop_dists.append(dist)
            name_stub += f'_dist{round(dist,4)}'

            options['color1']['artpop_model_fn'] = name_stub+f'_iter{iter}'+'_r.fits'
            options['color2']['artpop_model_fn'] = name_stub+f'_iter{iter}'+'_b.fits'

            if args.use_only_one_artpop and iter == 0:
                model1, model2, src = artpop_functions.run_artimager(options)
                options['color1']['artpop_mag'] = src.sp.total_mag('Bessell_R')
                options['color2']['artpop_mag'] = src.sp.total_mag('Bessell_B')
                options['color1']['artpop_sbfmag'] = src.sp.sbf_mag('Bessell_R')

        elif args.random_sersics : # Random sersic profiles

            # Create random sersic model
            name_stub = 'sersic_model'
            model_type = 'sersic'
            if type(vary_radius) is tuple:
                radius = random.uniform(vary_radius[0],vary_radius[1])
            else: radius = vary_radius
            print(f'Sersic Radius: {round(radius,2)} pixels')
            sersic_radii.append(radius)
            name_stub += f'_radius{round(radius,4)}'

            if type(vary_n) is tuple:
                n = random.uniform(vary_n[0],vary_n[1])
            else: n = vary_n
            print(f'Sersic Index: {round(n,4)}')
            n_true.append(n)
            name_stub += f'_n{round(n,4)}'

            if type(vary_theta) is tuple:
                theta = random.uniform(vary_theta[0],vary_theta[1])
            else: theta = vary_theta
            print(f'Sersic PA: {round(theta,4)}')
            thetas.append(theta)
            name_stub += f'_theta{round(theta,4)}'

            if type(vary_ellip) is tuple:
                ell = random.uniform(vary_ellip[0],vary_ellip[1])
            else: ell = vary_ellip
            print(f'Sersic Ellipticity: {round(ell,4)}')
            ells_true.append(ell)
            name_stub += f'_ell{round(ell,4)}'

            if type(vary_Ie_b) is tuple:
                I_e_b = random.uniform(vary_Ie_b[0],vary_Ie_b[1])
            else: I_e_b = vary_Ie_b
            print(f'Sersic I_e_b: {round(I_e_b,4)}')
            I_e_b_true.append(I_e_b)
            name_stub += f'_Ieb{round(I_e_b,4)}'

            if type(vary_Ie_r) is tuple:
                I_e_r = random.uniform(vary_Ie_r[0],vary_Ie_r[1])
            else: I_e_r = vary_Ie_r
            print(f'Sersic I_e_r: {round(I_e_r,4)}')
            I_e_r_true.append(I_e_r)
            name_stub += f'_Ier{round(I_e_r,4)}'


            params_b = {'I_e': I_e_b_true[-1],'r_e': sersic_radii[-1],  'n': n_true[-1], 'X0': 525.0, 'Y0': 525.0, 'ell': ells_true[-1],'PA': thetas[-1]}
            params_r = {'I_e': I_e_r_true[-1],'r_e': sersic_radii[-1],'n': n_true[-1],'X0': 525.0,'Y0': 525.0,'ell': ells_true[-1],'PA': thetas[-1]}

            zpt_b = options['color2']['zpt'] - 2.5 * np.log10(options['gain']/options['exposure_time']) #+ options['color2']['extinction']
            zpt_r = options['color1']['zpt'] - 2.5 * np.log10(options['gain']/options['exposure_time']) #+ options['color1']['extinction']
            sersic_b = pymfit.sersic.Sersic(params_b,zpt=zpt_b,pixscale=options['pixscale'])
            sersic_r = pymfit.sersic.Sersic(params_r,zpt=zpt_r,pixscale=options['pixscale'])

            zptb_corrected = zpt_b #+ options['color2']['color_term']*(sersic_b.m_tot-sersic_r.m_tot)
            zptr_corrected = zpt_r #+ options['color1']['color_term']*(sersic_b.m_tot-sersic_r.m_tot)
            sersic_b = pymfit.sersic.Sersic(params_b,zpt=zptb_corrected,pixscale=options['pixscale'])
            sersic_r = pymfit.sersic.Sersic(params_r,zpt=zptr_corrected,pixscale=options['pixscale'])

            psf_r = fits.getdata(os.path.join(options['out_dir'],options['color1']['psf']))
            psf_b = fits.getdata(os.path.join(options['out_dir'],options['color2']['psf']))

            sersic_r_file = fits.PrimaryHDU(data=sersic_r.array(shape=(1051,1051), psf=psf_r),header=fits.Header(params_r))
            sersic_b_file = fits.PrimaryHDU(data=sersic_b.array(shape=(1051,1051), psf=psf_b),header=fits.Header(params_b))
            sersic_r_file.writeto(os.path.join(options['out_dir'],name_stub+f'_iter{iter}'+'_r.fits'))
            sersic_b_file.writeto(os.path.join(options['out_dir'],name_stub+f'_iter{iter}'+'_b.fits'))
            options['color1']['artpop_model_fn'] = name_stub+f'_iter{iter}'+'_r.fits'
            options['color2']['artpop_model_fn'] = name_stub+f'_iter{iter}'+'_b.fits'

        # If creating random cutout from image, create cutout
        if options['random_inject_position']:

            regions = read_ds9(regions_fn)
            im = fits.getdata(os.path.join(options['image_dir'],options['color1']['image_fn']))

            star_x = 1294 #815
            star_y = 2210 #2080
            dists_from_star = []
            idx_used = []

            # Get random position from pre-selected regions
            idx = random.randrange(0, len(regions))
            while idx in idx_used:
                idx = random.randrange(0, len(regions))
            idx_used.append(idx)

            rand_x = int(regions[idx].center.x)
            rand_y = int(regions[idx].center.y)
            negone = np.argwhere(im == im[rand_y,rand_x])
            two = np.argwhere(im == im[star_y,star_x])
            dists_from_star.append(np.linalg.norm(negone - two))

            # Create cutout
            cutout_fn_r = os.path.join(options['image_dir'],options['color1']['image_fn'].replace('.fits','_cutout.fits'))
            cutout_fn_b = os.path.join(options['image_dir'],options['color2']['image_fn'].replace('.fits','_cutout.fits'))
            cutout_r = make_cutout(os.path.join(options['image_dir'],options['color1']['image_fn']), (rand_x,rand_y), (1051,1051), cutout_fn=cutout_fn_r)
            cutout_b = make_cutout(os.path.join(options['image_dir'],options['color2']['image_fn']), (rand_x,rand_y), (1051,1051), cutout_fn=cutout_fn_b)
            star_mask = make_cutout(full_star_mask_fn,(rand_x,rand_y), (1051,1051), cutout_fn=star_mask_fn)

            if use_skymodel:
                skymodel_cutout_r = make_cutout(skymodel_r_fn, (rand_x,rand_y), (1051,1051), cutout_fn=skymodel_cutout_r_fn)
                skymodel_cutout_b = make_cutout(skymodel_b_fn, (rand_x,rand_y), (1051,1051), cutout_fn=skymodel_cutout_b_fn)
                skymodel_cutouts = [skymodel_cutout_r,skymodel_cutout_b]


        # Run imfit
        bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, functions, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true = lbcred.modeling_imfit(config_filename=args.imfit_config, options=options, iter=iter, fn_stub=fn_stub, backmodel=skymodel_cutouts)

        options_sbf['color'] = color
        options_sbf['model_fn'] = model1_fn
        options_sbf['resid_fn'] = resid1_fn
        options_sbf['model_summary_fn'] = bestfit1_fn

        radii_recovered.append(bestfit2[b_comp]['r_e'])
        pa_recovered.append(bestfit2[b_comp]['PA'])
        ell_recovered.append(bestfit2[b_comp]['ell'])
        n_recovered.append(bestfit2[b_comp]['n'])
        xpos_recovered.append(bestfit2[b_comp]['X0'])
        ypos_recovered.append(bestfit2[b_comp]['Y0'])
        I_e_r_recovered.append(bestfit1[r_comp]['I_e'])
        I_e_b_recovered.append(bestfit2[b_comp]['I_e'])
        mags_r_recovered.append(mag1)
        mags_b_recovered.append(mag2)
        if args.random_sersics:
            sbf_mags_true.append(None)
            mags_r_true.append(sersic_r.m_tot - options['color1']['color_term']*(sersic_b.m_tot-sersic_r.m_tot) - options['color1']['extinction'])
            mags_b_true.append(sersic_b.m_tot - options['color2']['color_term']*(sersic_b.m_tot-sersic_r.m_tot) - options['color2']['extinction'])

        else:
            sbf_mags_true.append(sbf_mag_true)
            mags_r_true.append(mag1_true)
            mags_b_true.append(mag2_true)

        # Run SBF
        if args.run_sbf:
            sbf_mag, dist_a, dist_b = lbcred.modeling_sbf(config_filename=args.sbf_config, options=options_sbf, run_iter=iter, imfit_functions=functions)

            sbf_mags_recovered.append(sbf_mag)
            sbf_dist_a_recovered.append(dist_a)
            sbf_dist_b_recovered.append(dist_b)

        options['color1']['image_fn'] = orig_fn_r
        options['color2']['image_fn'] = orig_fn_b

        iter+=1

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Save results
    mags_r_true=np.asarray(mags_r_true)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_mag_r'), "wb")
    np.save(file, mags_r_true)
    file.close

    mags_b_true=np.asarray(mags_b_true)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_mag_b'), "wb")
    np.save(file, mags_b_true)
    file.close

    ages=np.asarray(ages)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_age'), "wb")
    np.save(file, ages)
    file.close

    fehs=np.asarray(fehs)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_feh'), "wb")
    np.save(file, fehs)
    file.close

    if args.random_artpops:
        artpop_radii=np.asarray(artpop_radii)
        file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_radius'), "wb")
        np.save(file, artpop_radii)
        file.close
    elif args.random_sersics:
        sersic_radii=np.asarray(sersic_radii)
        file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_radius'), "wb")
        np.save(file, sersic_radii)
        file.close

    n_true=np.asarray(n_true)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_n'), "wb")
    np.save(file, n_true)
    file.close

    thetas=np.asarray(thetas)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_theta'), "wb")
    np.save(file, thetas)
    file.close

    ells_true=np.asarray(ells_true)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_ellipticity'), "wb")
    np.save(file, ells_true)
    file.close

    masses=np.asarray(masses)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_mass'), "wb")
    np.save(file, masses)
    file.close

    artpop_dists=np.asarray(artpop_dists)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{args.num_iters}_distance'), "wb")
    np.save(file, artpop_dists)
    file.close

    if args.run_sbf:

        sbf_mags=np.asarray(sbf_mags_recovered)
        sbf_dist_a=np.asarray(sbf_dist_a_recovered)
        sbf_dist_b=np.asarray(sbf_dist_b_recovered)
        sbf_mags_true=np.asarray(sbf_mags_true)

        file = open(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{args.num_iters}_weighted_avg_sbf_mags'), "wb")
        np.save(file, sbf_mags)
        file.close
        file = open(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{args.num_iters}_weighted_avg_sbf_dist_a'), "wb")
        np.save(file, sbf_dist_a)
        file.close
        file = open(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{args.num_iters}_weighted_avg_sbf_dist_b'), "wb")
        np.save(file, sbf_dist_b)
        file.close
        sbf_mags_true=np.asarray(sbf_mags_true)
        file = open(os.path.join(options['out_dir'],f'artpop_parameters_{fn_stub}_num-iters{args.num_iters}_true_sbf_magnitude'), "wb")
        np.save(file, sbf_mags_true)
        file.close

        plt.clf()
        plt.hist(sbf_mags-sbf_mags_true,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean(sbf_mags-sbf_mags_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(sbf_mags-sbf_mags_true),3)}')
        plt.axvline(x=np.median(sbf_mags-sbf_mags_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(sbf_mags-sbf_mags_true),3)}')
        plt.xlabel('SBF Magnitude Error')
        plt.legend(prop={'size': 15})
        plt.savefig(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{args.num_iters}_sbfmagnitude.png'))

        plt.clf()
        plt.hist((sbf_dist_a/1e6)-artpop_dists,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((sbf_dist_a/1e6)-artpop_dists), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_a/1e6)-artpop_dists),3)}')
        plt.axvline(x=np.median((sbf_dist_a/1e6)-artpop_dists), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_a/1e6)-artpop_dists),3)}')
        plt.xlabel('SBF Distance Error A (Jerjen Eqn. 2)')
        plt.legend(prop={'size': 15})
        plt.savefig(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{args.num_iters}_distance_a.png'))

        plt.clf()
        plt.hist((sbf_dist_b/1e6)-artpop_dists,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((sbf_dist_b/1e6)-artpop_dists), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_b/1e6)-artpop_dists),3)}')
        plt.axvline(x=np.median((sbf_dist_b/1e6)-artpop_dists), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_b/1e6)-artpop_dists),3)}')
        plt.xlabel('SBF Distance Error B (Jerjen Eqn. 3)')
        plt.legend(prop={'size': 15})
        plt.savefig(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{args.num_iters}_distance_b.png'))
