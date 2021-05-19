import lbcred, os, yaml
from astropy.io import fits
from astropy.nddata import Cutout2D
import numpy as np
import matplotlib.pyplot as plt
import random
from lbcred.model import artpop_functions


def make_cutout(original_img_fn, position, shape, cutout_fn=None):
    img = fits.open(original_img_fn)
    cutout = Cutout2D(img[0].data, position, shape)
    img[0].data = cutout.data
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
    parser.add_argument('-i','--image_dir', type=str, help='name of directory containing all raw images to be used in analysis; if not running this program in the directory containing image_dir, make sure to include full path name')
    parser.add_argument('--out_dir', type=str, help='specify output directory name/path; default is located in the same directory as image_dir with name \'lbcreduce_modeling_<date>_<time>\'')
    parser.add_argument('--imfit_config', type=str, default=default_imfit_config, help='path of the .yml config file ')
    parser.add_argument('--sbf_config', type=str, default=default_sbf_config, help='path of the .yml config file ')
    parser.add_argument('--num-iters', type=int, default=1, help='number of positions to try placing artpop models')
    parser.add_argument('--random_inject_position', type=bool, default=False, help='inject artpop models at randomly chosen positions')
    parser.add_argument('--random_artpops', type=bool, default=False, help='do sbf measurement')
    parser.add_argument('--run_artpop_only_once', type=bool, default=True, help='do sbf measurement')
    parser.add_argument('--run_sbf', type=bool, default=True, help='do sbf measurement')

    fn_stub = 'blank_test' # 'random_position_single_artpop_sersic_fixed_on_b_fixed_all_struc_params_no_tp'
    b_comp = 'comp_1'
    r_comp = 'comp_1'
    exp_type = 'sci'
    vary_age = 10.1
    vary_feh = -1.05
    vary_radius = 894 # 600
    vary_n = 0.52078
    vary_theta = 0.
    vary_ellip = 0.
    vary_mass = 14900000.0 #8500000
    vary_dist = 3.7 # Mpc

    args = parser.parse_args()

    # Open and read config file
    with open(args.imfit_config, 'r') as filename:
        options = yaml.load(filename, Loader=yaml.FullLoader)

    # Open and read config file
    with open(args.sbf_config, 'r') as filename:
        options_sbf = yaml.load(filename, Loader=yaml.FullLoader)

    options['random_inject_position'] = args.random_inject_position

    options['run_artpop'] = args.random_artpops
    options['inject_artpop_model'] = True

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
    artpop_n = []
    thetas = []
    artpop_ells = []
    masses = []
    artpop_dists = []


    iter = 0

    # Run imfit on galaxy at various positions
    while iter < args.num_iters:

        if args.random_artpops or iter == 0:

            # Create random artpop model
            name_stub = 'artpop_model'
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
            artpop_n.append(n)
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
            artpop_ells.append(ell)
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

        if args.run_artpop_only_once and iter == 0:
            model1, model2, src = artpop_functions.run_artimager(options)
            options['color1']['artpop_mag'] = src.sp.total_mag('Bessell_R')
            options['color2']['artpop_mag'] = src.sp.total_mag('Bessell_B')
            options['color1']['artpop_sbfmag'] = src.sp.sbf_mag('Bessell_R')

        # If creating random cutout from image, create cutout
        if args.random_inject_position:

            from regions import read_ds9

            im_fn_r = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/mean_stack_r.fits'
            im_fn_b = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/mean_stack_b.fits'
            cutout_fn_r = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/random_position_r.fits'
            cutout_fn_b = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/random_position_b.fits'
            full_star_mask_fn = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/full_star_mask.fits'
            star_mask_fn = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/star_mask.fits'
            im = fits.getdata(im_fn_r)

            regions = read_ds9(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/inject_positions_careful.reg')

            star_x = 1294
            star_y = 2210
            dists_from_star = []
            idx_used = []

            # Get random position
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
            cutout_r = make_cutout(im_fn_r, (rand_x,rand_y), (1051,1051), cutout_fn=cutout_fn_r)
            cutout_b = make_cutout(im_fn_b, (rand_x,rand_y), (1051,1051), cutout_fn=cutout_fn_b)
            star_mask = make_cutout(full_star_mask_fn,(rand_x,rand_y), (1051,1051), cutout_fn=star_mask_fn)

        # Run imfit
        bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, functions, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true = lbcred.modeling_imfit(config_filename=args.imfit_config, options=options, iter=iter, fn_stub=fn_stub)

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

    mags_r_true=np.asarray(mags_r_true)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_mag_r', "wb")
    np.save(file, mags_r_true)
    file.close

    mags_b_true=np.asarray(mags_b_true)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_mag_b', "wb")
    np.save(file, mags_b_true)
    file.close

    ages=np.asarray(ages)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_age', "wb")
    np.save(file, ages)
    file.close

    fehs=np.asarray(fehs)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_feh', "wb")
    np.save(file, fehs)
    file.close

    artpop_radii=np.asarray(artpop_radii)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_radius', "wb")
    np.save(file, artpop_radii)
    file.close

    artpop_n=np.asarray(artpop_n)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_n', "wb")
    np.save(file, artpop_n)
    file.close

    thetas=np.asarray(thetas)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_theta', "wb")
    np.save(file, thetas)
    file.close

    artpop_ells=np.asarray(artpop_ells)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_ellipticity', "wb")
    np.save(file, artpop_ells)
    file.close

    masses=np.asarray(masses)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_mass', "wb")
    np.save(file, masses)
    file.close

    artpop_dists=np.asarray(artpop_dists)
    file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_distance', "wb")
    np.save(file, artpop_dists)
    file.close

    if args.run_sbf:

        sbf_mags=np.asarray(sbf_mags_recovered)
        sbf_dist_a=np.asarray(sbf_dist_a_recovered)
        sbf_dist_b=np.asarray(sbf_dist_b_recovered)
        sbf_mags_true=np.asarray(sbf_mags_true)

        file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_iters}_weighted_avg_sbf_mags', "wb")
        np.save(file, sbf_mags)
        file.close
        file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_iters}_weighted_avg_sbf_dist_a', "wb")
        np.save(file, sbf_dist_a)
        file.close
        file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_iters}_weighted_avg_sbf_dist_b', "wb")
        np.save(file, sbf_dist_b)
        file.close
        sbf_mags_true=np.asarray(sbf_mags_true)
        file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/artpop_parameters_{fn_stub}_num-iters{args.num_iters}_true_sbf_magnitude', "wb")
        np.save(file, sbf_mags_true)
        file.close

        plt.clf()
        plt.hist(sbf_mags-sbf_mags_true,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean(sbf_mags-sbf_mags_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(sbf_mags-sbf_mags_true),3)}')
        plt.axvline(x=np.median(sbf_mags-sbf_mags_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(sbf_mags-sbf_mags_true),3)}')
        plt.xlabel('SBF Magnitude Error')
        plt.legend(prop={'size': 15})
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_iters}_sbfmagnitude.png')

        plt.clf()
        plt.hist((sbf_dist_a/1e6)-artpop_dists,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((sbf_dist_a/1e6)-artpop_dists), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_a/1e6)-artpop_dists),3)}')
        plt.axvline(x=np.median((sbf_dist_a/1e6)-artpop_dists), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_a/1e6)-artpop_dists),3)}')
        plt.xlabel('SBF Distance Error A (Jerjen Eqn. 2)')
        plt.legend(prop={'size': 15})
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_iters}_distance_a.png')

        plt.clf()
        plt.hist((sbf_dist_b/1e6)-artpop_dists,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((sbf_dist_b/1e6)-artpop_dists), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_b/1e6)-artpop_dists),3)}')
        plt.axvline(x=np.median((sbf_dist_b/1e6)-artpop_dists), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_b/1e6)-artpop_dists),3)}')
        plt.xlabel('SBF Distance Error B (Jerjen Eqn. 3)')
        plt.legend(prop={'size': 15})
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_iters}_distance_b.png')
