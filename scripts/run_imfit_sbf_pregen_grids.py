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
    parser.add_argument('--run_sbf', type=bool, default=True, help='do sbf measurement')
    args = parser.parse_args()

    fn_stub = 'artpop_grid_tests_including-noise_fixed_on_b_restricted_ell'
    use_skymodel = False

    # Open and read imfit config file
    with open(args.imfit_config, 'r') as filename:
        options = yaml.load(filename, Loader=yaml.FullLoader)

    options['inject_artpop_model'] = True
    model_type = 'artpop'

    # Get imfit r- and b-band components given imfit steps in config file
    step_num = 1
    for step in options['imfit_steps']:
        if options['imfit_steps'][step]['color'] == 'b':
            b_step = step_num
        elif options['imfit_steps'][step]['color'] == 'r':
            r_step = step_num
        step_num+=1
    b_comp = misc.list_of_strings(options['imfit_steps'][f'step{b_step}']['functions']).index('Sersic')+1
    r_comp = misc.list_of_strings(options['imfit_steps'][f'step{r_step}']['functions']).index('Sersic')+1
    b_comp = f'comp_{b_comp}'
    r_comp = f'comp_{r_comp}'

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

    # Get pregen artpop models
    models_r = glob.glob(os.path.join(options['pregen_artpop_dir'],'artpop_model_*_r.fits'))
    models_b = glob.glob(os.path.join(options['pregen_artpop_dir'],'artpop_model_*_b.fits'))
    if len(models_r) is 0 : print('No ArtPop models found! Check pregen_artpop_dir in imfit configuration file!')
    models_r = np.sort(np.asarray(models_r))
    models_b = np.sort(np.asarray(models_b))
    options['run_artpop'] = False

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    iter = 0
    skymodel_cutouts = None

    # Run imfit on galaxy at various positions
    while iter < len(models_r):

        options['color1']['artpop_model_fn'] = models_r[iter]
        options['color2']['artpop_model_fn'] = models_b[iter]
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
        '''
        if args.random_sersics:
            sbf_mags_true.append(None)
            mags_r_true.append(sersic_r.m_tot - options['color1']['color_term']*(sersic_b.m_tot-sersic_r.m_tot) - options['color1']['extinction'])
            mags_b_true.append(sersic_b.m_tot - options['color2']['color_term']*(sersic_b.m_tot-sersic_r.m_tot) - options['color2']['extinction'])
        '''
        #else:
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
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_mag_r'), "wb")
    np.save(file, mags_r_true)
    file.close

    mags_b_true=np.asarray(mags_b_true)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_mag_b'), "wb")
    np.save(file, mags_b_true)
    file.close

    ages=np.asarray(ages)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_age'), "wb")
    np.save(file, ages)
    file.close

    fehs=np.asarray(fehs)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_feh'), "wb")
    np.save(file, fehs)
    file.close

    #if args.random_artpops:
    artpop_radii=np.asarray(artpop_radii)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_radius'), "wb")
    np.save(file, artpop_radii)
    file.close
    '''
    elif args.random_sersics:
        sersic_radii=np.asarray(sersic_radii)
        file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_radius'), "wb")
        np.save(file, sersic_radii)
        file.close
    '''

    n_true=np.asarray(n_true)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_n'), "wb")
    np.save(file, n_true)
    file.close

    thetas=np.asarray(thetas)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_theta'), "wb")
    np.save(file, thetas)
    file.close

    ells_true=np.asarray(ells_true)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_ellipticity'), "wb")
    np.save(file, ells_true)
    file.close

    masses=np.asarray(masses)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_mass'), "wb")
    np.save(file, masses)
    file.close

    artpop_dists=np.asarray(artpop_dists)
    file = open(os.path.join(options['out_dir'],f'{model_type}_parameters_{fn_stub}_num-iters{len(models_r)}_distance'), "wb")
    np.save(file, artpop_dists)
    file.close

    if args.run_sbf:

        sbf_mags=np.asarray(sbf_mags_recovered)
        sbf_dist_a=np.asarray(sbf_dist_a_recovered)
        sbf_dist_b=np.asarray(sbf_dist_b_recovered)
        sbf_mags_true=np.asarray(sbf_mags_true)

        file = open(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{len(models_r)}_weighted_avg_sbf_mags'), "wb")
        np.save(file, sbf_mags)
        file.close
        file = open(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{len(models_r)}_weighted_avg_sbf_dist_a'), "wb")
        np.save(file, sbf_dist_a)
        file.close
        file = open(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{len(models_r)}_weighted_avg_sbf_dist_b'), "wb")
        np.save(file, sbf_dist_b)
        file.close
        sbf_mags_true=np.asarray(sbf_mags_true)
        file = open(os.path.join(options['out_dir'],f'artpop_parameters_{fn_stub}_num-iters{len(models_r)}_true_sbf_magnitude'), "wb")
        np.save(file, sbf_mags_true)
        file.close

        plt.clf()
        plt.hist(sbf_mags-sbf_mags_true,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean(sbf_mags-sbf_mags_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(sbf_mags-sbf_mags_true),3)}')
        plt.axvline(x=np.median(sbf_mags-sbf_mags_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(sbf_mags-sbf_mags_true),3)}')
        plt.xlabel('SBF Magnitude Error')
        plt.legend(prop={'size': 15})
        plt.savefig(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{len(models_r)}_sbfmagnitude.png'))

        plt.clf()
        plt.hist((sbf_dist_a/1e6)-artpop_dists,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((sbf_dist_a/1e6)-artpop_dists), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_a/1e6)-artpop_dists),3)}')
        plt.axvline(x=np.median((sbf_dist_a/1e6)-artpop_dists), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_a/1e6)-artpop_dists),3)}')
        plt.xlabel('SBF Distance Error A (Jerjen Eqn. 2)')
        plt.legend(prop={'size': 15})
        plt.savefig(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{len(models_r)}_distance_a.png'))

        plt.clf()
        plt.hist((sbf_dist_b/1e6)-artpop_dists,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((sbf_dist_b/1e6)-artpop_dists), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_b/1e6)-artpop_dists),3)}')
        plt.axvline(x=np.median((sbf_dist_b/1e6)-artpop_dists), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_b/1e6)-artpop_dists),3)}')
        plt.xlabel('SBF Distance Error B (Jerjen Eqn. 3)')
        plt.legend(prop={'size': 15})
        plt.savefig(os.path.join(options['out_dir'],f'sbf_results_{fn_stub}_num-iters{len(models_r)}_distance_b.png'))
