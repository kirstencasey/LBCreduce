import lbcred, os, yaml
from astropy.io import fits
from astropy.nddata import Cutout2D
import numpy as np


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
    parser.add_argument('--num-positions', type=int, default=50, help='number of positions to try placing artpop models')
    parser.add_argument('--random_inject_position', type=bool, default=True, help='inject artpop models at randomly chosen positions')
    parser.add_argument('--run_sbf', type=bool, default=True, help='do sbf measurement')

    fn_stub = 'artpop_sersic_fixed_on_r_free_pa_ell'
    b_comp = 'comp_2'
    r_comp = 'comp_2'
    exp_type = 'sci'
    true_dist = 3.7 # Mpc

    args = parser.parse_args()

    options = {
        'image_dir' : args.image_dir,
        'out_dir' : args.out_dir,
        'random_inject_position' : args.random_inject_position,
        'run_sbf' : args.run_sbf
    }

    # Open and read config file
    with open(args.imfit_config, 'r') as filename:
        config = yaml.load(filename, Loader=yaml.FullLoader)

	# Replace options input via command line into config
    for key in options:
        if options[key] != None: config[key] = options[key]

    if args.random_inject_position:
        import random
        import matplotlib.pyplot as plt
        from regions import read_ds9

        im_fn_r = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/mean_stack_r.fits'
        im_fn_b = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/mean_stack_b.fits'
        cutout_fn_r = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/random_position_r.fits'
        cutout_fn_b = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/random_position_b.fits'
        full_star_mask_fn = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/full_star_mask.fits'
        star_mask_fn = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/star_mask.fits'
        im = fits.getdata(im_fn_r)
        '''
        mask_fn = f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/random_inject_mask.fits'
        mask = fits.getdata(mask_fn).astype(bool)
        goodvalues = [False]
        good_positions = np.isin(mask, goodvalues)
        y,x= np.where(good_positions)
        '''
        regions = read_ds9(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/inject_positions.reg')

        star_x = 1294
        star_y = 2210
        dists_from_star = []

        radii = []
        pa = []
        ell = []
        n = []
        I_e_r = []
        mags_r = []
        xpos = []
        ypos = []
        I_e_b = []
        mags_b = []
        idx_used = []
        sbf_mags = []
        sbf_dist_a = []
        sbf_dist_b = []

        iter = 0

        # Run imfit on galaxy at various positions
        while iter < args.num_positions:

            # Get random position
            '''
            pos_idx = random.randrange(len(y))
            rand_x = x[pos_idx]
            rand_y = y[pos_idx]
            '''
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

            bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true = lbcred.modeling_imfit(config_filename=args.imfit_config, options=options, iter=iter, fn_stub=fn_stub)

            radii.append(bestfit2[b_comp]['r_e'])
            pa.append(bestfit2[b_comp]['PA'])
            ell.append(bestfit2[b_comp]['ell'])
            n.append(bestfit2[b_comp]['n'])
            xpos.append(bestfit2[b_comp]['X0'])
            ypos.append(bestfit2[b_comp]['Y0'])
            I_e_r.append(bestfit1[r_comp]['I_e'])
            I_e_b.append(bestfit2[b_comp]['I_e'])
            mags_r.append(mag1)
            mags_b.append(mag2)
            '''
            if iter==0:
                options.pop('run_imfit')
                options.pop('run_artpop')
            '''
            options['color'] = color
            options['model_fn'] = model1_fn
            options['resid_fn'] = resid1_fn
            options['model_summary_fn'] = bestfit1_fn

            if args.run_sbf:
                sbf_mag, dist_a, dist_b = lbcred.modeling_sbf(config_filename=args.sbf_config, options=options, run_iter=iter)

                sbf_mags.append(sbf_mag)
                sbf_dist_a.append(dist_a)
                sbf_dist_b.append(dist_b)

            iter+=1

        if args.run_sbf:

            sbf_mags=np.asarray(sbf_mags)
            sbf_dist_a=np.asarray(sbf_dist_a)
            sbf_dist_b=np.asarray(sbf_dist_b)

            file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_positions}_weighted_avg_sbf_mags', "wb")
            np.save(file, sbf_mags)
            file.close
            file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_positions}_weighted_avg_sbf_dist_a', "wb")
            np.save(file, sbf_dist_a)
            file.close
            file = open(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_positions}_weighted_avg_sbf_dist_b', "wb")
            np.save(file, sbf_dist_b)
            file.close

            plt.clf()
            plt.hist(sbf_mags-sbf_mag_true,bins='auto',color='purple',alpha=0.8)
            plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=np.mean(sbf_mags-sbf_mag_true), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean(sbf_mags-sbf_mag_true),3)}')
            plt.axvline(x=np.median(sbf_mags-sbf_mag_true), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median(sbf_mags-sbf_mag_true),3)}')
            plt.xlabel('SBF Magnitude Error')
            plt.legend(prop={'size': 15})
            plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_positions}_sbfmagnitude.png')

            plt.clf()
            plt.hist((sbf_dist_a/1e6)-true_dist,bins='auto',color='purple',alpha=0.8)
            plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=np.mean((sbf_dist_a/1e6)-true_dist), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_a/1e6)-true_dist),3)}')
            plt.axvline(x=np.median((sbf_dist_a/1e6)-true_dist), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_a/1e6)-true_dist),3)}')
            plt.xlabel('SBF Distance Error A (Jerjen Eqn. 2)')
            plt.legend(prop={'size': 15})
            plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_positions}_distance_a.png')

            plt.clf()
            plt.hist((sbf_dist_b/1e6)-true_dist,bins='auto',color='purple',alpha=0.8)
            plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=np.mean((sbf_dist_b/1e6)-true_dist), color='firebrick', linestyle='dashed', linewidth=3, label=f'average: {round(np.mean((sbf_dist_b/1e6)-true_dist),3)}')
            plt.axvline(x=np.median((sbf_dist_b/1e6)-true_dist), color='gold', linestyle='dashed', linewidth=3, label=f'median: {round(np.median((sbf_dist_b/1e6)-true_dist),3)}')
            plt.xlabel('SBF Distance Error B (Jerjen Eqn. 3)')
            plt.legend(prop={'size': 15})
            plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/sbf_results_{fn_stub}_num-iters{args.num_positions}_distance_b.png')


        '''
        # Make histograms
        radii=np.asarray(radii)
        pa=np.asarray(pa)
        ell=np.asarray(ell)
        n=np.asarray(n)
        xpos=np.asarray(xpos)
        ypos=np.asarray(ypos)
        I_e_r=np.asarray(I_e_r)
        I_e_b=np.asarray(I_e_b)
        mags_r=np.asarray(mags_r)
        mags_b=np.asarray(mags_b)
        colors=mags_b-mags_r

        r_true = 224.133
        pa_true = 0
        ell_true = 0
        n_true = 0.568211
        I_e_r_true = 8.4 #89.7946
        mag_r_true = 14.531013184162031 #14.594793336161725
        xpos_true = 525.0
        ypos_true = 525.0
        I_e_b_true = 2.825 #23.825
        mag_b_true = 15.780570355571324 #15.861027188510267
        color_true = 1.2495571714092932 #1.2662338523485417

        plt.clf()
        num_bins = 40
        #bins = np.linspace((config['sersic_params']['r_e_min']-r_true)/r_true,(config['sersic_params']['r_e_max']-r_true)/r_true,num_bins)
        plt.hist((radii-r_true)/r_true,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((radii-r_true)/r_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Radius Error')
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_radii.png')
        plt.clf()
        #num_bins = 20
        ######################### Dont do this:
        bins = np.linspace((config['sersic_params']['PA_min']-pa_true)/pa_true,(config['sersic_params']['PA_max']-pa_true)/pa_true,num_bins)
        plt.hist((pa-pa_true)/pa_true,bins=bins,color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((pa-pa_true)/pa_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('PA Error')
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_pa.png')
        plt.clf()
        #num_bins = 20
        bins = np.linspace((config['sersic_params']['ell_min']-ell_true)/ell_true,(config['sersic_params']['ell_max']-ell_true)/ell_true,num_bins)
        plt.hist((ell-ell_true)/ell_true,bins=bins,color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((ell-ell_true)/ell_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Ellipticity Error')
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_ell.png')
        #######################
        plt.clf()
        #num_bins = 20
        #bins = np.linspace((config['sersic_params']['n_min']-n_true)/n_true,(config['sersic_params']['n_max']-n_true)/n_true,num_bins)
        plt.hist((n-n_true)/n_true,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((n-n_true)/n_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Sersic Index Error')
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_n.png')
        plt.clf()
        #num_bins = 20
        #bins = np.linspace((config['sersic_params']['I_e_min']-I_e_r_true)/I_e_r_true,(config['sersic_params']['I_e_max']-I_e_r_true)/I_e_r_true,num_bins)
        plt.hist((I_e_r-I_e_r_true)/I_e_r_true,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((I_e_r-I_e_r_true)/I_e_r_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Intensity Error (r-band)')
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_I_e_r.png')
        plt.clf()
        #num_bins = 20
        #bins = np.linspace((config['sersic_params']['I_e_min']-I_e_b_true)/I_e_b_true,(config['sersic_params']['I_e_max']-I_e_b_true)/I_e_b_true,num_bins)
        plt.hist((I_e_b-I_e_b_true)/I_e_b_true,bins='auto',color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((I_e_b-I_e_b_true)/I_e_b_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Intensity Error (b-band)')
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_I_e_b.png')
        plt.clf()
        #num_bins = 20
        #bins = np.linspace((14-mag_r_true)/mag_r_true,(17-mag_r_true)/mag_r_true,num_bins)
        plt.hist((mags_r-mag_r_true)/mag_r_true,bins='auto',range=(-0.5,0.5),color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((mags_r-mag_r_true)/mag_r_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Magnitude Error (r-band)')
        plt.savefig(f'/Users/kirstencasey/m8blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_mags_r.png')
        plt.clf()
        #num_bins = 20
        #bins = np.linspace((14-mag_b_true)/mag_b_true,(17-mag_b_true)/mag_b_true,num_bins)
        plt.hist((mags_b-mag_b_true)/mag_b_true,bins='auto',range=(-0.5,0.5),color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((mags_b-mag_b_true)/mag_b_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Magnitude Error (b-band)')
        plt.savefig(f'/Users/kirstencasey/m8blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_mags_b.png')
        plt.clf()
        #num_bins = 20
        #bins = np.linspace((0.5-color_true)/color_true,(2.0-color_true)/color_true,num_bins)
        plt.hist((colors-color_true)/color_true,bins='auto',range=(-0.7,0.7),color='purple',alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axvline(x=np.mean((colors-color_true)/color_true), color='firebrick', linestyle='dashed', linewidth=3, label='average')
        plt.xlabel('Color Error')
        plt.savefig(f'/Users/kirstencasey/m81blob_out/{exp_type}/lbc-reduce_testing/imfit_results_{fn_stub}_num-iters{args.num_positions}_colors.png')
        '''
    else: bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true = lbcred.modeling_imfit(config_filename=args.imfit_config, options=options)
