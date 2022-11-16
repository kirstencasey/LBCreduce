import lbcred, os, pandas
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from lbcred import sbf_functions
from lbcred.utils import misc, io
from astropy.io import fits

num_iters = 50
sbf_mags_measured = []
sbf_mags_true = []

r_mags_measured = []
r_mags_true = []

b_mags_measured = []
b_mags_true = []
save_fn = os.path.join('/fs/scratch/PCON0003/osu10713/fixed_artpop_blank_field_tests/sbf_result_with_blank_fields.png')

def create_sbf_plot(config, galaxy_results, residual_image, save_fn, k_min, k_max, cutout_size):
    blank_field_results = []
    blank_resids = []
    blank_masks = []
    sbf_mags = []
    
    blank_field_positions = Table.from_pandas(pandas.read_csv(config['pre_selected_fields']))
    psf = fits.open(os.path.join(config['image_dir'],config['psf']))[config['ext']].data
    psf = psf/np.sum(psf)
    
    if config['model_summary_fn'] != None:
        bestfit_params = io.read_results(os.path.join(config['image_dir'], config['model_summary_fn']), ['TiltedSkyPlane','Sersic'])
        comp = ['TiltedSkyPlane','Sersic'].index('Sersic')+1
        sersic_results = bestfit_params[f'comp_{comp}']

    for blank_num in range(len(blank_field_positions)):
        pos = blank_field_positions[blank_num]
        misc.make_cutout(os.path.join(config['image_dir'],config['stack_fn']), (pos['xpos'],pos['ypos']), (cutout_size,cutout_size), ext=config['ext'], cutout_fn=os.path.join(config['image_dir'],f'blank_field_{blank_num}_r.fits'))
        
        # Get resid, mask
        blank_resid, blank_mask = sbf_functions.get_sbf_mask_resid(os.path.join(config['image_dir'],config['model_fn']), os.path.join(config['image_dir'],f'blank_field_{blank_num}_r.fits'), sersic_results, config['masking_sbf']['grow_obj'], config['masking_sbf']['fixed_frac_of_radius'], config, blank_field=True)
        blank_resids.append(blank_resid)
        blank_masks.append(blank_mask)
        
        blank_results = sbf_functions.measure_sbf(blank_resid[config['ext']].data, psf, mask=blank_mask, k_range=[k_min, k_max],
                        fit_param_guess=[100, 50], num_radial_bins=config['num_radial_bins'],
                        use_sigma=config['use_sigma'])
                        
        blank_field_results.append(blank_results.ps_image)
                        
        sbf_mag, d_a, d_b = sbf_functions.get_sbf_distance(results, config['zpt'], config['color'], config['gain'], config['exposure_time'], colorterm = config['color_term'], extinction_correction=config['extinction'], blank_field_results=blank_results)
        sbf_mags.append(sbf_mag)
        
    # Get median blank field result
    median_blank_field = blank_results
    median_blank_field_ps_image = np.median(blank_field_results,axis=0)
    blank_results_k=blank_results.k
    blank_results_npix=blank_results.npix
    
    sbf_functions.sbf_results(galaxy_results, residual_image[config['ext']].data, subplots=None, xlabel=r'Spacial Frequency (pixel$^{-1}$)',
                    ylabel=f'Power (r-band)', xscale='linear', percentiles=[config['plot_percentiles_min'], config['plot_percentiles_max']],
                    yscale='log', plot_errors=False, ylim_factors=[0.5, 1.1],
                    cmap='gray_r',save_fn=save_fn, plot_blank_fields=True,blank_field_results_k=blank_results_k,blank_field_results_ps_image=median_blank_field_ps_image,blank_field_results_npix=blank_results_npix)
    
    return

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
    parser.add_argument('--run_imfit', type=bool, default=True, help='do imfit modelling')
    parser.add_argument('--run_artpop', type=bool, help='do artpop modelling before running imfit')
    parser.add_argument('--run_sbf', type=bool, default=True, help='do sbf measurement')
    args = parser.parse_args()

    for run_num in range(num_iters):
        options = {
            'image_dir' : args.image_dir,
            'out_dir' : args.out_dir,
            'run_imfit' : args.run_imfit,
            'run_artpop' : args.run_artpop
        }
    
        if args.run_imfit:
            bestfit1, bestfit2, bestfit1_fn, bestfit2_fn, mag1, mag2, color, model1_fn, resid1_fn, functions, model2_fn, resid2_fn, sbf_mag_true, mag1_true, mag2_true, orig_color = lbcred.modeling_imfit(config_filename=args.imfit_config, options=options, iter=run_num)
            sbf_mags_true.append(sbf_mag_true)
            r_mags_measured.append(mag1)
            r_mags_true.append(mag1_true)
            b_mags_measured.append(mag2)
            b_mags_true.append(mag2_true)
            
        options.pop('run_imfit')
        options.pop('run_artpop')
        #options['color'] = color
        options['color'] = orig_color
        options['model_fn'] = model1_fn
        options['resid_fn'] = resid1_fn
        options['model_summary_fn'] = bestfit1_fn
    
        if args.run_sbf:
            sbf_mag, wavg_dist_a, wavg_dist_b, uncert_mag, results, sbf_resid, k_min, k_max, config, cutout_size = lbcred.modeling_sbf(config_filename=args.sbf_config, options=options, imfit_functions=functions, run_iter=run_num, make_paper_plot=True)
            sbf_mags_measured.append(sbf_mag)
            
            #create_sbf_plot(config, results, sbf_resid, save_fn, k_min, k_max, cutout_size)
    
    sbf_mags_measured = np.asarray(sbf_mags_measured)
    sbf_mags_true = np.asarray(sbf_mags_true)
    plt.figure(figsize=(8,8))
    plt.hist(sbf_mags_true-sbf_mags_measured,bins='auto',color='navy',alpha=0.4, histtype='stepfilled', ec="k")
    plt.axvline(x=np.median(sbf_mags_true-sbf_mags_measured), color='navy', linestyle='dashed', linewidth=3, label=f'median error: {round(np.median(sbf_mags_true-sbf_mags_measured),3)}')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'SBF magnitude error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig('/fs/scratch/PCON0003/osu10713/fixed_artpop_blank_field_tests/background_test_results_position_sbfmag_4mpc_samemaskdistance-5_justmessing1.png')
    plt.clf()
    
    r_mags_measured = np.asarray(r_mags_measured)
    r_mags_true = np.asarray(r_mags_true)
    plt.figure(figsize=(8,8))
    plt.hist(r_mags_true-r_mags_measured,bins='auto',color='navy',alpha=0.4, histtype='stepfilled', ec="k")
    plt.axvline(x=np.median(r_mags_true-r_mags_measured), color='navy', linestyle='dashed', linewidth=3, label=f'median error: {round(np.median(r_mags_true-r_mags_measured),3)}')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'r-band magnitude error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig('/fs/scratch/PCON0003/osu10713/fixed_artpop_blank_field_tests/background_test_results_position_rmag_4mpc_samemaskdistance-5_justmessing1.png')
    plt.clf()
    
    b_mags_measured = np.asarray(b_mags_measured)
    b_mags_true = np.asarray(b_mags_true)
    plt.figure(figsize=(8,8))
    plt.hist(b_mags_true-b_mags_measured,bins='auto',color='navy',alpha=0.4, histtype='stepfilled', ec="k")
    plt.axvline(x=np.median(b_mags_true-b_mags_measured), color='navy', linestyle='dashed', linewidth=3, label=f'median error: {round(np.median(b_mags_true-b_mags_measured),3)}')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'b-band magnitude error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig('/fs/scratch/PCON0003/osu10713/fixed_artpop_blank_field_tests/background_test_results_position_bmag_4mpc_samemaskdistance-5_justmessing1.png')
    plt.clf()
    
    
    
    
    
    
    
