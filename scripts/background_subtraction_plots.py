import os
import matplotlib.pyplot as plt
import numpy as np
from lbcred.utils import misc

in_dir = '/Users/kirstencasey/real_image_tests/background_subtraction_tests/imfit_sbf'
save_dir = '/Users/kirstencasey/real_image_tests/background_subtraction_tests/imfit_sbf'
params = ['measured_mags_r','measured_mags_b','measured_sbfmags','measured_dists','measured_radii','measured_ellip','measured_n','measured_pa','measured_xpos','measured_ypos','measured_Ier',
            'measured_Ieb','true_mags_r','true_mags_b','true_sbfmags','true_dists','true_radii','true_ellip','true_n','true_pa','position_ids','artpop_ids','background_models']

back_models = ['polynomial'] #'SEsky','median',
colors_back = ['indigo'] #['gold','indigo','firebrick']
diff_param_vals_pos = [1] #[4,2,3,1]
colors_pos = ['tomato'] #['navy','darkgreen','gold','tomato']

# Get results from files
for model in back_models:
    temp_results = dict()
    file_dir = os.path.join(in_dir)
    for param in params:
        file = open(os.path.join(file_dir,f'background_subtraction_test_results_{param}'), "rb")
        temp_results[param] = np.load(file,allow_pickle=True)
        file.close
    if model == back_models[0]:
        results = temp_results
    else:
        for key in list(temp_results.keys()):
            results[key] = np.concatenate((results[key],temp_results[key]))

results['measured_color'] = results['measured_mags_b'] - results['measured_mags_r']
results['true_color'] = results['true_mags_b'] - results['true_mags_r']

# Make plots...

measured_params = ['measured_mags_b','measured_mags_r','measured_color','measured_sbfmags','measured_n','measured_ellip','measured_radii']
true_params = ['true_mags_b','true_mags_r','true_color','true_sbfmags','true_n','true_ellip','true_radii']
labels = ['b-band magnitude','r-band magnitude','color','SBF magnitude','Sersic index','Ellipticity','Fractional radius']
save_names = ['mag_b','mag_r','color','sbfmag','n','ellip','radius']

# as a function of background subtraction method
diff_param_name = 'background_models'
diff_param_vals_back = back_models
for measured, true, label, save in zip(measured_params,true_params,labels,save_names):

    for diff_param,color in zip(diff_param_vals_back,colors_back):

        measured_arr = misc.filter_dict(results,measured,{diff_param_name:diff_param})
        true_arr = misc.filter_dict(results,true,{diff_param_name:diff_param})

        if save == 'radius':
            plt.hist((true_arr-measured_arr)/true_arr,bins='auto',color=color,alpha=0.4, label=diff_param, histtype='stepfilled', ec="k", density=True)
            plt.axvline(x=np.median((true_arr-measured_arr)/true_arr), color=color, linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median((true_arr-measured_arr)/true_arr),3)}')

        else:
            plt.hist(true_arr-measured_arr,bins='auto',color=color,alpha=0.4, label=diff_param, histtype='stepfilled', ec="k", density=True)
            plt.axvline(x=np.median(true_arr-measured_arr), color=color, linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median(true_arr-measured_arr),3)}')

    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15})
    plt.savefig(os.path.join(save_dir,f'background_test_results_backgroundtype_{save}.png'))
    plt.clf()

# as a function of position
diff_param_name = 'position_ids'
for measured, true, label, save in zip(measured_params,true_params,labels,save_names):

    for diff_param,color in zip(diff_param_vals_pos,colors_pos):

        measured_arr = misc.filter_dict(results,measured,{diff_param_name:diff_param})
        true_arr = misc.filter_dict(results,true,{diff_param_name:diff_param})

        if save == 'radius':
            plt.hist((true_arr-measured_arr)/true_arr,bins='auto',color=color,alpha=0.4, label=diff_param, histtype='stepfilled', ec="k", density=True)
            plt.axvline(x=np.median((true_arr-measured_arr)/true_arr), color=color, linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median((true_arr-measured_arr)/true_arr),3)}')

        else:
            plt.hist(true_arr-measured_arr,bins='auto',color=color,alpha=0.4, label=diff_param, histtype='stepfilled', ec="k", density=True)
            plt.axvline(x=np.median(true_arr-measured_arr), color=color, linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median(true_arr-measured_arr),3)}')

    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15})
    plt.savefig(os.path.join(save_dir,f'background_test_results_position_{save}.png'))
    plt.clf()

xpos = misc.filter_dict(results,desired_param='measured_xpos',conditions={'position_ids':1})
ypos = misc.filter_dict(results,desired_param='measured_ypos',conditions={'position_ids':1})

plt.scatter(xpos,ypos)
plt.xlabel('x-pos')
plt.ylabel('y-pos')
plt.axvline(x=525, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=525, color='black', linestyle='dashed', linewidth=3)
plt.savefig(in_dir+f'imfit_results_position.png')
