import os
import matplotlib.pyplot as plt
import numpy as np
from lbcred.utils import misc
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

in_dirs = ['/Users/kirstencasey/real_image_tests/background_subtraction_tests/imfit_sbf_SEsky_TP_fix_all_struc','/Users/kirstencasey/real_image_tests/background_subtraction_tests/imfit_sbf_polynomial_TP_fix_all_struc']#,'/Users/kirstencasey/real_image_tests/background_subtraction_tests/imfit_sbf_median']
save_dir = '/Users/kirstencasey/real_image_tests/background_subtraction_tests/results_nostarsubtraction_TP_fix_all_struc'
params = ['measured_mags_r','measured_mags_b','measured_sbfmags','measured_dists','measured_radii','measured_ellip','measured_n','measured_pa','measured_xpos','measured_ypos','measured_Ier',
            'measured_Ieb','true_mags_r','true_mags_b','true_sbfmags','true_dists','true_radii','true_ellip','true_n','true_pa','position_ids','artpop_ids','background_models']

back_models = ['SEsky','polynomial']#,'median']
colors_back = ['indigo','gold']#['gold','indigo','firebrick']
diff_param_vals_pos = [1,2] #[4,2,3,1]
colors_pos = ['gold','tomato'] #['navy','darkgreen','gold','tomato']
shapes_pos = ['o','^']

norm = lambda x: x/len(back_models)
cmap = plt.cm.get_cmap('viridis')

# Get results from files
for model,in_dir in zip(back_models,in_dirs):
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
results['true_dists'] *= 1e6
results['true_radii_pixels'] = misc.parsecs_to_pixels(results['true_radii'], results['true_dists'], 0.224)
results['measured_ellip'][np.where(results['measured_ellip']<0)] = (results['measured_ellip'][np.where(results['measured_ellip']<0)])/(results['measured_ellip'][np.where(results['measured_ellip']<0)]-1)

# Make plots...

measured_params = ['measured_mags_b','measured_mags_r','measured_color','measured_sbfmags','measured_n','measured_ellip','measured_radii','measured_pa']
true_params = ['true_mags_b','true_mags_r','true_color','true_sbfmags','true_n','true_ellip','true_radii_pixels','true_pa']
labels = ['b-band magnitude','r-band magnitude','color','SBF magnitude','Sersic index','Ellipticity','Fractional radius','Position Angle']
save_names = ['mag_b','mag_r','color','sbfmag','n','ellip','radius', 'pa']

# as a function of background subtraction method
diff_param_name = 'background_models'
diff_param_vals_back = back_models
for measured, true, label, save in zip(measured_params,true_params,labels,save_names):

    for diff_param,color in zip(diff_param_vals_back,colors_back):

        measured_arr = misc.filter_dict(results,measured,{diff_param_name:diff_param})
        true_arr = misc.filter_dict(results,true,{diff_param_name:diff_param})

        if save == 'radius':
            plt.hist((true_arr-measured_arr)/true_arr,bins='auto',color=cmap(norm(back_models.index(diff_param))),alpha=0.4, label=diff_param, histtype='stepfilled', ec="k")
            plt.axvline(x=np.median((true_arr-measured_arr)/true_arr), color=cmap(norm(back_models.index(diff_param))), linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median((true_arr-measured_arr)/true_arr),3)}')

        elif save == 'pa':
            plt.hist(true_arr-measured_arr,bins='auto',color=cmap(norm(back_models.index(diff_param))),alpha=0.4, label=diff_param, histtype='stepfilled', ec="k")
            plt.axvline(x=np.median(true_arr-measured_arr), color=cmap(norm(back_models.index(diff_param))), linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median(true_arr-measured_arr),3)}')
            plt.axvline(x=180, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=-180, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=360, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=-360, color='black', linestyle='dashed', linewidth=3)

        else:
            plt.hist(true_arr-measured_arr,bins='auto',color=cmap(norm(back_models.index(diff_param))),alpha=0.4, label=diff_param, histtype='stepfilled', ec="k")
            plt.axvline(x=np.median(true_arr-measured_arr), color=cmap(norm(back_models.index(diff_param))), linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median(true_arr-measured_arr),3)}')

    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig(os.path.join(save_dir,f'background_test_results_backgroundtype_{save}.png'))
    plt.clf()

# as a function of position
diff_param_name = 'position_ids'
for measured, true, label, save in zip(measured_params,true_params,labels,save_names):

    for diff_param,color in zip(diff_param_vals_pos,colors_pos):

        measured_arr = misc.filter_dict(results,measured,{diff_param_name:diff_param})
        true_arr = misc.filter_dict(results,true,{diff_param_name:diff_param})

        if save == 'radius':
            plt.hist((true_arr-measured_arr)/true_arr,bins='auto',color=color,alpha=0.4, label=diff_param, histtype='stepfilled', ec="k")
            plt.axvline(x=np.median((true_arr-measured_arr)/true_arr), color=color, linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median((true_arr-measured_arr)/true_arr),3)}')

        elif save == 'pa':
            plt.hist(true_arr-measured_arr,bins='auto',color=color,alpha=0.4, label=diff_param, histtype='stepfilled', ec="k")
            plt.axvline(x=np.median(true_arr-measured_arr), color=color, linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median(true_arr-measured_arr),3)}')
            plt.axvline(x=180, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=-180, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=360, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=-360, color='black', linestyle='dashed', linewidth=3)

        else:
            plt.hist(true_arr-measured_arr,bins='auto',color=color,alpha=0.4, label=diff_param, histtype='stepfilled', ec="k")
            plt.axvline(x=np.median(true_arr-measured_arr), color=color, linestyle='dashed', linewidth=3, label=f'{diff_param} median: {round(np.median(true_arr-measured_arr),3)}')

    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig(os.path.join(save_dir,f'background_test_results_position_{save}.png'))
    plt.clf()

xpos = misc.filter_dict(results,desired_param='measured_xpos',conditions={'position_ids':1})
ypos = misc.filter_dict(results,desired_param='measured_ypos',conditions={'position_ids':1})


plt.scatter(xpos,ypos)
plt.xlabel('x-pos')
plt.ylabel('y-pos')
plt.axvline(x=525, color='black', linestyle='dashed', linewidth=3)
plt.axhline(y=525, color='black', linestyle='dashed', linewidth=3)
plt.savefig(os.path.join(save_dir,f'imfit_results_position.png'))
plt.clf()

################################################################################

# as a function of background subtraction method
diff_param_name = 'background_models'
diff_param_vals_back = back_models
used_param_combos = []
for param1,true_param1,label1,save1 in zip(measured_params,true_params,labels,save_names):
    for param2,true_param2,label2,save2 in zip(measured_params,true_params,labels,save_names):
        if param1 != param2:
            used_param_combos.append(param1+param2)
            if param2+param1 not in used_param_combos:
                for diff_param,color in zip(diff_param_vals_back,colors_back):

                    for diff_param_pos,shape in zip(diff_param_vals_pos,shapes_pos):

                        measured_arr1 = misc.filter_dict(results,param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                        true_arr1 = misc.filter_dict(results,true_param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                        measured_arr2 = misc.filter_dict(results,param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                        true_arr2 = misc.filter_dict(results,true_param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})

                        x_arr = true_arr1-measured_arr1
                        y_arr = true_arr2-measured_arr2

                        if save1 == 'radius':
                            x_arr = (true_arr1-measured_arr1)/true_arr1
                        elif save2 == 'radius':
                            y_arr = (true_arr2-measured_arr2)/true_arr2

                        plt.scatter(x_arr, y_arr,color=cmap(norm(back_models.index(diff_param))),alpha=0.4,marker=shape)

                if save1 == 'pa':
                    plt.axvline(x=180, color='black', linestyle='dashed', linewidth=3)
                    plt.axvline(x=-180, color='black', linestyle='dashed', linewidth=3)
                    plt.axvline(x=360, color='black', linestyle='dashed', linewidth=3)
                    plt.axvline(x=-360, color='black', linestyle='dashed', linewidth=3)
                elif save2 == 'pa':
                    plt.axhline(y=180, color='black', linestyle='dashed', linewidth=3)
                    plt.axhline(y=-180, color='black', linestyle='dashed', linewidth=3)
                    plt.axhline(y=360, color='black', linestyle='dashed', linewidth=3)
                    plt.axhline(y=-360, color='black', linestyle='dashed', linewidth=3)

                plt.xlabel(f'{label1} Error')
                plt.ylabel(f'{label2} Error')
                plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
                plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
                legend_elements = [Line2D([0], [0], marker='o', color='white', label='position 1',markerfacecolor='black', markersize=10),
                               Line2D([0], [0], marker='^', color='white',label='position 2',markerfacecolor='black', markersize=10),
                               Patch(facecolor=cmap(norm(0)),label=back_models[0]),
                               Patch(facecolor=cmap(norm(1)),label=back_models[1])]
                plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left',handles=legend_elements)
                plt.savefig(os.path.join(save_dir,f'background_test_results_backgroundtype_{save1}_{save2}.png'))
                plt.clf()

'''
# Look at radius outliers:

results['fractional_radius_error'] = (results['true_radii_pixels']-results['measured_radii'])/results['true_radii_pixels']
radius_outliers = misc.filter_dict(results, conditions={'fractional_radius_error' : 0.1}, condition_type={'fractional_radius_error': 'greater than'})
radius_outliers_small = misc.filter_dict(radius_outliers, conditions={'fractional_radius_error' : 0.3}, condition_type={'fractional_radius_error': 'less than'})
radius_outliers_large = misc.filter_dict(radius_outliers, conditions={'fractional_radius_error' : 0.3}, condition_type={'fractional_radius_error': 'greater than'})


# as a function of background subtraction method
diff_param_name = 'background_models'
diff_param_vals_back = back_models
for param1,true_param1,label1,save1 in zip(measured_params,true_params,labels,save_names):
    for param2,true_param2,label2,save2 in zip(measured_params,true_params,labels,save_names):

        for diff_param,color in zip(diff_param_vals_back,colors_back):

            for diff_param_pos,shape in zip(diff_param_vals_pos,shapes_pos):

                measured_arr1 = misc.filter_dict(radius_outliers_large,param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr1 = misc.filter_dict(radius_outliers_large,true_param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                measured_arr2 = misc.filter_dict(radius_outliers_large,param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr2 = misc.filter_dict(radius_outliers_large,true_param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})

                x_arr = true_arr1
                y_arr = true_arr2-measured_arr2

                if save2 == 'radius':
                    y_arr = (true_arr2-measured_arr2)/true_arr2

                plt.scatter(x_arr, y_arr,color=cmap(norm(back_models.index(diff_param))),alpha=0.4,marker=shape)

        plt.xlabel(f'{label1}')
        plt.ylabel(f'{label2} Error')
        legend_elements = [Line2D([0], [0], marker='o', color='white', label='position 1',markerfacecolor='black', markersize=10),
                       Line2D([0], [0], marker='^', color='white',label='position 2',markerfacecolor='black', markersize=10),
                       Patch(facecolor=cmap(norm(0)),label=back_models[0]),
                       Patch(facecolor=cmap(norm(1)),label=back_models[1])]
        plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left',handles=legend_elements)
        plt.savefig(os.path.join(save_dir,f'background_test_results_backgroundtype_large_radius_outliers_{save1}_{save2}error.png'))
        plt.clf()


#########

# as a function of background subtraction method
diff_param_name = 'background_models'
diff_param_vals_back = back_models
for param1,true_param1,label1,save1 in zip(measured_params,true_params,labels,save_names):
    for param2,true_param2,label2,save2 in zip(measured_params,true_params,labels,save_names):

        for diff_param,color in zip(diff_param_vals_back,colors_back):

            for diff_param_pos,shape in zip(diff_param_vals_pos,shapes_pos):

                measured_arr1 = misc.filter_dict(radius_outliers_small,param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr1 = misc.filter_dict(radius_outliers_small,true_param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                measured_arr2 = misc.filter_dict(radius_outliers_small,param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr2 = misc.filter_dict(radius_outliers_small,true_param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})

                x_arr = true_arr1
                y_arr = true_arr2-measured_arr2

                if save2 == 'radius':
                    y_arr = (true_arr2-measured_arr2)/true_arr2

                plt.scatter(x_arr, y_arr,color=cmap(norm(back_models.index(diff_param))),alpha=0.4,marker=shape)

        plt.xlabel(f'{label1}')
        plt.ylabel(f'{label2} Error')
        legend_elements = [Line2D([0], [0], marker='o', color='white', label='position 1',markerfacecolor='black', markersize=10),
                       Line2D([0], [0], marker='^', color='white',label='position 2',markerfacecolor='black', markersize=10),
                       Patch(facecolor=cmap(norm(0)),label=back_models[0]),
                       Patch(facecolor=cmap(norm(1)),label=back_models[1])]
        plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left',handles=legend_elements)
        plt.savefig(os.path.join(save_dir,f'background_test_results_backgroundtype_small_radius_outliers_{save1}_{save2}error.png'))
        plt.clf()


##############

good_radii = misc.filter_dict(results, conditions={'fractional_radius_error' : 0.1}, condition_type={'fractional_radius_error': 'less than'})

# as a function of background subtraction method
diff_param_name = 'background_models'
diff_param_vals_back = back_models
for param1,true_param1,label1,save1 in zip(measured_params,true_params,labels,save_names):
    for param2,true_param2,label2,save2 in zip(measured_params,true_params,labels,save_names):

        for diff_param,color in zip(diff_param_vals_back,colors_back):

            for diff_param_pos,shape in zip(diff_param_vals_pos,shapes_pos):

                measured_arr1 = misc.filter_dict(good_radii,param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr1 = misc.filter_dict(good_radii,true_param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                measured_arr2 = misc.filter_dict(good_radii,param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr2 = misc.filter_dict(good_radii,true_param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})

                x_arr = true_arr1
                y_arr = true_arr2-measured_arr2

                if save2 == 'radius':
                    y_arr = (true_arr2-measured_arr2)/true_arr2

                plt.scatter(x_arr, y_arr,color=cmap(norm(back_models.index(diff_param))),alpha=0.4,marker=shape)

        plt.xlabel(f'{label1}')
        plt.ylabel(f'{label2} Error')
        legend_elements = [Line2D([0], [0], marker='o', color='white', label='position 1',markerfacecolor='black', markersize=10),
                       Line2D([0], [0], marker='^', color='white',label='position 2',markerfacecolor='black', markersize=10),
                       Patch(facecolor=cmap(norm(0)),label=back_models[0]),
                       Patch(facecolor=cmap(norm(1)),label=back_models[1])]
        plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left',handles=legend_elements)
        plt.savefig(os.path.join(save_dir,f'background_test_results_backgroundtype_good_radii_{save1}_{save2}error.png'))
        plt.clf()

'''
'''
# as a function of background subtraction method
diff_param_name = 'background_models'
diff_param_vals_back = back_models
used_param_combos = []
for param1,true_param1,label1,save1 in zip(measured_params,true_params,labels,save_names):
    for param2,true_param2,label2,save2 in zip(measured_params,true_params,labels,save_names):


        for diff_param,color in zip(diff_param_vals_back,colors_back):

            for diff_param_pos,shape in zip(diff_param_vals_pos,shapes_pos):

                measured_arr1 = misc.filter_dict(results,param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr1 = misc.filter_dict(results,true_param1,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                measured_arr2 = misc.filter_dict(results,param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})
                true_arr2 = misc.filter_dict(results,true_param2,{diff_param_name:diff_param,'position_ids':diff_param_pos})

                x_arr = true_arr1
                y_arr = true_arr2-measured_arr2

                if save2 == 'radius':
                    y_arr = (true_arr2-measured_arr2)/true_arr2

                plt.scatter(x_arr, y_arr,color=color,label=diff_param,alpha=0.4,marker=shape)
                #plt.scatter(x_arr, y_arr,color=norm(1),alpha=0.4,marker=shape)

        if save1 == 'pa':
            plt.axvline(x=180, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=-180, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=360, color='black', linestyle='dashed', linewidth=3)
            plt.axvline(x=-360, color='black', linestyle='dashed', linewidth=3)
        elif save2 == 'pa':
            plt.axhline(y=180, color='black', linestyle='dashed', linewidth=3)
            plt.axhline(y=-180, color='black', linestyle='dashed', linewidth=3)
            plt.axhline(y=360, color='black', linestyle='dashed', linewidth=3)
            plt.axhline(y=-360, color='black', linestyle='dashed', linewidth=3)

        plt.xlabel(f'{label1} Error')
        plt.ylabel(f'{label2} Error')
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
        #legend_elements = [Line2D([0], [0], marker='o', color='black', label='position 1',markerfacecolor='black', markersize=10),
        #               Line2D([0], [0], marker='^', color='black', label='position 2',markerfacecolor='black', markersize=10),
        #               Patch(facecolor=cmap(norm(1)),label=back_models[0])]
        plt.legend(prop={'size': 15})
        plt.savefig(os.path.join(save_dir,f'background_test_results_backgroundtype1_{save1}_{save2}.png'))
        plt.clf()

# as a function of position
diff_param_name = 'position_ids'
used_param_combos = []
for param1,true_param1,label1,save1 in zip(measured_params,true_params,labels,save_names):
    for param2,true_param2,label2,save2 in zip(measured_params,true_params,labels,save_names):
        if param1 != param2:
            used_param_combos.append(param1+param2)
            if param2+param1 not in used_param_combos:
                for diff_param,color in zip(diff_param_vals_pos,colors_pos):

                    measured_arr1 = misc.filter_dict(results,param1,{diff_param_name:diff_param})
                    true_arr1 = misc.filter_dict(results,true_param1,{diff_param_name:diff_param})
                    measured_arr2 = misc.filter_dict(results,param2,{diff_param_name:diff_param})
                    true_arr2 = misc.filter_dict(results,true_param2,{diff_param_name:diff_param})

                    x_arr = true_arr1-measured_arr1
                    y_arr = true_arr2-measured_arr2

                    if save1 == 'radius':
                        x_arr = (true_arr1-measured_arr1)/true_arr1
                    elif save2 == 'radius':
                        y_arr = (true_arr2-measured_arr2)/true_arr2

                    plt.scatter(x_arr, y_arr,color=color,label=diff_param,alpha=0.4)

                if save1 == 'pa':
                    plt.axvline(x=180, color='black', linestyle='dashed', linewidth=3)
                    plt.axvline(x=-180, color='black', linestyle='dashed', linewidth=3)
                    plt.axvline(x=360, color='black', linestyle='dashed', linewidth=3)
                    plt.axvline(x=-360, color='black', linestyle='dashed', linewidth=3)
                elif save2 == 'pa':
                    plt.axhline(y=180, color='black', linestyle='dashed', linewidth=3)
                    plt.axhline(y=-180, color='black', linestyle='dashed', linewidth=3)
                    plt.axhline(y=360, color='black', linestyle='dashed', linewidth=3)
                    plt.axhline(y=-360, color='black', linestyle='dashed', linewidth=3)

                plt.xlabel(f'{label1} Error')
                plt.ylabel(f'{label2} Error')
                plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
                plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
                plt.legend(prop={'size': 15})
                plt.savefig(os.path.join(save_dir,f'background_test_results_position_{save1}_{save2}.png'))
                plt.clf()
'''
