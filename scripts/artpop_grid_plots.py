import glob, os, yaml
import matplotlib.pyplot as plt
import numpy as np
from lbcred.utils import io, misc
from lbcred.model import imfit
from artpop.stars import MISTIsochrone
import matplotlib.colors as colors

# Get results
in_dirs = ['/Users/kirstencasey/noise_image_tests/artpop_grids_output_6000000/']
save_dir = '/Users/kirstencasey/noise_image_tests/grid_results_6000000'
measured_fn_id_r = 'blank_mock_injected_Sersic_bestfit-params_*_r.txt'
measured_fn_id_b = 'blank_mock_injected_Sersic_bestfit-params_*_b.txt'
true_fn_id = 'artpop_parameters_*_num-iters*'
sbf_fn_id = 'sbf_results_*_num-iters*_weighted_avg_sbf*'
model_funcs = ['Sersic']
chisq_r = []
chisq_b = []
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
iters = []
fix_true_param_names = ['ellipticity','n','radius', 'theta', 'distance']
if len(model_funcs) > 1 : comp = 'comp_2'
else: comp = 'comp_1'
config_filename = os.path.join('../modeling-config_imfit.yml')
# Open and read config file
with open(config_filename, 'r') as filename:
    config = yaml.load(filename, Loader=yaml.FullLoader)

def sort_key(fn):
    return int(fn.split('_iter')[-1].split('_')[0])

for in_dir in in_dirs:
    # Get true parameters
    true_files = glob.glob(os.path.join(in_dir,true_fn_id))
    temp_results = dict()
    for true_file in true_files:
        file = open(true_file, "rb")
        param = true_file.split('_')[-1]
        if param in fix_true_param_names: param = 'true_'+param
        elif param == 'r' : param = 'true_mag_r'
        elif param == 'b' : param = 'true_mag_b'
        elif param == 'magnitude': param = 'true_sbf_magnitude'
        temp_results[param] = np.load(file,allow_pickle=True)
        file.close

    color_true =  temp_results['true_mag_b'] - temp_results['true_mag_r']
    temp_results['true_color'] = color_true

    # Get measured parameters
    sbf_files = glob.glob(os.path.join(in_dir,sbf_fn_id))
    for sbf_file in sbf_files:
        file = open(sbf_file, "rb")
        param = 'measured_' + sbf_file.split('avg_')[-1]
        temp_results[param] = np.load(file,allow_pickle=True)
        file.close

    measured_files_r = glob.glob(os.path.join(in_dir,measured_fn_id_r))
    measured_files_r.sort(key=sort_key)
    measured_files_b = glob.glob(os.path.join(in_dir,measured_fn_id_b))
    measured_files_b.sort(key=sort_key)
    for i in range(len(measured_files_r)):
        fn_r = measured_files_r[i]
        fn_b = measured_files_b[i]

        bestfit_r = io.read_results(fn_r, model_funcs)
        bestfit_b = io.read_results(fn_b, model_funcs)

        chisq_r.append(bestfit_r['reduced_chisq'])
        chisq_b.append(bestfit_b['reduced_chisq'])
        radii.append(bestfit_r[comp]['r_e'])
        pa.append(bestfit_r[comp]['PA'])
        ell.append(bestfit_r[comp]['ell'])
        n.append(bestfit_r[comp]['n'])
        xpos.append(bestfit_r[comp]['X0'])
        ypos.append(bestfit_r[comp]['Y0'])
        I_e_r.append(bestfit_r[comp]['I_e'])
        I_e_b.append(bestfit_b[comp]['I_e'])

        iters.append(i)

        mag_r, mag_b, color = imfit.summarize_results(config, bestfit_r[comp], bestfit_b[comp])
        mags_r.append(mag_r)
        mags_b.append(mag_b)

    temp_results['chisq_r']=np.asarray(chisq_r)
    temp_results['chisq_b']=np.asarray(chisq_b)
    temp_results['measured_radius']=np.asarray(radii)
    temp_results['measured_pa']=np.asarray(pa)
    temp_results['measured_ellipticity']=np.asarray(ell)
    temp_results['measured_n']=np.asarray(n)
    temp_results['measured_xpos']=np.asarray(xpos) - 1
    temp_results['measured_ypos']=np.asarray(ypos) - 1
    temp_results['measured_I_e_r']=np.asarray(I_e_r)
    temp_results['measured_I_e_b']=np.asarray(I_e_b)
    temp_results['measured_mag_r']=np.asarray(mags_r)
    temp_results['measured_mag_b']=np.asarray(mags_b)
    temp_results['model_num']=np.asarray(iters)
    temp_results['measured_color']=np.asarray(mags_b)-np.asarray(mags_r)

    if in_dir == in_dirs[0]:
        results = temp_results
    else:
        for key in list(temp_results.keys()):
            results[key] = np.concatenate((results[key],temp_results[key]))

# Unit conversions:
results['true_distance'] *= 1e6
results['true_radius_pixels'] = misc.parsecs_to_pixels(results['true_radius'], results['true_distance'], 0.224)
#print(results['measured_radius'])
#print(results['true_radius_pixels'])

results['mag_r_error'] = results['true_mag_r'] - results['measured_mag_r']
results['mag_b_error'] = results['true_mag_b'] - results['measured_mag_b']
results['color_error'] = results['measured_color'] - results['true_color']
results['frac_radius_error'] = (results['measured_radius'] - results['true_radius_pixels']) / results['true_radius_pixels']
results['radius_error'] = results['measured_radius'] - results['true_radius_pixels']
results['ellip_error'] = results['measured_ellipticity'] - results['true_ellipticity']
results['n_error'] = results['measured_n'] - results['true_n']
results['sbf_mag_error'] = results['true_sbf_magnitude'] - results['measured_sbf_mags']
results['dist_error'] = results['measured_sbf_dist_a'] - results['true_distance']
results['frac_dist_error'] = (results['measured_sbf_dist_a'] - results['true_distance']) / results['true_distance']
results['xpos_error'] = results['measured_xpos'] - 525.0
results['ypos_error'] = results['measured_ypos'] - 525.0

print('Finished collecting results.')
# Make plots!!! :
plot_params = ['mag_r_error','mag_b_error','color_error','frac_radius_error','radius_error','ellip_error','n_error','sbf_mag_error','dist_error','frac_dist_error']
labels = ['r-band magnitude','b-band magnitude','color','Fractional radius','Radius (pixels)','Ellipticity','Sersic index','SBF magnitude', 'Distance (pc)', 'Fractional distance']
save_names = ['mag_r','mag_b','color','frac_radius','radius','ellip','n','sbfmag','dist','frac_dist']

# As a function of ArtPop radius
rad300 = misc.filter_dict(results,conditions={'true_radius': 300})
rad600 = misc.filter_dict(results,conditions={'true_radius': 600})
rad900 = misc.filter_dict(results,conditions={'true_radius': 900})
# As a function of ArtPop sersic index
n2 = misc.filter_dict(results,conditions={'true_n': 0.2})
n5 = misc.filter_dict(results,conditions={'true_n': 0.5})
n8 = misc.filter_dict(results,conditions={'true_n': 0.8})
# As a function of ArtPop ellipticity
ell0 = misc.filter_dict(results,conditions={'true_ellipticity': 0.})
ell2 = misc.filter_dict(results,conditions={'true_ellipticity': 0.2})
ell4 = misc.filter_dict(results,conditions={'true_ellipticity': 0.4})
# As a function of stellar pop
pop1 = misc.filter_dict(results,conditions={'feh': -1.2})
pop2 = misc.filter_dict(results,conditions={'feh': -1.0})
pop3 = misc.filter_dict(results,conditions={'feh': -0.65})
age1 = 10.0
age2 = 9.9
age3 = 10.1
#'''
iso1 = MISTIsochrone(age1, -1.2, 'UBVRIplus', mist_path='/Users/kirstencasey/Documents/Resources/MIST')
iso2 = MISTIsochrone(age2, -1.0, 'UBVRIplus', mist_path='/Users/kirstencasey/Documents/Resources/MIST')
iso3 = MISTIsochrone(age3, -0.65, 'UBVRIplus', mist_path='/Users/kirstencasey/Documents/Resources/MIST')


for param,label,save in zip(plot_params,labels,save_names):
    print(f'Plotting {param}')
    # As a function of radius
    for filtered,color in zip([rad300,rad600,rad900],['lightpink','mediumvioletred','firebrick']):
        rad = filtered['true_radius'][0]
        plt.hist(filtered[param],bins='auto',color=color,alpha=0.4, label=f'{rad} pc radius', histtype='stepfilled', ec="k", density=False)
        plt.axvline(x=np.median(filtered[param]), color=color, linestyle='dashed', linewidth=3, label=f'{rad} pc radius median: {round(np.median(filtered[param]),3)}')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_test_radii_{save}.png'))
    plt.clf()

    # As a function of sersic index
    for filtered,color in zip([n2,n5,n8],['lightskyblue','dodgerblue','navy']):
        n = filtered['true_n'][0]
        plt.hist(filtered[param],bins='auto',color=color,alpha=0.4, label=f'{n} sersic index', histtype='stepfilled', ec="k", density=False)
        plt.axvline(x=np.median(filtered[param]), color=color, linestyle='dashed', linewidth=3, label=f'{n} sersic index median: {round(np.median(filtered[param]),3)}')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_test_n_{save}.png'))
    plt.clf()

    # As a function of ellipticity
    for filtered,color in zip([ell0,ell2,ell4],['forestgreen','mediumseagreen','yellowgreen']):
        ell = filtered['true_ellipticity'][0]
        plt.hist(filtered[param],bins='auto',color=color,alpha=0.4, label=f'{ell} ellip', histtype='stepfilled', ec="k", density=False)
        plt.axvline(x=np.median(filtered[param]), color=color, linestyle='dashed', linewidth=3, label=f'{ell} ellip median: {round(np.median(filtered[param]),3)}')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_test_ell_{save}.png'))
    plt.clf()

    # As a function of stellar pop
    for filtered,color,age in zip([pop1,pop2,pop3],['mediumslateblue','rebeccapurple','darkmagenta'],[age1,age2,age3]):
        pop = 'feh ' + str(filtered['feh'][0]) + f' age {age}'
        plt.hist(filtered[param],bins='auto',color=color,alpha=0.4, label=pop, histtype='stepfilled', ec="k", density=False)
        plt.axvline(x=np.median(filtered[param]), color=color, linestyle='dashed', linewidth=3, label=f'{pop} ellip median: {round(np.median(filtered[param]),3)}')
    plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
    plt.xlabel(f'{label} error')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_test_stellarpop_{save}.png'))
    plt.clf()

for pop,iso,age in zip([pop1,pop2,pop3],[iso1,iso2,iso3],[age1,age2,age3]):
    pop_name = 'feh ' + str(pop['feh'][0]) + f' age {age}'
    iso_sbf = np.zeros(len(pop['true_sbf_magnitude']))+ iso.ssp_sbf_mag('Bessell_R') #- 0.190278
    jerjen_sbfmags = (6.09 * pop['true_color']) - 8.81

    plt.scatter(jerjen_sbfmags,iso_sbf,color='steelblue',label='Jerjen empirical relation, `true\' color from ArtPop model')
    iso_color = np.zeros(len(pop['true_color']))+ iso.ssp_color('Bessell_B','Bessell_R') #+ 0.29779
    jerjen_sbfmags = (6.09 * iso_color) - 8.81
    plt.scatter(jerjen_sbfmags,iso_sbf,marker='*',color='tomato',s=200,label='Jerjen empirical relation, MISTIsochrone\'s ssp color')
    plt.scatter(pop['true_sbf_magnitude']-5*np.log10(3.7e6)+5,iso_sbf,color='dimgray',label='`true\' ArtPop SBF mag and color')
    plt.ylabel('MISTIsochrone\'s ssp sbf mag')
    plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
    plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_test_stellarpop_compare_sbfmags_{pop_name}_mistiso.png'))
    plt.clf()

    jerjen_sbfmags = (6.09 * pop['true_color']) - 8.81
    artpop_sbfmags = pop['true_sbf_magnitude']-5*np.log10(3.7e6)+5
    hline = (6.09 * (iso.ssp_color('Bessell_B','Bessell_R'))) - 8.81
    vline = iso.ssp_sbf_mag('Bessell_R')
    plt.scatter(artpop_sbfmags,jerjen_sbfmags,color='dimgray')
    plt.ylabel('Jerjen SBF mag using ArtPop\'s `true\' color')
    plt.xlabel('SBF mag using ArtPop\'s `true\' apparent SBF and distance')
    plt.vlines(vline,min(jerjen_sbfmags),max(jerjen_sbfmags),label='MISTIsochrone\'s ssp sbf mag',color='steelblue')
    plt.hlines(hline,min(artpop_sbfmags),max(artpop_sbfmags),color='tomato',label='MISTIsochrone\'s ssp color')
    plt.legend(prop={'size': 10},loc='lower left')
    plt.title(pop_name+f' | Jerjen - MIST SBF mag : {round(hline - vline,3)}')
    plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_test_stellarpop_compare_sbfmags_{pop_name}_jerjen_artpop.png'))
    plt.clf()

for pop,iso,age,color in zip([pop1,pop2,pop3],[iso1,iso2,iso3],[age1,age2,age3],['mediumslateblue','rebeccapurple','darkmagenta']):
    label = 'feh ' + str(pop['feh'][0]) + f' age {age}'
    dist_mod_mist = pop['measured_sbf_mags'] - iso.ssp_sbf_mag('Bessell_R')
    d_a = 10**((dist_mod_mist + 5)/5) # parsecs
    true_dist = pop['true_distance']
    plt.hist((d_a - pop['true_distance'])/1e3,bins='auto',color=color,alpha=0.4, label=label, histtype='stepfilled', ec="k", density=False)
    plt.axvline(x=np.median((d_a - true_dist)/1e3), color=color, linestyle='dashed', linewidth=3, label=f'{label} median: {round(np.median((d_a - true_dist)/1e3),3)}')

plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
plt.xlabel(f'Distance error (kpc)')
plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_test_stellarpop_dist_mist.png'))
plt.clf()
#'''
# Examine outliers:
color_outliers = misc.filter_dict(results,conditions={'color_error': -0.015},condition_type={'color_error':'less than'})
ellip_outliers = misc.filter_dict(results,conditions={'ellip_error': -0.1},condition_type={'ellip_error':'less than'})
radius_outliers, radius_mask = misc.filter_dict(results,conditions={'radius_error': -10},condition_type={'radius_error':'less than'}, return_mask=True)
b_mag_outliers = misc.filter_dict(results,conditions={'mag_b_error': -0.03},condition_type={'mag_b_error':'less than'})
r_mag_outliers = misc.filter_dict(results,conditions={'mag_r_error': -0.03},condition_type={'mag_r_error':'less than'})
n_outliers = misc.filter_dict(results,conditions={'n_error': -0.03},condition_type={'n_error':'less than'})

b_mag_outliers = misc.merge_dicts([b_mag_outliers, misc.filter_dict(results,conditions={'mag_b_error': 0.03},condition_type={'mag_b_error':'greater than'})])
r_mag_outliers = misc.merge_dicts([r_mag_outliers, misc.filter_dict(results,conditions={'mag_r_error': 0.03},condition_type={'mag_r_error':'greater than'})])
n_outliers = misc.merge_dicts([n_outliers, misc.filter_dict(results,conditions={'n_error': 0.03},condition_type={'n_error':'greater than'})])
#'''
for outliers,outlier_type in zip([color_outliers,ellip_outliers,radius_outliers,b_mag_outliers,r_mag_outliers,n_outliers],['color','ellip','radius','bmag','rmag','n']):
    pop1 = misc.filter_dict(outliers,conditions={'feh': -1.2})
    pop2 = misc.filter_dict(outliers,conditions={'feh': -1.0})
    pop3 = misc.filter_dict(outliers,conditions={'feh': -0.65})

    for param,label,save in zip(plot_params,labels,save_names):
        for filtered,color,age,feh in zip([pop1,pop2,pop3],['mediumslateblue','rebeccapurple','darkmagenta'],[age1,age2,age3],[-1.2,-1.0,-0.65]):
            #print(age, filtered)
            pop = f'feh {feh} age {age}'
            plt.hist(filtered[param],bins='auto',color=color,alpha=0.4, label=pop, histtype='stepfilled', ec="k", density=False)
            plt.axvline(x=np.median(filtered[param]), color=color, linestyle='dashed', linewidth=3, label=f'{pop} ellip median: {round(np.median(filtered[param]),3)}')
        plt.axvline(x=0, color='black', linestyle='dashed', linewidth=3)
        plt.xlabel(f'{label} error')
        plt.title(f'{outlier_type} outliers')
        plt.legend(prop={'size': 15},bbox_to_anchor=(1.05, 1),loc='upper left')
        plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_{outlier_type}_outliers_stellarpop_{save}.png'))
        plt.clf()
#'''
#'''
# Ellip outliers
plt.scatter(ellip_outliers['radius_error'], ellip_outliers['mag_b_error'], c=ellip_outliers['ellip_error'], cmap='inferno')
plt.xlabel('Radius error (pixels)')
plt.ylabel('b-band magnitude error')
plt.title('Ellipticity outliers')
plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=min(ellip_outliers['ellip_error']), vmax=max(ellip_outliers['ellip_error']), clip=False),cmap='inferno'),label='Ellipticity error')
plt.hlines(0,min(ellip_outliers['radius_error']),max(ellip_outliers['radius_error']),color='black',linestyle='dashed')
plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_ellip_outliers_radius_bmag.png'))
plt.clf()

plt.scatter(ellip_outliers['radius_error'], ellip_outliers['n_error'], c=ellip_outliers['ellip_error'], cmap='inferno')
plt.xlabel('Radius error (pixels)')
plt.ylabel('Sersic index error')
plt.title('Ellipticity outliers')
plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=min(ellip_outliers['ellip_error']), vmax=max(ellip_outliers['ellip_error']), clip=False),cmap='inferno'),label='Ellipticity error')
plt.hlines(0,min(ellip_outliers['radius_error']),max(ellip_outliers['radius_error']),color='black',linestyle='dashed')
plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_ellip_outliers_radius_n.png'))
plt.clf()

# Radius outliers
plt.scatter(radius_outliers['mag_b_error'], radius_outliers['n_error'], c=radius_outliers['frac_radius_error'], cmap='inferno')
plt.xlabel('b-band magnitude error')
plt.ylabel('Sersic index error')
plt.title('Radius outliers')
plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=min(radius_outliers['frac_radius_error']), vmax=max(radius_outliers['frac_radius_error']), clip=False),cmap='inferno'),label='Fractional radius error (pixels)')
#plt.hlines(0,min(ellip_outliers['radius_error']),max(ellip_outliers['radius_error']),color='black',linestyle='dashed')
plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_radius_outliers_bmag_n.png'))
plt.clf()

plt.scatter(radius_outliers['ellip_error'], radius_outliers['n_error'], c=radius_outliers['frac_radius_error'], cmap='inferno')
plt.xlabel('Ellipticity error')
plt.ylabel('Sersic index error')
plt.title('Radius outliers')
plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=min(radius_outliers['frac_radius_error']), vmax=max(radius_outliers['frac_radius_error']), clip=False),cmap='inferno'),label='Fractional radius error (pixels)')
plt.hlines(0,min(radius_outliers['ellip_error']),max(radius_outliers['ellip_error']),color='black',linestyle='dashed')
plt.vlines(0,min(radius_outliers['n_error']),max(radius_outliers['n_error']),color='black',linestyle='dashed')
plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_radius_outliers_ellip_n.png'))
plt.clf()

x = np.linspace(-.23,0.08)
plt.scatter(radius_outliers['mag_b_error'], radius_outliers['mag_r_error'], c=radius_outliers['frac_radius_error'], cmap='inferno')
plt.plot(x,x,color='black')
plt.xlabel('b-band magnitude error')
plt.ylabel('r-band magnitude error')
plt.title('Radius outliers')
plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=min(radius_outliers['frac_radius_error']), vmax=max(radius_outliers['frac_radius_error']), clip=False),cmap='inferno'),label='Fractional radius error')
plt.hlines(0,min(radius_outliers['mag_b_error']),max(radius_outliers['mag_b_error']),color='black',linestyle='dashed')
plt.vlines(0,min(radius_outliers['mag_r_error']),max(radius_outliers['mag_r_error']),color='black',linestyle='dashed')
plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_radius_outliers_bmag_rmag.png'))
plt.clf()

plt.scatter(radius_outliers['mag_r_error'], radius_outliers['sbf_mag_error'], c=radius_outliers['frac_radius_error'], cmap='inferno')
plt.xlabel('r-band magnitude error')
plt.ylabel('sbf magnitude error')
plt.title('Radius outliers')
plt.colorbar(plt.cm.ScalarMappable(norm=colors.Normalize(vmin=min(radius_outliers['frac_radius_error']), vmax=max(radius_outliers['frac_radius_error']), clip=False),cmap='inferno'),label='Fractional radius error')
plt.hlines(0,min(radius_outliers['mag_r_error']),max(radius_outliers['mag_r_error']),color='black',linestyle='dashed')
plt.vlines(0,min(radius_outliers['sbf_mag_error']),max(radius_outliers['sbf_mag_error']),color='black',linestyle='dashed')
plt.savefig(os.path.join(save_dir,f'artpop_grid_test_results_radius_outliers_rmag_sbfmag.png'))
plt.clf()
#'''
print(np.where(radius_mask==True))
