import lbcred, os, yaml, artpop
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.collections import PatchCollection

# Acceptable color ranges
carlsten_gi_range_inBR = [0.7,1.2] # B-R | [0.3,0.8] # g-i old
carlsten_gr_range_inBR = [0.8,1.3] # B-R | [0.3,0.6]  # g-r 
cantiello_gi_range_inBR = [1.2,1.5] # B-R | [0.8,1.05] # g-i old
acceptable_color = 'navy'
unacceptable_color = 'darkred'

def get_sbf_distance(config_fn, model_feh=None) : 
    
    with open(config_fn, 'r') as filename:
        config = yaml.load(filename, Loader=yaml.FullLoader)
        if model_feh != None: config['get_distance']['feh_model'] = model_feh
    
    metals = np.linspace(config['get_distance']['feh_min'], config['get_distance']['feh_max'], num=config['get_distance']['num_isochrones'])
    ages = np.linspace(config['get_distance']['age_min'], config['get_distance']['age_max'], num=config['get_distance']['num_isochrones'])
        
    closest_metal = metals[np.where(abs(metals-config['get_distance']['feh_model'])==min(abs(metals-config['get_distance']['feh_model'])))[0]][0]
    print(f'\nUsing closest target feh: {round(closest_metal,3)}\n')
    
    mag2_targetgal = np.array([15.96])
    mag1_targetgal = np.array([14.76])
    sbf_measured_targetgal = np.array([26.34])
    sbf_uncertainty_targetgal = np.array([0.07])
    if config['get_distance']['save_id'] is not None: save_id = config['get_distance']['save_id']
    else: save_id = ''
    if config['get_distance']['correct_for_bias']:
        sbf_bias = config['get_distance']['sbf_mag_bias']
        color_bias = config['get_distance']['color_bias']
    else:
        sbf_bias = 0.
        color_bias = 0.

    color_targetgal = mag2_targetgal - mag1_targetgal
    color_targetgal_inrange = []
    for color_range in [carlsten_gi_range_inBR,carlsten_gr_range_inBR,cantiello_gi_range_inBR]:
        if color_targetgal[0]+color_bias < color_range[0] or color_targetgal[0]+color_bias > color_range[1] : color_targetgal_inrange.append(unacceptable_color)
        else: color_targetgal_inrange.append(acceptable_color)
    color_targetgal_inrange = np.asarray(color_targetgal_inrange)

    metal_arrs = []
    for metal in metals: 
        metal_arr = np.zeros(config['get_distance']['num_isochrones']) + metal
        metal_arrs.append(metal_arr)
        
    print('Getting isochrones...')
    sbf_mags_mist_arrs = []
    colors_mist_arrs = []
    sbf_mags_cfht_ginew_arrs = []
    colors_cfht_ginew_arrs = []
    sbf_mags_cfht_giold_arrs = []
    colors_cfht_giold_arrs = []
    sbf_mags_cfht_gr_arrs = []
    colors_cfht_gr_arrs = []
    
    for metal_arr in metal_arrs:
        sbf_mags_mist = []
        colors_mist = []
        sbf_mags_cfht_ginew = []
        colors_cfht_ginew = []
        sbf_mags_cfht_giold = []
        colors_cfht_giold = []
        sbf_mags_cfht_gr = []
        colors_cfht_gr = []
        print(f'Working on feh = {metal_arr[0]}')
        for metal,age in zip(metal_arr,ages):
            iso_mist = artpop.stars.MISTIsochrone(age, metal, 'UBVRIplus', mist_path=config['get_distance']['mist_path'])
            sbf_mags_mist.append(iso_mist.ssp_sbf_mag('Bessell_R') - 0.190278)
            colors_mist.append(iso_mist.ssp_color('Bessell_B','Bessell_R') + 0.107512 + 0.190278)
            
            iso_cfht = artpop.stars.MISTIsochrone(age, metal, 'CFHTugriz', mist_path=config['get_distance']['mist_path'])
            
            sbf_mags_cfht_ginew.append(iso_cfht.ssp_sbf_mag('CFHT_i_new'))
            colors_cfht_ginew.append(iso_cfht.ssp_color('CFHT_g','CFHT_i_new'))
            
            sbf_mags_cfht_giold.append(iso_cfht.ssp_sbf_mag('CFHT_i_old'))
            colors_cfht_giold.append(iso_cfht.ssp_color('CFHT_g','CFHT_i_old'))
            
            sbf_mags_cfht_gr.append(iso_cfht.ssp_sbf_mag('CFHT_r'))
            colors_cfht_gr.append(iso_cfht.ssp_color('CFHT_g','CFHT_r'))

        sbf_mags_mist_arrs.append(np.asarray(sbf_mags_mist))
        colors_mist_arrs.append(np.asarray(colors_mist)) 
        sbf_mags_cfht_ginew_arrs.append(np.asarray(sbf_mags_cfht_ginew))
        colors_cfht_ginew_arrs.append(np.asarray(colors_cfht_ginew)) 
        sbf_mags_cfht_giold_arrs.append(np.asarray(sbf_mags_cfht_giold))
        colors_cfht_giold_arrs.append(np.asarray(colors_cfht_giold))
        sbf_mags_cfht_gr_arrs.append(np.asarray(sbf_mags_cfht_gr))
        colors_cfht_gr_arrs.append(np.asarray(colors_cfht_gr))
        
    
    ## Get conversions between filters 
    print('Converting between filters...')
    
    # Find conversion between B-R (Bessell/Vega) and g-i (CFHT/AB)
    br_gi_new_conversion = np.asarray(colors_mist_arrs) - np.asarray(colors_cfht_ginew_arrs)
    br_gi_old_conversion = np.asarray(colors_mist_arrs) - np.asarray(colors_cfht_giold_arrs)
    
    # Find conversion between R-band (Bessell/Vega) and i-band (CFHT/AB) SBF magnitudes
    r_i_ginew_conversion = np.asarray(sbf_mags_mist_arrs) - np.asarray(sbf_mags_cfht_ginew_arrs)
    r_i_giold_conversion = np.asarray(sbf_mags_mist_arrs) - np.asarray(sbf_mags_cfht_giold_arrs)
    
    # For each sp in Bessell MIST isos, calculate the g-i color and Mi SBF magnitude using the conversions, Carlsten relation
    gi_new_calculated = np.asarray(colors_mist_arrs) - br_gi_new_conversion
    gi_old_calculated = np.asarray(colors_mist_arrs) - br_gi_old_conversion
    
    Mi_ginew_calculated = -3.17+2.15*gi_new_calculated
    Mi_giold_calculated = -3.17+2.15*gi_old_calculated
    
    Mi_giold_cantiello = -0.93+3.25*(gi_old_calculated-0.95)
    
    # Convert Mi to Mr using conversion
    Mr_ginew_calculated = Mi_ginew_calculated + r_i_ginew_conversion
    Mr_giold_calculated = Mi_giold_calculated + r_i_giold_conversion
    Mr_giold_cantiello = Mi_giold_cantiello + r_i_giold_conversion
    
    # --------------
    
    # Find conversion between B-R (Bessell/Vega) and g-r (CFHT/AB)
    BR_gr_conversion = np.asarray(colors_mist_arrs) - np.asarray(colors_cfht_gr_arrs)
    
    # Find conversion between R-band (Bessell/Vega) and r-band (CFHT/AB) SBF magnitudes
    R_r_gr_conversion = np.asarray(sbf_mags_mist_arrs) - np.asarray(sbf_mags_cfht_gr_arrs)
    
    # For each sp in Bessell MIST isos, calculate the g-r color and Mi SBF magnitude using the conversions, Carlsten relation
    gr_calculated = np.asarray(colors_mist_arrs) - BR_gr_conversion
    Mr_gr_calculated = 4.21*gr_calculated - 3.0
    
    # Convert Mi to Mr using conversion
    MR_gr_calculated = Mr_gr_calculated + R_r_gr_conversion

    ## Get results
    # Combine colors and SBF mags 
    norm = lambda x: x/len(metals)
    cmap = plt.cm.get_cmap('inferno')
    min_mag = 999
    max_mag = -999
    if config['get_distance']['plot_mag_color_relation']: 
        print('Plotting SBF magnitude--color relation...')
        fig,ax = plt.subplots(1,1,figsize=(15,10))
    for i in range(len(metal_arrs)):
        sbf_mags_mist = sbf_mags_mist_arrs[i]
        colors_mist = colors_mist_arrs[i]
        sbf_mags_carlsten_ginew = Mr_ginew_calculated[i]
        sbf_mags_carlsten_giold = Mr_giold_calculated[i]
        sbf_mags_carlsten_gr = MR_gr_calculated[i]
        sbf_mags_cantiello_giold = Mr_giold_cantiello[i]
    
        # Fits for SBF mags as a function of color
        fit_carlsten_ginew = np.polyfit(colors_mist,sbf_mags_carlsten_ginew, 1)
        p_carlsten_ginew = np.poly1d(fit_carlsten_ginew)
        fit_carlsten_giold = np.polyfit(colors_mist,sbf_mags_carlsten_giold, 1)
        p_carlsten_giold = np.poly1d(fit_carlsten_giold)
        
        fit_carlsten_gr = np.polyfit(colors_mist,sbf_mags_carlsten_gr, 1)
        p_carlsten_gr = np.poly1d(fit_carlsten_gr)
        fit_carlsten_gr = np.polyfit(colors_mist,sbf_mags_carlsten_gr, 1)
        p_carlsten_gr = np.poly1d(fit_carlsten_gr)
        
        fit_cantiello_giold = np.polyfit(colors_mist,sbf_mags_cantiello_giold, 1)
        p_cantiello_giold = np.poly1d(fit_cantiello_giold)
        
        fit_mist = np.polyfit(colors_mist,sbf_mags_mist, 1)
        p_mist = np.poly1d(fit_mist)
    
        # Plot
        if config['get_distance']['plot_mag_color_relation']: 
            ax.scatter(colors_mist,sbf_mags_mist,color=cmap(norm(i+1)),marker='o')
            ax.plot(colors_mist,p_mist(np.asarray(colors_mist)),'-',label=f'MIST; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
        
            ax.scatter(colors_mist,sbf_mags_carlsten_ginew,color=cmap(norm(i+1)),marker='^')
            ax.plot(colors_mist,p_carlsten_ginew(np.asarray(colors_mist)),'-.',label=f'cfht $g-i$ new; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
        
            ax.scatter(colors_mist,sbf_mags_carlsten_giold,color=cmap(norm(i+1)),marker='v')
            ax.plot(colors_mist,p_carlsten_giold(np.asarray(colors_mist)),'-.',label=f'cfht $g-i$ old; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
        
            ax.scatter(colors_mist,sbf_mags_cantiello_giold,color=cmap(norm(i+1)),marker='x')
            ax.plot(colors_mist,p_cantiello_giold(np.asarray(colors_mist)),':',label=f'cantiello; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
            
            ax.scatter(colors_mist,sbf_mags_carlsten_gr,color=cmap(norm(i+1)),marker='*')
            ax.plot(colors_mist,p_carlsten_gr(np.asarray(colors_mist)),'-.',label=f'cfht $g-r$; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
            
            if min(sbf_mags_cantiello_giold) < min_mag: min_mag = min(sbf_mags_cantiello_giold) 
            if max(sbf_mags_cantiello_giold) > max_mag: max_mag = max(sbf_mags_cantiello_giold)

        if metals[i] == closest_metal:
            save_p_mist = p_mist
            save_p_carlsten_ginew = p_carlsten_ginew
            save_p_carlsten_giold = p_carlsten_giold
            save_p_cantiello_giold = p_cantiello_giold
            save_p_carlsten_gr = p_carlsten_gr

    if config['get_distance']['plot_mag_color_relation']: 
        ax.set_title('SBF Magnitude -- Color Calibrations (Vega)')
        ax.set_xlabel('ArtPop model $B-R$ Color')
        ax.set_ylabel('Absolute SBF Magnitude ($R$-band)')
        ax.set_xlim(0.7)
    
        legend_elements = [Line2D([0], [0], marker='o', color='black', ls='-', label='MIST relation',markerfacecolor='black', markersize=10),
                           Line2D([0], [0], marker='^', color='black', ls='-.', label='Carlsten $g-i$ (new) relation',markerfacecolor='black', markersize=10),
                           Line2D([0], [0], marker='v', color='black', ls='-.', label='Carlsten $g-i$ (old) relation',markerfacecolor='black', markersize=10),
                           Line2D([0], [0], marker='*', color='black', ls='-.', label='Carlsten $g-r$ relation',markerfacecolor='black', markersize=10),
                           Line2D([0], [0], marker='x', color='black', ls=':', label='Cantiello $g-i$ (old) relation',markerfacecolor='black', markersize=10)]
        i=0
        for metal in metals:
            legend_elements.append(Patch(facecolor=cmap(norm(i+1)),label=f'metallicity: {round(metal,2)}'))
            i+=1
        ax.set_ylim(min_mag-0.1,max_mag+0.1)
        ax.vlines(color_targetgal[0],min_mag-0.2,max_mag+0.2,colors='black',linestyles='--')
        plt.legend(prop={'size': 15},handles=legend_elements)#,bbox_to_anchor=(1.05, 1),loc='upper left'
        plt.savefig(os.path.join(config['out_dir'],f'color_to_absSBFmag_{save_id}.png'))
        plt.clf()
    
    
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------

    median_dists_targetgal = []
    dist_mods_targetgal = []
    median_dist_errors_targetgal_syst = None
    
    for abs_sbf_fit,calibration_id in zip([save_p_mist,save_p_carlsten_ginew,save_p_carlsten_giold,save_p_carlsten_gr,save_p_cantiello_giold],['mist','carlsten_ginew','carlsten_giold','carlsten_gr','cantiello_giold']):
        print('\n\nCalculating distance using: ',calibration_id)
        #### TargetGal SBF
        dist_mod_targetgal = sbf_measured_targetgal + sbf_bias - abs_sbf_fit(color_targetgal+color_bias)
        print(f'{calibration_id} absolute magnitude: ',abs_sbf_fit(color_targetgal+color_bias))
        dist_mods_targetgal.append(dist_mod_targetgal[0])
        dist_targetgal = 10**((dist_mod_targetgal + 5)/5) # parsecs
        print(f'Measured distance: {round(dist_targetgal[0]/1e6,2)} Mpc')
        
        # Calculate maximum error - systematic
        if calibration_id == 'cantiello_giold':
            dist_mod_max_syst = sbf_measured_targetgal + sbf_bias - abs_sbf_fit(color_targetgal+color_bias) + 0.12
        elif calibration_id == 'mist':
            dist_mod_max_syst = dist_mod_targetgal
        else:
            dist_mod_max_syst = sbf_measured_targetgal + sbf_bias - abs_sbf_fit(color_targetgal+color_bias) + 0.26
    
        dist_mod_max_error_syst = abs(dist_mod_max_syst - dist_mod_targetgal)
        
        # Calculate maximum error - statistical
        if calibration_id == 'mist': 
            dist_mod_max_stat = dist_mod_targetgal
        else: 
            dist_mod_max_stat = sbf_measured_targetgal + sbf_bias + sbf_uncertainty_targetgal + config['get_distance']['additional_sbf_mag_uncertainty'] - abs_sbf_fit(color_targetgal + color_bias - config['get_distance']['additional_color_uncertainty'])
        dist_mod_max_error_stat = abs(dist_mod_max_stat - dist_mod_targetgal)
        
        
        # Convert distance modulus errors to distance errors
        dist_err_stat_lower = dist_targetgal - 10**((dist_mod_targetgal-dist_mod_max_error_stat + 5)/5) 
        dist_err_stat_upper = 10**((dist_mod_targetgal+dist_mod_max_error_stat + 5)/5) - dist_targetgal
        dist_err_syst_lower = dist_targetgal - 10**((dist_mod_targetgal-dist_mod_max_error_syst + 5)/5) 
        dist_err_syst_upper = 10**((dist_mod_targetgal+dist_mod_max_error_syst + 5)/5) - dist_targetgal
        print(f'Systematic distance error (lower): {round(dist_err_syst_lower[0]/1e3,2)} kpc')
        print(f'Systematic distance error (upper): {round(dist_err_syst_upper[0]/1e3,2)} kpc')
        print(f'Statistical distance error (lower): {round(dist_err_stat_lower[0]/1e3,2)} kpc')
        print(f'Statistical distance error (upper): {round(dist_err_stat_upper[0]/1e3,2)} kpc')
        print('------------------------------------------------------------------')
        if median_dist_errors_targetgal_syst is None and calibration_id != 'mist' and calibration_id != 'carlsten_ginew': 
            median_dist_errors_targetgal_stat = np.reshape([dist_err_stat_lower[0],dist_err_stat_upper[0]],(2,-1))
            median_dist_errors_targetgal_syst = np.reshape([dist_err_syst_lower[0],dist_err_syst_upper[0]],(2,-1))
        elif median_dist_errors_targetgal_syst is not None:
            median_dist_errors_targetgal_stat = np.concatenate((median_dist_errors_targetgal_stat,np.reshape([dist_err_stat_lower[0],dist_err_stat_upper[0]],(2,-1))),axis=1)
            median_dist_errors_targetgal_syst = np.concatenate((median_dist_errors_targetgal_syst,np.reshape([dist_err_syst_lower[0],dist_err_syst_upper[0]],(2,-1))),axis=1)
        median_dists_targetgal.append(dist_targetgal[0])

    if config['get_distance']['plot_distances']:
        print('Plotting distance...')
        # Get only the relations we want to plot
        median_dists_targetgal_mist = median_dists_targetgal[0]
        print('MIST: ',median_dists_targetgal_mist)
        median_dists_targetgal = median_dists_targetgal[2:]
        print('others: ',median_dists_targetgal)
        plt.clf()
        fig,ax = plt.subplots(1,1,figsize=(10,10))
        
        good_idxs = [i for i in range(len(color_targetgal_inrange)) if color_targetgal_inrange[i] == acceptable_color]
        bad_idxs = [i for i in range(len(color_targetgal_inrange)) if color_targetgal_inrange[i] == unacceptable_color]
        
        x=np.array([0,1,2])
        my_ticks = ['Carlsten ($g-i$)', 'Carlsten ($g-r$)', 'Cantiello ($g-i$)']  # (uses g-i old)
        ax.set_xticks(x)
        ax.set_xticklabels(labels=my_ticks,rotation=-15)
        ax.set_xlim(-0.5,2.5)
        
        # Systematic uncertainties
        if len(good_idxs) != 0 : 
            ax.errorbar(x[good_idxs],np.asarray(median_dists_targetgal)[good_idxs]/1e6,yerr=np.asarray(median_dist_errors_targetgal_syst)[0:,good_idxs]/1e6,xerr=None,fmt='o',markersize=24,alpha=0.8,capsize=5.0,capthick=3.,c=acceptable_color)
        if len(bad_idxs) != 0: 
            ax.errorbar(x[bad_idxs],np.asarray(median_dists_targetgal)[bad_idxs]/1e6,yerr=np.asarray(median_dist_errors_targetgal_syst)[0:,bad_idxs]/1e6,xerr=None,fmt='o',markersize=24,alpha=0.8,capsize=5.0,capthick=3.,c=unacceptable_color)
        
        # Statistical uncertainties
        if len(good_idxs) != 0 : 
            ax.errorbar(x[good_idxs],np.asarray(median_dists_targetgal)[good_idxs]/1e6,yerr=np.asarray(median_dist_errors_targetgal_stat)[0:,good_idxs]/1e6,xerr=None,fmt='o',markersize=24,alpha=0.8,capsize=10.0,capthick=5.,c=acceptable_color)
        if len(bad_idxs) != 0: 
            ax.errorbar(x[bad_idxs],np.asarray(median_dists_targetgal)[bad_idxs]/1e6,yerr=np.asarray(median_dist_errors_targetgal_stat)[0:,bad_idxs]/1e6,xerr=None,fmt='o',markersize=24,alpha=0.8,capsize=10.0,capthick=5.,c=unacceptable_color)
        
        ax.axhline(median_dists_targetgal_mist/1e6,color='k',ls='--',label='MIST ($B-R$) prediction')
        ax.legend(bbox_to_anchor=(0, 1.1),loc='upper left',fontsize=20)
        ax.set_ylabel('Distance (Mpc)',fontsize=24)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(24)
        plt.savefig(os.path.join(config['out_dir'],f'measured_distance_{save_id}.png'),dpi=500);

    return


#####################################################################################

# Run 'get_sbf_distance' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Measure distance from SBF magnitude results.')
    default_config = os.path.join('/users/PCON0003/osu10713/Blobby_results', 'modeling-config_sbf.yml')
    parser.add_argument('--config_fn', type=str, default=default_config, help='path of the SBF .yml config file' )
    parser.add_argument('--feh', type=float)
    args = parser.parse_args()

    get_sbf_distance(config_fn=args.config_fn, model_feh=args.feh)