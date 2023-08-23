import lbcred, os, yaml, artpop
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.collections import PatchCollection


def get_sbf_distance(in_dir, out_dir, mist_path, model_feh) : 

    num_isochrones = 10
    metals = np.linspace(-2.5, -0.5, num=num_isochrones)
    ages = np.linspace(8.8, 10.2, num=num_isochrones)
    
    closest_metal = metals[np.where(abs(metals-model_feh)==min(abs(metals-model_feh)))[0]]
    
    mag2_targetgal = np.load(open(os.path.join(in_dir,'background_subtraction_test_results_measured_mags_b'), "rb"))
    mag1_targetgal = np.load(open(os.path.join(in_dir,'background_subtraction_test_results_measured_mags_r'), "rb"))
    sbf_measured_targetgal = np.load(open(os.path.join(in_dir,'background_subtraction_test_results_measured_sbfmags'), "rb"))
    sbf_uncertainty_targetgal = np.load(open(os.path.join(in_dir,'background_subtraction_test_results_sbf_mag_uncertainty'), "rb"))
    
    color_targetgal = mag2_targetgal - mag1_targetgal
    
    metal_arrs = []
    for metal in metals: 
        metal_arr = np.zeros(num_isochrones) + metal
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
        print(metal_arr)
        for metal,age in zip(metal_arr,ages):
            iso_mist = artpop.stars.MISTIsochrone(age, metal, 'UBVRIplus', mist_path=mist_path)
            sbf_mags_mist.append(iso_mist.ssp_sbf_mag('Bessell_R'))
            colors_mist.append(iso_mist.ssp_color('Bessell_B','Bessell_R'))
            
            iso_cfht_ginew = artpop.stars.MISTIsochrone(age, metal, 'CFHTugriz', mist_path=mist_path)
            sbf_mags_cfht_ginew.append(iso_cfht_ginew.ssp_sbf_mag('CFHT_i_new'))
            colors_cfht_ginew.append(iso_cfht_ginew.ssp_color('CFHT_g','CFHT_i_new'))
            
            iso_cfht_giold = artpop.stars.MISTIsochrone(age, metal, 'CFHTugriz', mist_path=mist_path)
            sbf_mags_cfht_giold.append(iso_cfht_giold.ssp_sbf_mag('CFHT_i_old'))
            colors_cfht_giold.append(iso_cfht_giold.ssp_color('CFHT_g','CFHT_i_old'))
            
            iso_cfht_gr = artpop.stars.MISTIsochrone(age, metal, 'CFHTugriz', mist_path=mist_path)
            sbf_mags_cfht_gr.append(iso_cfht_gr.ssp_sbf_mag('CFHT_r'))
            colors_cfht_gr.append(iso_cfht_gr.ssp_color('CFHT_g','CFHT_r'))
    
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
    
    # For each sp in Bessell MIST isos, calculate the g-i color and Mi SBF magnitude using the conversions, Carleston relation
    gi_new_calculated = np.asarray(colors_mist_arrs) - br_gi_new_conversion
    gi_old_calculated = np.asarray(colors_mist_arrs) - br_gi_old_conversion
    
    #Mi_ginew_calculated = -3.48+2.48*gi_new_calculated
    #Mi_giold_calculated = -3.48+2.48*gi_old_calculated
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
    
    # For each sp in Bessell MIST isos, calculate the g-r color and Mi SBF magnitude using the conversions, Carleston relation
    gr_calculated = np.asarray(colors_mist_arrs) - BR_gr_conversion
    Mr_gr_calculated = 4.21*gr_calculated - 3.0
    
    # Convert Mi to Mr using conversion
    MR_gr_calculated = Mr_gr_calculated + R_r_gr_conversion
    
    
    ## Get results
    print('Plotting results...')
    # Combine colors and SBF mags - Carleston g-i relation
    norm = lambda x: x/len(metals)
    cmap = plt.cm.get_cmap('inferno')
    
    fig,ax = plt.subplots(1,1,figsize=(20,10))
    for i in range(len(metal_arrs)):
        sbf_mags_mist = sbf_mags_mist_arrs[i]
        colors_mist = colors_mist_arrs[i]
        sbf_mags_carleston_ginew = Mr_ginew_calculated[i]
        sbf_mags_carleston_giold = Mr_giold_calculated[i]
        sbf_mags_carleston_gr = MR_gr_calculated[i]
        sbf_mags_cantiello_giold = Mr_giold_cantiello[i]
    
        # Fits for SBF mags as a function of color
        fit_carleston_ginew = np.polyfit(colors_mist,sbf_mags_carleston_ginew, 1)
        p_carleston_ginew = np.poly1d(fit_carleston_ginew)
        fit_carleston_giold = np.polyfit(colors_mist,sbf_mags_carleston_giold, 1)
        p_carleston_giold = np.poly1d(fit_carleston_giold)
        
        fit_carleston_gr = np.polyfit(colors_mist,sbf_mags_carleston_gr, 1)
        p_carleston_gr = np.poly1d(fit_carleston_gr)
        fit_carleston_gr = np.polyfit(colors_mist,sbf_mags_carleston_gr, 1)
        p_carleston_gr = np.poly1d(fit_carleston_gr)
        
        fit_cantiello_giold = np.polyfit(colors_mist,sbf_mags_cantiello_giold, 1)
        p_cantiello_giold = np.poly1d(fit_cantiello_giold)
        
        fit_mist = np.polyfit(colors_mist,sbf_mags_mist, 1)
        p_mist = np.poly1d(fit_mist)
    
        # Plot
        ax.scatter(colors_mist,sbf_mags_mist,color=cmap(norm(i+1)),marker='o')
        ax.plot(colors_mist,p_mist(np.asarray(colors_mist)),'--',label=f'MIST; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
    
        ax.scatter(colors_mist,sbf_mags_carleston_ginew,color=cmap(norm(i+1)),marker='^')
        ax.plot(colors_mist,p_carleston_ginew(np.asarray(colors_mist)),'-.',label=f'cfht $g-i$ new; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
    
        ax.scatter(colors_mist,sbf_mags_carleston_giold,color=cmap(norm(i+1)),marker='v')
        ax.plot(colors_mist,p_carleston_giold(np.asarray(colors_mist)),'-.',label=f'cfht $g-i$ old; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
    
        ax.scatter(colors_mist,sbf_mags_cantiello_giold,color=cmap(norm(i+1)),marker='x')
        ax.plot(colors_mist,p_cantiello_giold(np.asarray(colors_mist)),'-.',label=f'cantiello; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)
        
        ax.scatter(colors_mist,sbf_mags_carleston_gr,color=cmap(norm(i+1)),marker='*')
        ax.plot(colors_mist,p_carleston_gr(np.asarray(colors_mist)),'-.',label=f'cfht $g-r$; metallicity: {metals[i]}',color=cmap(norm(i+1)),alpha=0.7)

        if metals[i] == closest_metal:
            save_p_mist = p_mist
            save_p_carleston_ginew = p_carleston_ginew
            save_p_carleston_giold = p_carleston_giold
            save_p_cantiello_giold = p_cantiello_giold
            save_p_carleston_gr = p_carleston_gr

    
    ax.set_title('SBF Magnitudes (Vega)')
    ax.set_xlabel('ArtPop model $B-R$ Color')
    ax.set_ylabel('Absolute SBF Magnitude ($R$-band)')
    
    #ax[1].set_title('Measured SBF Magnitudes (Vega)')
    #ax[1].set_xlabel('Measured $B-R$ Color')
    #ax[1].set_ylabel('Absolute SBF Magnitude ($R$-band)')
    
    
    legend_elements = [Line2D([0], [0], marker='o', color='black', ls='--', label='MIST relation',markerfacecolor='black', markersize=10),
                       Line2D([0], [0], marker='^', color='black', ls='-.', label='Carleston $g-i$ (new) relation',markerfacecolor='black', markersize=10),
                       Line2D([0], [0], marker='v', color='black', ls='-.', label='Carleston $g-i$ (old) relation',markerfacecolor='black', markersize=10),
                       Line2D([0], [0], marker='*', color='black', ls='-.', label='Carleston $g-r$ relation',markerfacecolor='black', markersize=10),
                       Line2D([0], [0], marker='x', color='black', ls=':', label='Cantiello $g-i$ relation',markerfacecolor='black', markersize=10)]
    i=0
    for metal in metals:
        legend_elements.append(Patch(facecolor=cmap(norm(i+1)),label=f'metallicity: {round(metal,2)}'))
        i+=1
        
    plt.legend(prop={'size': 15},handles=legend_elements)#,bbox_to_anchor=(1.05, 1),loc='upper left'
    plt.savefig(os.path.join(out_dir,'color_to_absSBFmag.png'))
    #plt.show()
    plt.clf()
    
    median_dists_targetgal = []
    dist_mods_targetgal = []
    
    remove_bias_from_targetgal = False
    i=0
    j=0
    
    for abs_sbf_fit,calibration_id in zip([save_p_mist,save_p_carleston_ginew,save_p_carleston_giold,save_p_carleston_gr,save_p_cantiello_giold],['mist','carleston_ginew','carleston_giold','carleston_gr','cantiello_giold']):
        print('\n',calibration_id,j)
        #### TargetGal SBF
        if remove_bias_from_targetgal:
            dist_mod_targetgal = sbf_measured_targetgal + sbf_bias - abs_sbf_fit(color_targetgal+color_bias)
    
        else: dist_mod_targetgal = sbf_measured_targetgal - abs_sbf_fit(color_targetgal)
        dist_mods_targetgal.append(dist_mod_targetgal[0])
        dist_targetgal = 10**((dist_mod_targetgal + 5)/5) # parsecs
        
        # Calculate maximum error - systematic
        if calibration_id == 'cantiello_giold':
            dist_mod_max_syst = sbf_measured_targetgal - abs_sbf_fit(color_targetgal) + 0.12
            if remove_bias_from_targetgal: dist_mod_max_syst = sbf_measured_targetgal + sbf_bias - abs_sbf_fit(color_targetgal+color_bias) + 0.12
        else:
            #dist_mod_max_syst = sbf_measured_targetgal+std_sbf - abs_sbf_fit(color_targetgal) + 0.26 # Old way, possibly overcounting intrinsic scatter
            dist_mod_max_syst = sbf_measured_targetgal - abs_sbf_fit(color_targetgal) + 0.26
            if remove_bias_from_targetgal: dist_mod_max_syst = sbf_measured_targetgal + sbf_bias - abs_sbf_fit(color_targetgal+color_bias) + 0.26
    
        dist_mod_max_error_syst = abs(dist_mod_max_syst - dist_mod_targetgal)
        
        # Calculate maximum error - statistical

        dist_mod_max_error_stat = dist_mod_stat_uncerts[j]
        
        if include_background_type_difference:
            dist_mod_max_error_stat = dist_mod_max_error_stat + abs(dist_mod_targetgal[0]-SE_dist_mod_measurement[i])
            print('dist_mod_max_error_stat: ',dist_mod_max_error_stat)
            dist_SE = 10**((SE_dist_mod_measurement[i] + 5)/5) # parsecs
            dist_diff = abs(dist_SE-dist_targetgal)
            print('dist difference (kpc):',dist_diff/1e3)
            i+=1
        
        # Convert distance modulus errors to distance errors
        dist_err_stat_lower = [dist_uncerts[j]*1e6+dist_diff]
        dist_err_stat_upper = [dist_uncerts[j]*1e6+dist_diff]
        dist_err_syst_lower = dist_targetgal - 10**((dist_mod_targetgal-dist_mod_max_error_syst + 5)/5) 
        dist_err_syst_upper = 10**((dist_mod_targetgal+dist_mod_max_error_syst + 5)/5) - dist_targetgal
        if j==0: 
            median_dist_errors_targetgal_stat = np.reshape([dist_err_stat_lower[0],dist_err_stat_upper[0]],(2,-1))
            median_dist_errors_targetgal_syst = np.reshape([dist_err_syst_lower[0],dist_err_syst_upper[0]],(2,-1))
        else:
            median_dist_errors_targetgal_stat = np.concatenate((median_dist_errors_targetgal_stat,np.reshape([dist_err_stat_lower[0],dist_err_stat_upper[0]],(2,-1))),axis=1)
            median_dist_errors_targetgal_syst = np.concatenate((median_dist_errors_targetgal_syst,np.reshape([dist_err_syst_lower[0],dist_err_syst_upper[0]],(2,-1))),axis=1)
        median_dists_targetgal.append(dist_targetgal[0])
    
        j+=1

    
    return


#####################################################################################

# Run 'reduce' if called directly
if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Measure distance from SBF magnitude results.')
    parser.add_argument('--in_dir', type=str)
    parser.add_argument('--out_dir', type=str)
    parser.add_argument('--mist_path', type=str)
    parser.add_argument('--feh',type=float)
    args = parser.parse_args()
    
    in_dir = '/fs/scratch/PCON0003/osu10713/Halogas/prelim_testing/polynomial'
    out_dir = '/users/PCON0003/osu10713/Halogas'
    mist_path = '/users/PCON0003/osu10713/ArtPop_tests/MIST'
    get_sbf_distance(in_dir=in_dir, out_dir=out_dir, mist_path=mist_path,model_feh=args.feh)