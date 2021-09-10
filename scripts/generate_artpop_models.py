from astropy.io import fits
import numpy as np
import random
from lbcred.model import artpop_functions

num_models_to_generate = 1 #00
save_dir = '/Users/kirstencasey/'#artpop_models_9000000' #/background_subtraction_tests/single-exposure_input/artpop_models'
in_dir = '/Users/kirstencasey/background_subtraction_tests'

vary_age = 10. #(9.9,10.15)
vary_feh = -1.2 #(-1.2,-1.0)
vary_radius = 300 #(800,900)
vary_n = 0.2 #(0.3,0.7)
vary_theta = 0.  #(0.,180)
vary_ellip = 0. #-0.0926788  #(0.,0.3)
vary_mass = 15000000.0
vary_dist = 3.7 # Mpc
random_state = None
zpt_r = 28. #02708203737313
zpt_b = 28. #26591947732506
include_readnoise = True
include_sky_sb = True


config={'image_dir' : in_dir,
        'color1' : {'psf' : 'lbcr.median_psf_298.fits',
                    'zpt' : zpt_r,
                    'artpop_band' : 'Bessell_R',
                    'extinction' : 0.0,
                    'color_term' : 0,
                    'artpop_model_fn' : 'artpop_model_r.fits',
                    'sky_sb' : 21.0103},
        'color2' : {'psf' : 'lbcb.median_psf_298.fits',
                    'zpt' : zpt_b,
                    'artpop_band' : 'Bessell_B',
                    'extinction' : 0.0,
                    'color_term' : 0,
                    'artpop_model_fn' : 'artpop_model_b.fits',
                    'sky_sb' : 22.7025},
        'include_sky_sb' : include_sky_sb,
        'exposure_time' : 300.281,
        'gain': 1.75,
        'out_dir' : save_dir,
        'add_remnants' : False,
        'pixscale' : 0.224,
        'mist_path' : '/Users/kirstencasey/Documents/Resources/MIST',
        'random_state' : random_state,
        'include_readnoise' : include_readnoise,
        'readnoise' : 12.0,
        'use_src_counts' : False,
        'artpop_model_params' : {'imf' : 'kroupa',
                                'phot_system' : 'UBVRIplus',
                                'xy_dim' : 1051}}

# Use the following for generating random models
'''
params = ['log_age','feh','r_eff','n','theta','ellip','total_mass','distance']
param_values = [vary_age,vary_feh,vary_radius,vary_n,vary_theta,vary_ellip,vary_mass,vary_dist]

for model_num in range(num_models_to_generate):

    print(f'\nGenerating ArtPop model #{model_num}...')

    for param,values in zip(params,param_values):

        if type(values) is tuple: config['artpop_model_params'][param] = random.uniform(values[0],values[1])
        else : config['artpop_model_params'][param] = values

        print(param,':',config['artpop_model_params'][param])

    config['color1']['artpop_model_fn'] = f'artpop_model_{model_num}_r.fits'
    config['color2']['artpop_model_fn'] = f'artpop_model_{model_num}_b.fits'
    model1, model2, src = artpop_functions.run_artimager(config)
'''


# Use the following for generating models in grids:
masses = [6000000] #9000000,6000000,12000000,15000000]
radii = [300]#300,600,900]
ns = [0.5]#[0.2,0.5,0.8]
ells = [0.0]#[0.0, 0.2, 0.4]
ages = [10.1]#10.,9.9,10.1]
fehs = [-0.65]#-1.2,-1.0,-0.65]
theta = 0.
distance = 3.7 # Mpc

config['artpop_model_params']['theta'] = theta
config['artpop_model_params']['distance'] = distance

for mass in masses:
    config['artpop_model_params']['total_mass'] = mass
    for age,feh in zip(ages,fehs):
        config['artpop_model_params']['log_age'] = age
        config['artpop_model_params']['feh'] = feh
        for radius in radii:
            config['artpop_model_params']['r_eff'] = radius
            for n in ns:
                config['artpop_model_params']['n'] = n
                for ell in ells:
                    config['artpop_model_params']['ellip'] = ell
                    for model in range(num_models_to_generate):
                        config['color1']['artpop_model_fn'] = f'artpop_model_mass{mass}_age{age}_feh{feh}_radius{radius}_n{n}_ellip{ell}_iter{model}_r.fits'
                        config['color2']['artpop_model_fn'] = f'artpop_model_mass{mass}_age{age}_feh{feh}_radius{radius}_n{n}_ellip{ell}_iter{model}_b.fits'
                        print(f'Generating ArtPop model {model+1}/{num_models_to_generate}.')
                        print(config['artpop_model_params'])
                        model1, model2, src = artpop_functions.run_artimager(config)
