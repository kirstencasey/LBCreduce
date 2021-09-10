import artpop, os
import numpy as np
from astropy.io import fits
from astropy import units as u

def run_artimager(config, options={}):

    # Replace config with options
    for key in options:
        config[key] = options[key]

    psf1 = fits.getdata(os.path.join(config['image_dir'],config['color1']['psf']))
    psf2 = fits.getdata(os.path.join(config['image_dir'],config['color2']['psf']))

    if config['include_readnoise']: readnoise = config['readnoise']
    else: readnoise = 0

    zpt_inst = {config['color1']['artpop_band'] : config['color1']['zpt'] + 0.190278,
                config['color2']['artpop_band'] : config['color2']['zpt'] - 0.107512}

    imager = artpop.ArtImager(phot_system=None, zpt_inst=zpt_inst, read_noise = readnoise, random_state = np.random.RandomState(config['random_state']))

    src = artpop.MISTSersicSSP(
        log_age = config['artpop_model_params']['log_age'],             # log of age in years
        feh = config['artpop_model_params']['feh'],                     # metallicity [Fe/H]
        r_eff = config['artpop_model_params']['r_eff'] * u.pc,          # effective radius, 879
        n = config['artpop_model_params']['n'],                         # Sersic index
        theta = config['artpop_model_params']['theta'] * u.deg,         # position angle
        ellip = config['artpop_model_params']['ellip'],                 # ellipticity
        total_mass = config['artpop_model_params']['total_mass'],       # stellar mass
        imf = config['artpop_model_params']['imf'],                     # initial mass function
        phot_system = config['artpop_model_params']['phot_system'],     # photometric system
        distance = config['artpop_model_params']['distance'] * u.Mpc,   # distance to system
        xy_dim = config['artpop_model_params']['xy_dim'],               # image dimension
        random_state = np.random.RandomState(config['random_state']),   # random state for reproducibility
        mist_path = config['mist_path'],
        pixel_scale = config['pixscale'],                               # pixel scale in arcsec / pixel
        add_remnants = config['add_remnants']
    )

    #colors = src.mags[config['color2']['artpop_band']] - src.mags[config['color1']['artpop_band']]
    color_inst = (src.sp.total_mag(config['color2']['artpop_band'])-src.sp.total_mag(config['color1']['artpop_band'])+config['color2']['extinction']-config['color1']['extinction'])/(config['color1']['color_term']-config['color2']['color_term']+1)
    #color_inst = (src.sp.total_mag(config['color2']['artpop_band'])-0.107512-src.sp.total_mag(config['color1']['artpop_band'])-0.190278+config['color2']['extinction']-config['color1']['extinction']-config['color2']['zpt']+config['color1']['zpt'])/(config['color1']['color_term']-config['color2']['color_term']+1)
    term1 = color_inst*config['color1']['color_term']
    term2 = color_inst*config['color2']['color_term']

    src.mags[config['color1']['artpop_band']] = src.mags[config['color1']['artpop_band']] + 0.190278 + config['color1']['extinction'] + term1    # Vega -> AB mags, color/extinction correction
    src.mags[config['color2']['artpop_band']] = src.mags[config['color2']['artpop_band']] - 0.107512 + config['color2']['extinction'] + term2   # Vega -> AB mags, color/extinction correction

    if config['include_sky_sb']:
        sky_sb1 = config['color1']['sky_sb']
        sky_sb2 = config['color2']['sky_sb']
    else:
        sky_sb1 = None
        sky_sb2 = None

    obs1 = imager.observe(src, bandpass=config['color1']['artpop_band'], psf=psf1, exptime=config['exposure_time'], sky_sb=sky_sb1)
    obs2 = imager.observe(src, bandpass=config['color2']['artpop_band'], psf=psf2, exptime=config['exposure_time'], sky_sb=sky_sb2)
    model1 = obs1.raw_counts
    model2 = obs2.raw_counts

    if config['include_sky_sb']:
        print('SKY COUNTS: ',obs1.sky_counts)
        print('SKY COUNTS: ',obs2.sky_counts)
        model1 -= obs1.sky_counts
        model2 -= obs2.sky_counts

    model1 = model1 / config['gain']
    model2 = model2 / config['gain']

    header = fits.Header(config['artpop_model_params'])
    header[config['color1']['artpop_band']+'_mag'] = src.sp.total_mag(config['color1']['artpop_band'])
    header[config['color2']['artpop_band']+'_mag'] = src.sp.total_mag(config['color2']['artpop_band'])
    header[config['color1']['artpop_band']+'_sbfmag'] = src.sp.sbf_mag(config['color1']['artpop_band'])
    header[config['color1']['artpop_band']+'_zpt'] = config['color1']['zpt']
    header[config['color2']['artpop_band']+'_zpt'] = config['color2']['zpt']
    header[config['color1']['artpop_band']+'_extinction'] = config['color1']['extinction']
    header[config['color2']['artpop_band']+'_extinction'] = config['color2']['extinction']
    header[config['color1']['artpop_band']+'_color_term'] = config['color1']['color_term']
    header[config['color2']['artpop_band']+'_color_term'] = config['color2']['color_term']
    header['sky_sb_included'] = config['include_sky_sb']
    header['readnoise_included'] = config['include_readnoise']
    if config['include_sky_sb']:
        header[config['color1']['artpop_band']+'_sky_sb'] = config['color1']['sky_sb']
        header[config['color2']['artpop_band']+'_sky_sb'] = config['color2']['sky_sb']
    if config['include_readnoise']:
        header['readnoise'] = config['readnoise']
    header['color'] = src.sp.total_mag(config['color2']['artpop_band']) - src.sp.total_mag(config['color1']['artpop_band'])

    hdu0 = fits.PrimaryHDU(header=header)
    header['EXTNAME'] = 'raw_counts'
    hdu10 = fits.ImageHDU(model1, header=header)
    hdu20 = fits.ImageHDU(model2, header=header)

    header['EXTNAME'] = 'src_counts'
    header['NOTE'] = 'Poisson noise added to original src_counts.'

    if config['include_sky_sb'] or config['include_readnoise'] :
        rng = artpop.util.check_random_state(config['random_state'])
        obs1.src_counts[obs1.src_counts < 0] = 0.
        obs2.src_counts[obs2.src_counts < 0] = 0.
        model1_src = rng.poisson(obs1.src_counts) / config['gain']
        model2_src = rng.poisson(obs2.src_counts) / config['gain']
    else :
        model1_src = model1
        model2_src = model2

    hdu11 = fits.ImageHDU(model1_src, header=header)
    new_hdul1 = fits.HDUList([hdu0, hdu10, hdu11])
    new_hdul1.writeto(os.path.join(config['out_dir'], config['color1']['artpop_model_fn']), overwrite=True)
    hdu21 = fits.ImageHDU(model2_src, header=header)
    new_hdul2 = fits.HDUList([hdu0, hdu20, hdu21])
    new_hdul2.writeto(os.path.join(config['out_dir'], config['color2']['artpop_model_fn']), overwrite=True)

    if config['use_src_counts'] : return model1_src, model2_src, src

    else : return model1, model2, src

def run_idealimager(config, image_fns, psf_fns):


    return
