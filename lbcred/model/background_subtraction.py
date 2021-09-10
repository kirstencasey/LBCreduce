from lbcred.detection import sextractor_sky_model
import lbcred, yaml
import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.stats import mad_std, sigma_clip

__all__ = ['background_subtraction']

def generate_mask_for_background_estimation(filename, thresh=20.0):

    mask_fn = filename.replace('backsub','weightimage')
    im = fits.getdata(filename)
    median = np.median(im.flatten())
    std = mad_std(im.flatten())

    mask = im < median + (thresh*std)
    mask &= im > median - (thresh*std)

    mask_hdu = fits.PrimaryHDU(mask.astype(float))
    mask_hdu.writeto(mask_fn, overwrite=True)

    return mask_fn

def plane_skymodel(filename, config, SEsky = None): # Untested

    if SEsky is not None:
        sky = sextractor_sky_model(filename, **SEsky)
    else:
        sky = fits.getdata(filename)

    p_init = models.Planar2D()
    fit_p = fitting.LinearLSQFitter()
    y, x = np.mgrid[:sky.shape[0], :sky.shape[1]]
    print(sky.shape)
    y = np.arange(sky.shape[1])
    x = range(sky.shape[0])
    p = fit_p(p_init,x,y,sky)

    with fits.open(filename, mode='update') as hdul:
        hdul[0].data -= p(x,y)
        hdul[0].header['BACKSUB_TYPE'] = 'Planar2D'
        for i in range(len(p.param_names)):
            hdul[0].header[f'BACKSUB_PARAMS_{i}'] = p.param_names[i]
            hdul[0].header[f'BACKSUB_PARAMVALS_{i}'] = p.parameters[i]
        hdul[0].header['BACKSUB_VAL'] = np.median(p(x,y))
        hdul.close()

    return

def polynomial_skymodel(filename, config, save_models=False, SEsky = None, fn_id=None):

    if SEsky is not None:
        sky = sextractor_sky_model(filename, **SEsky)
    else:
        sky = fits.getdata(filename)

    degree = config['degree']
    p_init = models.Polynomial2D(degree=degree)
    fit_p = fitting.LinearLSQFitter()
    y, x = np.mgrid[:sky.shape[0], :sky.shape[1]]
    p = fit_p(p_init,x,y,sky)

    with fits.open(filename, mode='update') as hdul:
        hdul[0].data -= p(x,y)
        hdul[0].header['BACKSUB_TYPE'] = f'Polynomial2D, degree={degree}'
        for i in range(len(p.param_names)):
            hdul[0].header[f'BACKSUB_PARAMS_{i}'] = p.param_names[i]
            hdul[0].header[f'BACKSUB_PARAMVALS_{i}'] = p.parameters[i]
            hdul[0].header['BACKSUB_VAL'] = np.median(p(x,y))
        hdul.close()

    if save_models:
        model = fits.PrimaryHDU(p(x,y))
        model.header['BACKSUB_TYPE'] = f'Polynomial2D, degree={degree}'
        for i in range(len(p.param_names)):
            model.header[f'BACKSUB_PARAMS_{i}'] = p.param_names[i]
            model.header[f'BACKSUB_PARAMVALS_{i}'] = p.parameters[i]
            model.header['BACKSUB_VAL'] = np.median(p(x,y))
        if fn_id is not None: back_id = f'backmodel_{fn_id}'
        else: back_id = 'backmodel'
        model.writeto(filename.replace('backsub',back_id),overwrite=True)
        if SEsky is not None:
            model = fits.PrimaryHDU(sky)
            model.header['BACKSUB_TYPE'] = f'SE skymodel'
            model.header['BOX_SIZE'] = SEsky['box_size']
            model.writeto(filename.replace('backsub','SEbackmodel'),overwrite=True)

    return


def median_skymodel(filename, config, SEsky=None):

    if SEsky is not None:
        sky = sextractor_sky_model(filename, **SEsky)
    else:
        sky = fits.getdata(filename)

    median = np.median(sky)
    with fits.open(filename, mode='update') as hdul:
        hdul[0].data -= median
        hdul[0].header['BACKSUB_TYPE'] = 'median'
        hdul[0].header['BACKSUB_VAL'] = median
        hdul.close()

    return

def clipped_mean_skymodel(filename, config, SEsky = None): # Untested

    if SEsky is not None:
        sky = source_extractor_skymodel(filename, **SEsky)
    else:
        sky = fits.getdata(filename)
        sky = sigma_clip(sky, sigma=config['thresh'], stdfunc=mad_std)

    mean = np.mean(sky)
    with fits.open(filename, mode='update') as hdul:
        hdul[0].data -= mean
        hdul[0].header['BACKSUB_TYPE'] = 'clipped median'
        hdul[0].header['BACKSUB_THRESH'] = config['thresh']
        hdul[0].header['BACKSUB_VAL'] = mean
        hdul.close()

    return



def background_subtraction(filename, config, fn_id=None):

    SEsky = None

    if config['background_model'] == 'SEsky':

        if config['generate_image_mask']:
            mask_fn = generate_mask_for_background_estimation(filename, thresh=config['image_mask_thresh'])
            config['SEsky']['WEIGHT_IMAGE'] = mask_fn
            config['SEsky']['WEIGHT_TYPE'] = 'MAP_WEIGHT'
        if fn_id is not None: back_id = f'SEbackmodel_{fn_id}'
        else: back_id = 'SEbackmodel'

        sky = sextractor_sky_model(filename, sky_fn=filename.replace('backsub',back_id), **config['SEsky'])

        with fits.open(filename, mode='update') as hdul:
            hdul[0].data -= sky
            hdul[0].header['BACKSUB_TYPE'] = f'SE skymodel'
            hdul[0].header['BOX_SIZE'] = config['SEsky']['box_size']
            hdul[0].header['BACKSUB_VAL'] = np.median(sky)
            hdul.close()

    elif config['background_model'] == 'median':
        if config['median']['initialize_with_SEsky'] :

            if config['generate_image_mask']:
                mask_fn = generate_mask_for_background_estimation(filename, thresh=config['image_mask_thresh'])
                config['SEsky']['WEIGHT_IMAGE'] = mask_fn
                config['SEsky']['WEIGHT_TYPE'] = 'MAP_WEIGHT'

            if fn_id is not None and config['save_background_models']: config['SEsky']['sky_fn'] = filename.replace('backsub',f'SEbackmodel_{fn_id}')
            SEsky = config['SEsky']
        median_skymodel(filename, config['median'], SEsky = SEsky)

    elif config['background_model'] == 'clipped_mean':
        if config['clipped_mean']['initialize_with_SEsky'] :
            if config['generate_image_mask']:
                mask_fn = generate_mask_for_background_estimation(filename, thresh=config['image_mask_thresh'])
                config['SEsky']['WEIGHT_IMAGE'] = mask_fn
                config['SEsky']['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            if fn_id is not None and config['save_background_models']: config['SEsky']['sky_fn'] = filename.replace('backsub',f'SEbackmodel_{fn_id}')
            SEsky = config['SEsky']
        clipped_mean_skymodel(filename, config['clipped_mean'], SEsky = SEsky)

    elif config['background_model'] == 'plane':
        if config['plane']['initialize_with_SEsky'] :
            if config['generate_image_mask']:
                mask_fn = generate_mask_for_background_estimation(filename, thresh=config['image_mask_thresh'])
                config['SEsky']['WEIGHT_IMAGE'] = mask_fn
                config['SEsky']['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            if fn_id is not None and config['save_background_models']: config['SEsky']['sky_fn'] = filename.replace('backsub',f'SEbackmodel_{fn_id}')
            SEsky = config['SEsky']
        plane_skymodel(filename, config['plane'], SEsky = SEsky)

    elif config['background_model'] == 'polynomial':
        if config['polynomial']['initialize_with_SEsky'] :
            if config['generate_image_mask']:
                mask_fn = generate_mask_for_background_estimation(filename, thresh=config['image_mask_thresh'])
                config['SEsky']['WEIGHT_IMAGE'] = mask_fn
                config['SEsky']['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            if fn_id is not None and config['save_background_models']: config['SEsky']['sky_fn'] = filename.replace('backsub',f'SEbackmodel_{fn_id}')
            SEsky = config['SEsky']
        polynomial_skymodel(filename, config['polynomial'], save_models=config['save_background_models'], SEsky = SEsky, fn_id=fn_id)


    return
