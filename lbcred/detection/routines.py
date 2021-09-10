"""
Pre-cooked source detection routines.
"""
from . import sextractor
import os
from astropy.io import fits
from astropy.table import vstack, Table
import numpy as np


__all__ = ['extract_bright_stars','sextractor_sky_model']


def extract_bright_stars(path_or_pixels,catalog_path=None):
    cat = sextractor.run(path_or_pixels,
                         DETECT_MINAREA=3,
                         DETECT_THRESH=10,
                         PIXEL_SCALE=0.224,
                         catalog_path=catalog_path)
    star_query = 'FLAGS==0 and ISOAREA_IMAGE > 5 and \
                  FWHM_IMAGE > 1 and FWHM_IMAGE < 26'
    cat = cat[cat.to_pandas().query(star_query).index.values]
    return cat


def sextractor_sky_model(path_or_pixels, run_label=None, tmp_path='/tmp',
                         sky_fn=None, box_size=512, **sextractor_options):
    label = '' if run_label is None else '_' + run_label
    if sky_fn is not None:
        created_tmp = False
    else:
        sky_fn = os.path.join(tmp_path, f'skymodel{label}.fits')
        created_tmp = True
    options = {}
    for k, v in sextractor_options.items():
        options[k.upper()] = v
    config = dict(CHECKIMAGE_TYPE='BACKGROUND',
                  CHECKIMAGE_NAME=sky_fn,
                  BACK_SIZE=box_size,
                  tmp_path=tmp_path,
                  run_label=run_label)
    config['DETECT_THRESH'] = options.pop('DETECT_THRESH', 2)
    config.update(options)
    sextractor.run(path_or_pixels, **config)
    sky = fits.getdata(sky_fn)
    '''
    if created_tmp:
        os.remove(sky_fn)
    '''
    return sky

def create_mask_from_cat(cat, mask_shape, mask_radius_factor=1.0, min_radius=1., max_radius=20., mask_fn=None):

    # Create blank mask
    mask = np.zeros(mask_shape)

    if len(cat) == 0: return mask.astype(bool)

    # Filter catalog based on input radius features
    filtered_cat = None
    for obj in cat:

        if obj['FLUX_RADIUS'] > min_radius and obj['FLUX_RADIUS'] < max_radius:
            if filtered_cat is None: filtered_cat = Table(obj)
            else: filtered_cat = vstack([filtered_cat, Table(obj)])

    if len(filtered_cat) == 0: return mask.astype(bool)

    filtered_cat['FLUX_RADIUS'] *= mask_radius_factor

    # Add circles to mask
    ii, jj = np.mgrid[:mask.shape[0], :mask.shape[1]]

    for obj in filtered_cat:

        r_circ = obj['FLUX_RADIUS']
        i_c, j_c = (obj['Y_IMAGE'], obj['X_IMAGE'])
        circ_mask = ((ii - i_c)**2 + (jj - j_c)**2) - r_circ**2 < 0

        mask = mask.astype(bool) | circ_mask

    if mask_fn is not None:
        mask = mask.astype(float)
        mask_hdu = fits.PrimaryHDU(mask)
        mask_hdu.writeto(mask_fn,overwrite=True)

    return mask.astype(bool)
