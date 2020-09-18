"""
Pre-cooked source detection routines.
"""
from . import sextractor
import os
from astropy.io import fits


__all__ = ['extract_bright_stars','sextractor_sky_model']


def extract_bright_stars(path_or_pixels,catalog_path=None):
    cat = sextractor.run(path_or_pixels,
                         DETECT_MINAREA=3,
                         DETECT_THRESH=10,
                         PIXEL_SCALE=0.225,
                         catalog_path=catalog_path)
    star_query = 'FLAGS==0 and ISOAREA_IMAGE > 5 and \
                  FWHM_IMAGE > 1 and FWHM_IMAGE < 26'
    cat = cat[cat.to_pandas().query(star_query).index.values]
    return cat


def sextractor_sky_model(path_or_pixels, run_label=None, tmp_path='/tmp',
                         sky_fn=None, **sextractor_options):
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
                  tmp_path=tmp_path,
                  run_label=run_label)
    config['DETECT_THRESH'] = options.pop('DETECT_THRESH', 2)
    config.update(options)
    sextractor.run(path_or_pixels, **config)
    sky = fits.getdata(sky_fn)
    if created_tmp:
        os.remove(sky_fn)
    return sky
