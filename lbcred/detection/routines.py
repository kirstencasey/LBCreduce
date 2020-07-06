"""
Pre-cooked source detection routines. 
"""
from . import sextractor


__all__ = ['extract_bright_stars']


def extract_bright_stars(path_or_pixels):
    cat = sextractor.run(path_or_pixels,
                         DETECT_MINAREA=3,
                         DETECT_THRESH=10, 
                         PIXEL_SCALE=0.225)
    star_query = 'FLAGS==0 and ISOAREA_IMAGE > 5 and \
                  FWHM_IMAGE > 1 and FWHM_IMAGE < 26'
    cat = cat[cat.to_pandas().query(star_query).index.values]
    return cat
