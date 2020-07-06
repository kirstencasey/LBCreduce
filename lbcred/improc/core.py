"""
Core functions for processing and/or manipulating images.
"""
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from ..utils import io
from ..log import logger
from ..astrometry import pixel_area_map
from .. import ResultStruct


try:
    import fitsio
    from astrometry.util.util import Tan, Sip
    from astrometry.util.resample import resample_with_wcs, OverlapError
    DUSTIN_MAGIC = True
except:
    logger.warning('The magic of Dustin not found -- no resampling')
    DUSTIN_MAGIC = False


__all__ = ['crop_image',
           'replace_dead_pixels', 
           'resample_image']


def crop_image(path_or_pixels, new_shape, save_fn=None, save_header=None):
    image = io.load_path_or_pixels(path_or_pixels)

    o_r, o_c = image.shape # old_row, old_col
    n_r, n_c = new_shape   # new_row, new_col
    odd = (np.array(new_shape) % 2).astype(int)

    s_ = np.s_[
        int(o_r / 2.) - int(n_r / 2.): int(o_r / 2. ) + int(n_r / 2.) + odd[0],
        int(o_c / 2.) - int(n_c / 2.): int(o_c / 2. ) + int(n_c / 2.) + odd[1]
    ]

    cropped_image = image[s_]
    
    if save_fn is not None:
        if save_header is not None:
            msg = 'Image cropped from {} to {}'.format(image.shape, new_shape)
            kw = dict(msg=msg, NAXIS1=n_r, NAXIS2=n_c)
            save_header = io.update_header(save_header, **kw)
        io.write_pixels(save_fn, cropped_image, header=save_header)
    
    return cropped_image


def replace_dead_pixels(image_pixels, padding=1, dead_value=0):
    """ 
    Replace dead pixels in an image with the median of their surrounding
    pixel values. (Watch out for edges!)

    Parameters
    ----------
    image_pixels : ndarray
        Image as a numpy array.
    padding : int
        Sets the size of the region in which to compute the median value of
        "nearby" pixels.

    Returns
    -------
    image_pixels : ndarray
        The image, with dead pixels replaced.
    num_dead_pixels : int
        The number of dead pixels replaced.
    """
    padding = int(padding)

    dead_pixel_indices = np.argwhere(image_pixels==dead_value)
    num_dead_pixels = len(dead_pixel_indices)

    for ind in dead_pixel_indices:
        row, col = ind

        min_row = np.max([0, row - padding])
        max_row = np.min([row + padding + 1, image_pixels.shape[0]])
        min_col = np.max([0, col - padding])
        max_col = np.min([col+padding+1, image_pixels.shape[1]])

        med_value = np.median(image_pixels[min_row:max_row, min_col:max_col])
        image_pixels[row, col] = med_value

    results = ResultStruct(pixels=image_pixels, num_dead_pixels=num_dead_pixels)

    return results


def resample_image(path_or_pixels, ra_c, dec_c, pixscale, width, height,
                   fitsio_header=None, apply_pam=True, normalize_at='peak',
                   interp_type='Lanczos'):
    """
    Resample image using Dustin's magic.

    Thanks Dustin <3

    Parameters
    ----------
    frame_path : str
        Fits file name.
    ra_c : float
        Central RA of resampled image in degrees.
    dec_c : float
        Central DEC of resampled image in degrees.
    pixscale : float
        The desired pixel scale.
    width : int
        The width of the image in pixels.
    height : int
        The height of the image in pixels.
    interp_type : str
        Interpolation method (lanczos or nearest)

    Return
    ------
    results : ResultStruct
        results.resampled_image : ndarray
            The resampled image.
        results.wcs : astropy.wcs.WCS
            Astropy WCS object.
    """
    if not DUSTIN_MAGIC:
        raise Exception('You do not have the magic of Dustin.')

    header = fitsio_header
    if not (type(path_or_pixels) == str or type(path_or_pixels) == np.str_):
        msg = 'Must provide *fitsio* header if pixels are given!'
        assert header is not None, msg
        assert type(header) == fitsio.header.FITSHDR, msg
    elif header is None:
        header = fitsio.read_header(path_or_pixels)
    pixels = io.load_path_or_pixels(path_or_pixels)

    ps = pixscale / 3600.
    x_c = width / 2 + 0.5
    y_c = height / 2 + 0.5
    target_wcs = Tan(ra_c, dec_c, x_c, y_c, -ps, 0., 0.,
                     ps, float(width), float(height))
    resampled = np.zeros((height, width), np.float32)
    _wcs = Sip(header)

    # TODO: find a better way to handle different header types
    astropy_header = fits.Header()
    for k, v in dict(fitsio_header).items():
        astropy_header[k] = v
    astropy_header['EPOCH'] = 2000.0
    pam = pixel_area_map(astropy_header, normalize_at=normalize_at)

    if apply_pam:
        pixels = pixels.astype(np.float64)
        pixels = pixels / pam
    try:
        if interp_type.lower() == 'lanczos':
            Yo, Xo, Yi, Xi, (rim,) = resample_with_wcs(target_wcs,
                                                       _wcs, [pixels])
            resampled[Yo, Xo] += rim
        elif interp_type.lower() == 'nearest':
            Yo, Xo, Yi, Xi, rims = resample_with_wcs(target_wcs, _wcs)
            resampled[Yo, Xo] += pixels[Yi, Xi]
        else:
            raise Exception(interp_type + ' is not a valid interp type')
    except OverlapError:
        logger.critical('WCS frames do not overlap!')
        return False

    # create an astropy wcs object
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = target_wcs.crval
    w.wcs.crpix = target_wcs.crpix
    w.wcs.cdelt = target_wcs.cd[0], target_wcs.cd[-1]
    w.pixel_shape = width, height
    results = ResultStruct(pixels=resampled, wcs=w, pam=pam)

    return results
