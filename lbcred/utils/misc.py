"""
Functions for doing random (but useful!) things.
"""
import logging
import numpy as np
import pandas as pd
from time import time
from functools import wraps
from astropy import units as u
from ..log import logger
from artpop.util import embed_slices
from lbcred.utils import io
from copy import deepcopy
from astropy.nddata import Cutout2D
from astropy.io import fits


__all__ = [
    'check_astropy_units',
    'command_line_override',
    'func_timer',
    'is_list_like',
    'list_of_strings',
    'make_list_like_if_needed',
    'parse_dates_to_list',
    'reverse_dict',
    'setup_logger',
    'slice_image_center',
    'fetch_psf',
    'fetch_cutout',
    'inject_model',
]


def check_astropy_units(value, default_unit):
    t = type(default_unit)
    if type(value) == u.Quantity:
        quantity = value
    elif (t == u.IrreducibleUnit) or (t == u.Unit):
        quantity = value * default_unit
    elif t == str:
        quantity = value * getattr(u, default_unit)
    else:
        raise Exception('default_unit must be an astropy unit or string')
    return quantity


def command_line_override(args, config):
    """
    Overwrite settings in config dictionary that were passed as
    command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments parsed by argparse.
    config : dict
        Pipeline configuration.

    Returns
    -------
    config :  dict
        Updated pipeline configuration.

    Notes
    -----
    The configuration dictionary is modified in place. Set the default
    command-line argument for a parameter to None if you want to use the value
    in the configuration dictionary.
    """
    for k, v in args._get_kwargs():
        if k in config.keys() and v is not None:
            config[k] = v
    return config


def func_timer(f):
    """
    A function decorator to time how long it takes to execute. The time is
    printed as INFO using the logger.
    """
    @wraps(f)
    def wrapper(*args, **kwargs):
        start = time()
        result = f(*args, **kwargs)
        end = time()
        dt = end - start
        if dt > 120:
            dt /= 60
            unit = 'min'
        else:
            unit = 'sec'
        logger.info('{} completed in {:.2f} {}'.format(f.__name__, dt, unit))
        return result
    return wrapper


def is_list_like(check):
    t = type(check)
    c = t == list or t == np.ndarray or t == pd.Series or t == pd.Int64Index
    return c


def list_of_strings(str_or_list):
    """
    Return a list of strings from a single string of comma-separated values.

    Parameters
    ----------
    str_or_list : str or list-like
        Single string of comma-separated values or a list of strings. If it's
        the latter, then the inpits list is simply returned.

    Examples
    --------

            INPUT                                 OUTPUT
    'flag_1,flag_2,flag_3'         --> ['flag_1', 'flag_2', 'flag_3']
    'flag_1, flag_2, flag_3'       --> ['flag_1', 'flag_2', 'flag_3']
    ['flag_1', 'flag_2', 'flag_3'] --> ['flag_1', 'flag_2', 'flag_3']
    """
    if is_list_like(str_or_list):
        ls_str = str_or_list
    elif type(str_or_list) == str:
        ls_str = str_or_list.replace(' ', '').split(',')
    else:
        Exception('{} is not correct type for list of str'.format(str_or_list))
    return ls_str


def make_list_like_if_needed(obj):
    if is_list_like(obj):
        list_like_obj = obj
    else:
        list_like_obj = [obj]
    return list_like_obj


def parse_dates_to_list(dates=None, dates_fn=None):
    if dates_fn is not None:
        dates =  np.loadtxt(dates_fn, dtype=str)
    if not is_list_like(dates):
        dates = [dates]
    return dates


def reverse_dict(d):
    return dict(map(reversed, d.items()))


def setup_logger(level, log_fn=None):
    """
    Setup the pipeline logger.

    Parameters
    ----------
    level : str
        The log level (debug, info, warn, or error).
    log_fn : str (optional)
       Log file name.
    """
    if log_fn is not None:
        fh = logging.FileHandler(log_fn)
        formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s: %(message)s',
            '%Y-%m-%d | %H:%M:%S')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    logger.setLevel(level.upper())


def slice_image_center(shape):
    """
    Generate the 2D slice for the center of an image accounting for whether
    its dimensions are even or odd.
    """
    nrow, ncol = shape

    if nrow % 2 == 0:
        row_c = [nrow//2-1, nrow//2]
    else:
        row_c = nrow // 2

    if ncol % 2 == 0:
        col_c = [ncol//2-1, ncol//2]
    else:
        col_c = ncol // 2

    row_c, col_c = np.meshgrid(row_c, col_c)
    row_c = row_c.flatten()
    col_c = col_c.flatten()

    slice_c = np.s_[row_c, col_c]

    return slice_c


def fetch_psf(ra, dec, layer='dr8', bands='grz', save_files=False, fn_root='legacy_survey', out_dir=''):
    url = 'https://www.legacysurvey.org/viewer/'
    url += f'coadd-psf/?ra={ra}&dec={dec}&layer={layer}&bands={bands}'
    session = requests.Session()
    resp = session.get(url)
    hdulist = fits.open(BytesIO(resp.content))
    if save_file:
        # write to file
        psf_dict = {}
        for i in range(len(bands)):
            band = bands[i]
            filename = os.path.join(out_dir,fn_root+f'_psf_{band}.fits')
            hdulist[i].writeto(filename)
            psf_dict[band] = filename
    else:
        psf_dict = {'grz'[i]: hdulist[i].data for i in range(len(bands))}

    return psf_dict

def fetch_cutout(ra, dec, bands='grz', size = 1500, save_files=False, fn_root='legacy_survey', out_dir='', get_cutout_only = False):

    url_prefix = 'https://www.legacysurvey.org/viewer/'
    url = url_prefix + f'fits-cutout?ra={ra}&dec={dec}&size={size}&'
    url += 'layer=ls-dr8&pixscale=0.262&bands=grz&subimage'
    session = requests.Session()
    resp = session.get(url)
    if get_cutout_only:
        cutout = fits.getdata(BytesIO(resp.content))
        image = {bands[i]: cutout[i, :, :] for i in range(len(bands))}
        return image

    hdulist = fits.open(BytesIO(resp.content))
    if save_files:
        # write to file
        image_dict = {}
        invvar_dict = {}
        for i in range(len(bands)):
            band = bands[i]
            image_fn = os.path.join(out_dir,fn_root+f'_{band}.fits')
            invvar_fn = os.path.join(out_dir,fn_root+f'_invvar_{band}.fits')
            hdulist[[1, 3, 5][i]].writeto(image_fn,overwrite=True)
            hdulist[[2, 4, 6][i]].writeto(invvar_fn,overwrite=True)
            image_dict[band] = image_fn
            invvar_dict[band] = invvar_fn
        return image_dict, invvar_dict
    else:
        image_dict = {}
        invvar_dict = {}
        image_headers = {}
        invvar_headers = {}
        for i in range(len(bands)):
            image_dict[bands[i]] = hdulist[[1, 3, 5][i]].data
            image_headers[bands[i]] = hdulist[[1, 3, 5][i]].header
            invvar_dict[bands[i]] = hdulist[[2, 4, 6][i]].data
            invvar_headers[bands[i]] = hdulist[[2, 4, 6][i]].header
        return image_dict, image_headers, invvar_dict, invvar_headers


def inject_model(image, model, xpos, ypos):

    image = io.load_path_or_pixels(image)
    model = io.load_path_or_pixels(model)

    # work on deep copy in case we want to make adjustments
    mock_image = deepcopy(image)

    # get slices to inject source
    img_slice, arr_slice = embed_slices((xpos,ypos), model.shape, model.shape)

    # inject source into image
    mock_image[img_slice] += model[arr_slice]

    return mock_image

def save_image_odd_shape(img_fn,ext):

    img = fits.open(img_fn)
    shape = min(img[ext].data.shape)

    if shape%2 == 0:
        shape -= 1
        cutout = Cutout2D(img[ext].data, (shape/2,shape/2), shape)
        img[ext].data = cutout.data
        img.writeto(img_fn, overwrite=True)

    return shape


def make_cutout(original_img_fn, position, shape, ext, cutout_fn=None):
    img = fits.open(original_img_fn)
    cutout = Cutout2D(img[ext].data, position, shape)
    img[ext].data = cutout.data
    if cutout_fn is not None:
        img.writeto(cutout_fn, overwrite=True)

    return cutout


def get_chisquare(f_obs, f_exp):
    return np.sum((f_obs-f_exp)**2)
