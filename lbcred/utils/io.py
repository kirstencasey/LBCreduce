"""
Functions for reading input and writing output.
"""
import os
import yaml
import getpass
import numpy as np
from astropy.io import fits
from ..log import logger
from .coordinates import get_today
try:
    import cPickle as pickle
except:
    import pickle
from collections import OrderedDict
from pymfit.model import param_names


__all__ = [
    'default_header_history',
    'mkdir_if_needed',
    'load_config',
    'load_path_or_header',
    'load_path_or_pixels',
    'load_pixels_and_header',
    'load_pickled_data',
    'pickle_data',
    'update_header',
    'temp_fits_file',
    'write_pixels',
]


def default_header_history(prefix=None):
    """
    Return default header line with the date and user name.
    """
    if prefix == None:
        msg = 'Created/updated with DFReduce by {} on {}'
    else:
        msg = prefix + ' by {} on {}'
    return msg.format(getpass.getuser(), get_today())


def load_config(config_fn):
    """
    Load  a yaml configuration file into a dictionary.

    Parameters
    ----------
    config_fn : str
        Configuration file name.

    Returns
    -------
    config : dict
        The configuration.
    """
    with open(config_fn, 'r') as fn:
        logger.debug('Loading configuration file: ' + config_fn)
        config = yaml.load(fn, Loader=yaml.FullLoader)
    return config


def load_path_or_header(path_or_header):
    _type = type(path_or_header)
    if _type == str or _type == np.str_:
        header = fits.getheader(path_or_header)
    elif _type == fits.Header:
        header = path_or_header
    else:
        logger.critical('{} is not a valid path or header type!'.format(_type))
        header = None
    return header


def load_path_or_pixels(path_or_pixels, dtype=None, extname=None):
    """
    Check if the input is a file or numpy array and return a numpy array.

    Parameters
    ----------
    path_or_pixels : str or ndarray or list of one of these
        An image file name or numpy array of its pixels.
    dtype : type (optional)
        If not None, change the image to this data type.

    Returns
    -------
    image_pixels : ndarray
        The image pixels.

    Notes
    -----
    This function is useful for making it possible to pass either an image
    file name or the pixels.
    """
    data = path_or_pixels
    if type(data) == str or type(data) == np.str_:
        data = fits.getdata(data, extname=extname)
    if dtype is not None:
        data = data.astype(dtype)
    return data


def load_pixels_and_header(path_or_pixels, path_or_header=None):
    if not (type(path_or_pixels) == str or type(path_or_pixels) == np.str_):
        msg = 'Must provide header if pixels are given!'
        assert path_or_header is not None, msg
        header = load_path_or_header(path_or_header)
    else:
        if path_or_header is not None:
            header = load_path_or_header(path_or_header)
        else:
            header = load_path_or_header(path_or_pixels)
    pixels = load_path_or_pixels(path_or_pixels)
    return pixels, header


def load_pickled_data(filename):
    """
    Load pickled data

    input
    -----
    filename : string. name of file that
        contains pickled data

    output
    ------
    the unpickled data
    """
    pkl_file = open(filename, 'rb')
    data = pickle.load(pkl_file)
    pkl_file.close()
    return data


def mkdir_if_needed(directory):
    """"
    Create directory if it does not exist.
    """
    if not os.path.isdir(directory):
        logger.info('Creating directory: ' + directory)
        os.makedirs(directory)


def pickle_data(filename, data, warn_overwrite=True):
    """
    Pickle data

    input
    -----
    filename : string. name of file to
        output with full path name.
    data : any data construct. the data
        to be pickled.
    """
    if os.path.isfile(filename) and warn_overwrite:
        logger.warning(filename + ' exists. Will overwrite.')
    pkl_file = open(filename, 'wb')
    logger.debug('Pickling data to ' + filename)
    pickle.dump(data, pkl_file)
    pkl_file.close()


def temp_fits_file(path_or_pixels, tmp_path='/tmp', run_label=None,
                   prefix='tmp',  header=None):
    is_str = type(path_or_pixels) == str or type(path_or_pixels) == np.str_
    if is_str and header is None:
        path = path_or_pixels
        created_tmp = False
    else:
        if is_str:
            path_or_pixels = fits.getdata(path_or_pixels)
        label = '' if run_label is None else '_' + run_label
        fn = '{}{}.fits'.format(prefix, label)
        path = os.path.join(tmp_path, fn)
        logger.debug('Writing temporary fits file {}'.format(path))
        fits.writeto(path, path_or_pixels, header=header, overwrite=True)
        created_tmp = True
    return path, created_tmp


def update_header(path_or_header, msg=None, **kwargs):
    header = load_path_or_header(path_or_header)
    for k, v in kwargs.items():
        header[k] = v
    header.add_history(default_header_history(msg))
    return header


def write_pixels(file_name, pixels, header=None):
    """
    Write pixels to fits file.

    Parameters
    ----------
    file_name : str
        The file name.
    pixels : ndarray
        Image pixels.
    header : astropy.io.fits.Header (optional)
        Image header.

    """
    if os.path.isfile(file_name):
        logger.warning(file_name + ' exists -- will overwrite')
    logger.debug('Writing ' + file_name)
    fits.writeto(file_name, pixels, header=header, overwrite=True)

def rename_fits_file(old_fn, new_fn, delete_old=False, overwrite=False):
    hdu = fits.open(old_fn)
    hdu.writeto(new_fn, overwrite=overwrite)
    hdu.close()
    if delete_old: os.remove(old_fn)
    return

# Adapted from Johnny's pymfitter read_results function
def read_results(filename, models):
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    comments = [l for l in lines if l[0]=='#']
    params = [l for l in lines if l[0] != '#' if l[:2] != '\n'\
                               if l[0] != 'F' if l[:2] != 'X0'\
                               if l[:2] != 'Y0']
    cen_text = [l for l in lines if l[0] != '#'\
                                 if (l[:2] == 'X0' or l[:2] == 'Y0')]
    centers = []

    for i in range(len(cen_text)//2):
        j=i*2
        _, x0, _, _, xerr = cen_text[j].split()
        _, y0, _, _, yerr = cen_text[j+1].split()
        pos_list = [float(x0), float(y0), float(xerr), float(yerr)]
        centers.append(pos_list)

    model_params = OrderedDict()
    for m in models:
        model_params[m] = param_names[m]

    bestfit_params = OrderedDict()
    par_num = 0
    for i in range(len(models)):
        bestfit_params[f'comp_{i+1}'] = {}
        bestfit_params[f'comp_{i+1}']['function'] = models[i]
        bestfit_params[f'comp_{i+1}']['X0'] = centers[i][0]
        bestfit_params[f'comp_{i+1}']['Y0'] = centers[i][1]
        bestfit_params[f'comp_{i+1}']['X0_err'] = centers[i][2]
        bestfit_params[f'comp_{i+1}']['Y0_err'] = centers[i][3]

        for param_name in model_params[models[i]]:
            name, val = params[par_num].split()[:2]
            err = params[par_num].split()[-1]
            assert name == param_name
            bestfit_params[f'comp_{i+1}'].update({param_name: float(val)})
            bestfit_params[f'comp_{i+1}'].update({param_name+'_err': float(err)})
            par_num += 1

    reduced_chisq = [c for c in comments if
                     c.split()[1] == 'Reduced'][0].split()[-1]
    if reduced_chisq != 'none':
        bestfit_params['reduced_chisq'] = float(reduced_chisq)

    return bestfit_params
    
def write_results(results, filename, comments=[]):
    '''
    Only works with Sersic functions right now
    '''
    
    # Open file
    f = open(filename, 'w')
    
    # Write params line by line
    for comment in comments:
        f.write('# '+comment)
        f.write('\n')
    
    f.write('X0\t\t'+str(results['X0'])+' # +/- 9999.9999\n')
    f.write('Y0\t\t'+str(results['Y0'])+' # +/- 9999.9999\n')
    f.write('FUNCTION Sersic\n')
    f.write('PA\t\t' + str(results['PA'])+' # +/- 9999.9999\n')
    f.write('ell\t\t'+str(results['ell'])+' # +/- 9999.9999\n')
    f.write('n\t\t'+str(results['n'])+' # +/- 9999.9999\n')
    f.write('I_e\t\t'+str(results['I_e'])+' # +/- 9999.9999\n')
    f.write('r_e\t\t'+str(results['r_e'])+' # +/- 9999.9999\n')
        
    # Close the file
    f.close()
    
    return
    
