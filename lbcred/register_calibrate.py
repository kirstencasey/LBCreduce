import os, yaml, lbcred, random, ccdproc, warnings, sys, glob
import numpy as np
from shutil import copyfile
from astropy.io import fits, ascii
from astropy.table import  Column, Table, vstack
from astropy.stats import sigma_clip
from astropy import log
from astropy.modeling import models, fitting
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import Longitude, Latitude, AltAz, SkyCoord, EarthLocation
from astropy.wcs import FITSFixedWarning, WCS
from astropy.io.fits.verify import VerifyWarning
from astropy.nddata import CCDData, Cutout2D
from astropy.stats import mad_std

from lbcred.detection import extract_bright_stars
from lbcred.astrometry import solve_field, check_astrometric_solution, TweakWCS
from lbcred import improc, utils, logger, image, interactive, tools, reduce, detection
from lbcred.utils import misc, coordinates, io
from lbcred.model.background_subtraction import background_subtraction
from lbcred.model import artpop_functions, imfit

from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import chisquare
#import astromatic_wrapper as aw
from matplotlib.lines import Line2D
from photutils import CircularAperture
from photutils import aperture_photometry
from matplotlib import colors

lbt = EarthLocation.of_site('Large Binocular Telescope')
warnings.simplefilter('ignore', category=FITSFixedWarning)
warnings.simplefilter('ignore', category=VerifyWarning)


def fetch_files(path, band='R', exptime_range=[0, 50],
                name_must_contain='M81blob.proc'):
    files = glob.glob(os.path.join(path, f'lbc{band.lower()}*{name_must_contain}*'))

    files = [f for f in files if\
             float(fits.getheader(f)['EXPTIME']) > exptime_range[0] if\
             float(fits.getheader(f)['EXPTIME']) < exptime_range[1]]
    # sort files by time
    # NOTE: This assumes the file
    # names contain the exposure number
    files.sort()
    return files


def get_date_time(cat):
    date = cat['filename'][0].split('.')[1]
    time = int(cat['filename'][0].split('.')[2].split('-')[0])

    return date, time

def find_closest_cat(given_cat,cats_to_search):

    '''
    This function finds the catalog that corresponds to the closest matching time stamp of the given catalog.
    Ex. If you give it a catalog as given_cat, it finds the closest matching time stamp of all the catalogs in
    cats_to_search.
    Assumes there is at least one catalog in cats_to_search that has the same date as the given_cat
    '''

    given_date, given_time = get_date_time(given_cat)

    # Loop through cats_to_search
    closest_diff = 10000000
    for cat in cats_to_search:
        date, time = get_date_time(cat)

        if date == given_date:
            diff = abs(given_time - time)
            if diff < closest_diff:
                closest_diff = diff
                closest_cat = cat

    return closest_cat

def calculate_air_mass(skycoord, utc_time, location):
    """
    Assumes approximation from Hardie (1962).
    """
    utc_time = Time(utc_time, format='isot')
    nms_altaz = AltAz(obstime=utc_time, location=location)
    target_nms = skycoord.transform_to(nms_altaz)
    secz = target_nms.secz.value
    air_mass = secz -  0.0018167 * (secz - 1) -\
                       0.0028750 * (secz - 1)**2 -\
                       0.0008083 * (secz - 1)**3
    return air_mass


def find_largest_airmass_diff(cat,airmass_name):
    end = len(cat)
    diff = 0
    ind1 = 0
    ind2 = 1
    for i in range(len(cat)):
        j=i+1
        while j < end:
            if abs(cat[airmass_name][i]-cat[airmass_name][j]) > diff:
                diff = abs(cat[airmass_name][i]-cat[airmass_name][j])
                ind1 = i
                ind2 = j
            j+=1
    return cat[np.array([ind1,ind2])]

def make_region_files(cat, ftype, out_dir, coordtype='physical', plot_panstarrs = False, panstarrs_cat = None):

    cat['id'] = np.zeros(len(cat),dtype=int)
    if plot_panstarrs: panstarrs_cat['id'] = np.zeros(len(cat),dtype=int)

    files = np.array([])
    filenames = np.array([])

    if ftype == 'cali': color = 'blue'
    else: color = 'purple'

    num = 0
    for star in cat:
        cat[num]['id'] = num

        if plot_panstarrs:
            panstarrs_cat[num]['id'] = num
            x_p = panstarrs_cat[num]['ra']
            y_p = panstarrs_cat[num]['dec']
            line_p = f'J2000; circle {x_p}d {y_p}d 3p # color=red\n'
        if coordtype == 'J2000':
            x = star['ALPHA_J2000']
            y = star['DELTA_J2000']
            line = f'J2000; circle {x}d {y}d 3p # text={{{num}}}\n'
        elif coordtype == 'physical':
            x = star['X_IMAGE']
            y = star['Y_IMAGE']
            line = f'circle {x} {y} 3p # text={{{num}}}\n'

        if star['filename'] in filenames:
            # Add to existing region file
            file = files[np.where(filenames==star['filename'])][0]
            file.write(line)
            if plot_panstarrs:
                file.write(line_p)
        else:
            # Make new region file
            fname = star['filename']
            filenames = np.append(filenames,fname)

            file = open(os.path.join(out_dir,fname.replace('.fits','.reg')),'w')
            file.write(f'global color={color}\n')
            file.write(line)
            if plot_panstarrs:
                file.write(line_p)
            files = np.append(files,file)
        num+=1

    # Save all files
    for file in files:
        file.close()

    return cat, panstarrs_cat

def apply_aper_color_corrections(file_dir,cat_r,cat_b,cat_p,mag_name_in,mag_name_out,aper_radius,aper_corrections_r,aper_corrections_b,aper_correction_files_r,aper_correction_files_b,gain,texp,saturation_limit,model='linear',degree=None,use_SE_fluxes=False):

    size = (51,51)
    aper_sums_r = []
    aper_sums_b = []
    aper_mags_r = []
    aper_mags_b = []
    corrected_fluxes_r = []
    corrected_fluxes_b = []
    mask_r = []
    mask_b = []
    zp = Table()

    for idx in range(len(cat_p)):
        file_r = cat_r[idx]['filename']
        file_b = cat_b[idx]['filename']
        if len(aper_corrections_r) == 1:
            exp_median_r = aper_corrections_r[0]
        if len(aper_corrections_b) == 1:
            exp_median_b = aper_corrections_r[0]
        else:
            exp_median_r = aper_corrections_r[aper_correction_files_r.index(file_r)]
            exp_median_b = aper_corrections_b[aper_correction_files_b.index(file_b)]

        if not use_SE_fluxes:
            im_hdul_r = fits.open(os.path.join(file_dir,file_r),ignore_blank=True)
            position_r = (cat_r[idx]['X_IMAGE']-1,cat_r[idx]['Y_IMAGE']-1)
            cutout_r = Cutout2D(im_hdul_r[0].data, position_r, size)
            im_hdul_b = fits.open(os.path.join(file_dir,file_b),ignore_blank=True)
            position_b = (cat_b[idx]['X_IMAGE']-1,cat_b[idx]['Y_IMAGE']-1)
            cutout_b = Cutout2D(im_hdul_b[0].data, position_b, size)

            if np.sum(cutout_r.data >= saturation_limit) != 0:
                mask_r.append(False)
            else: mask_r.append(True)
            if np.sum(cutout_b.data >= saturation_limit) != 0:
                mask_b.append(False)
            else: mask_b.append(True)

            x_corrected_r = position_r[0]-cutout_r.origin_original[0]
            y_corrected_r = position_r[1]-cutout_r.origin_original[1]
            position_corrected_r = (x_corrected_r,y_corrected_r)
            aperture_r = CircularAperture(position_corrected_r, r=aper_radius)
            phot_table_r = aperture_photometry(cutout_r.data, aperture_r)
            for phot_col in phot_table_r.colnames:
                 phot_table_r[phot_col].info.format = '%.8g'  # for consistent table output

            aper_sums_r.append(phot_table_r[0][f'aperture_sum'])

            x_corrected_b = position_b[0]-cutout_b.origin_original[0]
            y_corrected_b = position_b[1]-cutout_b.origin_original[1]
            position_corrected_b = (x_corrected_b,y_corrected_b)
            aperture_b = CircularAperture(position_corrected_b, r=aper_radius)
            phot_table_b = aperture_photometry(cutout_b.data, aperture_b)
            for phot_col in phot_table_b.colnames:
                 phot_table_b[phot_col].info.format = '%.8g'  # for consistent table output

            aper_sums_b.append(phot_table_b[0][f'aperture_sum'])
            im_hdul_r.close()
            im_hdul_b.close()

            corrected_flux_r = phot_table_r[0][f'aperture_sum'] / (1-exp_median_r)
            corrected_fluxes_r.append(corrected_flux_r)
            corrected_flux_b = phot_table_b[0][f'aperture_sum'] / (1-exp_median_b)
            corrected_fluxes_b.append(corrected_flux_b)

        else:
            corrected_flux_r = cat_r[idx][f'FLUX_APER_{aper_radius-1}']
            corrected_flux_b = cat_b[idx][f'FLUX_APER_{aper_radius-1}']
            corrected_fluxes_r.append(corrected_flux_r)
            corrected_fluxes_b.append(corrected_flux_b)


        # Calculate magnitude, get panstarrs magnitude
        aper_mags_r.append(-2.5 * np.log10(corrected_flux_r*gain/texp))
        aper_mags_b.append(-2.5 * np.log10(corrected_flux_b*gain/texp))

    cat_r['aper_corrected_flux'] = corrected_fluxes_r
    cat_b['aper_corrected_flux'] = corrected_fluxes_b
    cat_r['aper_corrected_mag'] = aper_mags_r
    cat_b['aper_corrected_mag'] = aper_mags_b

    # Remove saturated stars
    if not use_SE_fluxes:
        full_mask = ~(~np.asarray(mask_r)+(~np.asarray(mask_b)))
        cat_r = cat_r[full_mask]
        cat_b = cat_b[full_mask]
        cat_p = cat_p[full_mask]

    fit = fitting.LinearLSQFitter()
    if model == 'linear':
        line_init = models.Linear1D()
    elif model == 'polynomial':
        line_init = models.Polynomial1D(degree=degree,domain=(0,3))

    fig,ax = plt.subplots(2,2,figsize=(15,15))
    files_r = np.asarray(np.unique(cat_r['filename']))
    files_b = []
    zp_arr_r = []
    zp_arr_b = []
    for im_r in files_r:

        #label = im_r.split('.')[2].split('-')[0]
        mask = cat_r['filename']==im_r
        im_cat_r = cat_r[mask]
        im_cat_b = cat_b[mask]
        im_cat_p = cat_p[mask]
        files_b.append(im_cat_b['filename'][0])

        x_p = im_cat_p['B-BESSEL_linear'] - im_cat_p['R-BESSEL_linear']
        x_inst = im_cat_b['aper_corrected_mag'] - im_cat_r['aper_corrected_mag']
        y_r = im_cat_r['aper_corrected_mag'] - im_cat_p['R-BESSEL_linear']
        y_b = im_cat_b['aper_corrected_mag'] - im_cat_p['B-BESSEL_linear']

        fitted_line_inst_b = fit(line_init, x_inst, y_b)
        fitted_line_inst_r = fit(line_init, x_inst, y_r)

        zp_arr_r.append(fitted_line_inst_r(0.))
        zp_arr_b.append(fitted_line_inst_b(0.))

        ax[1][0].scatter(x_inst,y_b,color='darkblue')
        ax[1][0].plot(x_inst, fitted_line_inst_b(x_inst))

        ax[1][1].scatter(x_inst,y_r,color='darkred')
        ax[1][1].plot(x_inst, fitted_line_inst_r(x_inst))


    x_p = cat_p['B-BESSEL_linear'] - cat_p['R-BESSEL_linear']
    x_inst = cat_b['aper_corrected_mag'] - cat_r['aper_corrected_mag']
    y_r = cat_r['aper_corrected_mag'] - cat_p['R-BESSEL_linear']
    y_b = cat_b['aper_corrected_mag'] - cat_p['B-BESSEL_linear']

    fitted_line_pan_b = fit(line_init, x_p, y_b)
    fitted_line_pan_r = fit(line_init, x_p, y_r)
    fitted_line_inst_b = fit(line_init, x_inst, y_b)
    fitted_line_inst_r = fit(line_init, x_inst, y_r)
    zp_arr_r.append(fitted_line_inst_r(0.))
    zp_arr_b.append(fitted_line_inst_b(0.))
    files_r = np.append(files_r,'all_exposures')
    files_b = np.append(files_b,'all_exposures')
    zp_r = Table()
    zp_r['filename'] = files_r
    zp_r['zp'] = zp_arr_r
    zp_b = Table()
    zp_b['filename'] = files_b
    zp_b['zp'] = zp_arr_b

    ax[0][0].scatter(x_p,y_b,color='darkblue')
    ax[0][0].set_xlabel(f'Pan-STARRS B-R')
    ax[0][0].set_ylabel(f'Instrumental B - Pan-STARRS B')
    ax[0][0].plot(x_p, fitted_line_pan_b(x_p), 'k-', label='Fitted Model')
    #ax[0][0].text(1.0,-27.7,f'slope: {round(fitted_line_inst_all_b.slope.value, 2)}\nintercept: {round(fitted_line_inst_all_b.intercept.value, 2)}',fontsize=14)
    ax[0][1].scatter(x_p,y_r,color='darkred')
    ax[0][1].set_xlabel(f'Pan-STARRS B-R')
    ax[0][1].set_ylabel(f'Instrumental R - Pan-STARRS R')
    ax[0][1].plot(x_p, fitted_line_pan_r(x_p), 'k-', label='Fitted Model')

    #ax[1][0].scatter(x_inst,y_b,color='darkblue')
    ax[1][0].set_xlabel(f'Instrumental B-R')
    ax[1][0].set_ylabel(f'Instrumental B - Pan-STARRS B')
    ax[1][0].plot(x_inst, fitted_line_inst_b(x_inst), 'k-', label='Fitted Model')

    #ax[1][1].scatter(x_inst,y_r,color='darkred')
    ax[1][1].set_xlabel(f'Instrumental B-R')
    ax[1][1].set_ylabel(f'Instrumental R - Pan-STARRS R')
    ax[1][1].plot(x_inst, fitted_line_inst_r(x_inst), 'k-', label='Fitted Model')
    plt.savefig(os.path.join(file_dir,'diagnostic-plots','compare_instrumental_and_catalog_mags.png'))

    mag_color_corrected_pan_r = cat_r['aper_corrected_mag'] - fitted_line_pan_r(x_p)
    mag_color_corrected_pan_b = cat_b['aper_corrected_mag'] - fitted_line_pan_b(x_p)
    mag_color_corrected_inst_r = cat_r['aper_corrected_mag'] - fitted_line_inst_r(x_inst)
    mag_color_corrected_inst_b = cat_b['aper_corrected_mag'] - fitted_line_inst_b(x_inst)

    cat_r[mag_name_out+'_pan'] = mag_color_corrected_pan_r
    cat_r[mag_name_out+'_inst'] = mag_color_corrected_inst_r
    cat_b[mag_name_out+'_pan'] = mag_color_corrected_pan_b
    cat_b[mag_name_out+'_inst'] = mag_color_corrected_inst_b

    return cat_r, cat_b, cat_p, zp_r, zp_b, fitted_line_inst_r.slope.value, fitted_line_inst_b.slope.value

def plot_residuals(cat_r,cat_b,cat_p,mag_name,save_figs=False,fig_name=None):#,separate_exposures=False,files_r=None,files_b=None):
    '''
    if separate_exposures:
        num_plots = len(files_r)+1
        fig,ax = plt.subplots(num_plots,2,figsize=(15,num_plots*8))
    '''

    fig,ax = plt.subplots(2,2,figsize=(15,15))
    x = cat_p['B-BESSEL_linear']
    y = cat_b[mag_name+'_pan'] - cat_p['B-BESSEL_linear']
    y_rms = np.sqrt(np.mean(np.asarray(y)**2))
    ylim = max([abs(min(y)),abs(max(y))])
    ax[0][0].scatter(x,y,color='darkblue')
    ax[0][0].plot(np.arange(30),np.zeros(30),linestyle='dashed', color='black')
    ax[0][0].set_xlim(min(x)-0.5,max(x)+0.5)
    ax[0][0].set_ylim(-ylim-0.1,ylim+0.1)
    ax[0][0].text(19,-0.3,f'rms: {round(y_rms,4)}',fontsize=20)
    ax[0][0].set_xlabel('Pan-STARRS B (mag)')
    ax[0][0].set_ylabel('Calibrated B mag - Pan-STARRS B mag')
    ax[0][0].title.set_text('Residuals Using Pan-STARRS Colors for Color Terms')

    x = cat_p['R-BESSEL_linear']
    y = cat_r[mag_name+'_pan'] - cat_p['R-BESSEL_linear']
    y_rms = np.sqrt(np.mean(np.asarray(y)**2))
    ax[0][1].scatter(x,y,color='darkred')
    ax[0][1].plot(np.arange(30),np.zeros(30),linestyle='dashed', color='black')
    ax[0][1].set_xlim(min(x)-0.5,max(x)+0.5)
    ax[0][1].set_ylim(-ylim-0.1,ylim+0.1)
    ax[0][1].text(17.5,-0.3,f'rms: {round(y_rms,4)}',fontsize=20)
    ax[0][1].set_xlabel('Pan-STARRS R (mag)')
    ax[0][1].set_ylabel('Calibrated R mag - Pan-STARRS R mag')
    ax[0][1].title.set_text('Residuals Using Pan-STARRS Colors for Color Terms')

    x = cat_p['B-BESSEL_linear']
    y = cat_b[mag_name+'_inst'] - cat_p['B-BESSEL_linear']
    y_rms = np.sqrt(np.mean(np.asarray(y)**2))
    ax[1][0].scatter(x,y,color='darkblue')
    ax[1][0].plot(np.arange(30),np.zeros(30),linestyle='dashed', color='black')
    ax[1][0].set_xlim(min(x)-0.5,max(x)+0.5)
    ax[1][0].set_ylim(-ylim-0.1,ylim+0.1)
    ax[1][0].text(19,-0.3,f'rms: {round(y_rms,4)}',fontsize=20)
    ax[1][0].set_xlabel('Pan-STARRS B (mag)')
    ax[1][0].set_ylabel('Calibrated B mag - Pan-STARRS B mag')
    ax[1][0].title.set_text('Residuals Using Calibrated Inst Colors for Color Terms')

    x = cat_p['R-BESSEL_linear']
    y = cat_r[mag_name+'_inst'] - cat_p['R-BESSEL_linear']
    y_rms = np.sqrt(np.mean(np.asarray(y)**2))
    ax[1][1].scatter(x,y,color='darkred')
    ax[1][1].plot(np.arange(30),np.zeros(30),linestyle='dashed', color='black')
    ax[1][1].set_xlim(min(x)-0.5,max(x)+0.5)
    ax[1][1].set_ylim(-ylim-0.1,ylim+0.1)
    ax[1][1].text(17.5,-0.3,f'rms: {round(y_rms,4)}',fontsize=20)
    ax[1][1].set_xlabel('Pan-STARRS R (mag)')
    ax[1][1].set_ylabel('Calibrated R mag - Pan-STARRS R mag')
    ax[1][1].title.set_text('Residuals Using Calibrated Inst Colors for Color Terms')

    if save_figs:
        plt.savefig(fig_name)
    return

def replace_dead_pixels(image_pixels, padding=5, dead_value=0, mask=None):
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
    dead_value : int
        Value of the pixels considered "dead"
    mask : ndarray
        Array the same shape as image_pixels where pixels to be masked = True
        and unmasked pixels = False
    Returns
    -------
    image_pixels : ndarray
        The image, with dead pixels replaced.
    num_dead_pixels : int
        The number of dead pixels replaced.
    """
    padding = int(padding)

    image_pixels[mask] = dead_value+1
    dead_pixel_indices = np.argwhere(image_pixels==dead_value)
    image_pixels[mask] = dead_value#-1
    num_dead_pixels = len(dead_pixel_indices)
    for ind in dead_pixel_indices:
        row, col = ind
        min_row = np.max([0, row - padding])
        max_row = np.min([row + padding + 1, image_pixels.shape[0]])
        min_col = np.max([0, col - padding])
        max_col = np.min([col+padding+1, image_pixels.shape[1]])
        med_value = np.median(image_pixels[min_row:max_row, min_col:max_col])
        image_pixels[row, col] = med_value
    return image_pixels, num_dead_pixels


@utils.func_timer
def register_images(tmp_path, bandpass, index_path=None, ref_cat=None, make_plots=True, config_fn=None, back_fn_id=None):

    if type(config_fn) is not dict:
        with open(config_fn, 'r') as config:
            config = yaml.load(config, Loader=yaml.FullLoader)
        # Save copy of config to output directory
        copyfile(config_fn,os.path.join(config['out_dir'],config_fn.split('/')[-1]))
    else:
        config = config_fn
        # Save copy of config to output directory
        with open(os.path.join(config['out_dir'],'register_calibration-config.yml'), 'w') as outfile:
            yaml.dump(config, outfile, default_flow_style=False)

    glob_select = config['glob_select']
    data_path = config['image_dir']
    out_path = config['out_dir']
    center = [config['center_ra'],config['center_dec']]

    # fetch reference catalog if necessary
    if ref_cat is None:
        ref_cat = config['reference_catalog']
    if ref_cat is None:
        fn =  f'panstarrs-{center[0]:.1f}-{center[1]:.1f}.dat'
        fn = os.path.join(out_path, fn)
        print('PANSTARRS FN: ',fn)
        # if has been fetched previously, load local copy
        if os.path.isfile(fn):
            logger.info('reading reference catalog from ' + fn)
            ref_cat = Table.read(fn)
        else:
            # otherwise fetch using pyvo
            from lbcred.tap import PanstarrsTAP
            panstarrs_tap = PanstarrsTAP()
            logger.warning('no reference catalog. fetching PanSTARRS catalog')
            query_template = panstarrs_tap.\
                             default_query_template.\
                             replace('10', '16').\
                             replace('23', '22')
            ref_cat = panstarrs_tap.quey_region(ra_c=config['center_ra'], dec_c=config['center_dec'], radius=0.25,
                                                query_template=query_template)
            ref_cat.rename_column('raMean', 'ra')
            ref_cat.rename_column('decMean', 'dec')
            ref_cat.write(fn, overwrite=True)
    else:
        # if you passed a file name, load that
        ref_cat = Table.read(ref_cat)

    if config['model_subtract_star']:
        back_out = os.path.join(out_path, 'star_subtracted')
        utils.mkdir_if_needed(back_out)
        all_files = glob.glob(os.path.join(data_path, f'lbc{bandpass.lower()}*{glob_select}'))

        print('~~~~~~LOOOOOOOOOOOOOOOK HHHHEREEEEEEEEEEEE!!!!!!!!!!!!!!!!!!!!!!!')
        print('all_files: \n',all_files)
        print('data_path: \n', data_path)
        print('glob_select: \n', glob_select)
        print(f'bandpass: \nlbc{bandpass.lower()}')
        print('back_out: \n',back_out)
        print('CONFIG::::::: \n',config)

        for fi in all_files:
            fi_base = fi.split('/')[-1].replace(glob_select,f'starsub_{glob_select}')
            copyfile(fi, os.path.join(back_out,fi_base))
        data_path = back_out
        glob_select = f'starsub_{glob_select}'

    elif config['subtract_background']:
        back_out = os.path.join(out_path, 'back_subtracted')
        utils.mkdir_if_needed(back_out)
        all_files = glob.glob(os.path.join(data_path, f'lbc{bandpass.lower()}*{glob_select}'))
        for fi in all_files:
            fi_base = fi.split('/')[-1].replace(glob_select,f'backsub_{glob_select}')
            copyfile(fi, os.path.join(back_out,fi_base))
        data_path = back_out
        glob_select = f'backsub_{glob_select}'

    # short exposures for finding the astrometric solution
    files_cali = fetch_files(data_path, bandpass, [0, 50], name_must_contain=glob_select)
    num_cali = len(files_cali)

    # long exposures for making science images
    files_sci = fetch_files(data_path, bandpass, [200, 500], name_must_contain=glob_select)
    num_sci = len(files_sci)

    # loop over file types (calibration and science)
    for ftype, files, num in zip(['cali', 'sci'], [files_cali, files_sci], [num_cali, num_sci]):

        if ftype == 'cali' and not config['run_on_cali_files']: continue
        if ftype == 'sci' and not config['run_on_sci_files']: continue
        if ftype == 'cali' : texp = config['texp_cali']
        else: texp = config['texp_sci']

        logger.info(f'Solving astrometry for {bandpass}-band {ftype} frames')

        frame_out = os.path.join(out_path, ftype.lower())
        utils.mkdir_if_needed(frame_out)

        astrom = []
        sky_pos = []
        logger.start_tqdm()

        if config['construct_artpop']:
            options = {'color1' : {'psf' : config['psf_fn_r'], 'zpt' : config['zpt_r'], 'artpop_band' : config['artpop_band_r'], 'artpop_model_fn' : os.path.join(config['artpop_model_dir'],config['artpop_model_fn_r']), 'extinction' : config['extinction_r'], 'color_term' : config['colorterm_r']},
            'color2' : {'psf' : config['psf_fn_b'], 'zpt' : config['zpt_b'], 'artpop_band' : config['artpop_band_b'], 'artpop_model_fn' : os.path.join(config['artpop_model_dir'],config['artpop_model_fn_b']), 'extinction' : config['extinction_b'], 'color_term' : config['colorterm_b']},
            'include_sky_sb': False, 'exposure_time' : texp, 'out_dir' : os.path.join(config['out_dir'],ftype), 'image_dir' : os.path.join(config['out_dir'],ftype)}
        if config['construct_artpop'] and bandpass == 'R':
            logger.info('Constructing ArtPop model...')
            model1, model2, src = artpop_functions.run_artimager(config, options)
        elif not config['construct_artpop']:
            src = None

        # run solve-field on each image
        files_updated = []
        for fn in tqdm(files):

            fn_updated = fn
            if config['model_subtract_star']:
                logger.info('Subtracting star for ' + fn)
                imfit.subtract_bright_star(fn, config, bandpass.lower(), glob_select, use_cutout=True)
                if config['subtract_background']:
                    back_out = os.path.join(out_path, 'back_subtracted')
                    utils.mkdir_if_needed(back_out)
                    fn_base = fn.split('/')[-1].replace(glob_select,f'backsub_{glob_select}')
                    copyfile(fn, os.path.join(back_out,fn_base))
                    data_path = back_out
                    fn_updated = os.path.join(data_path,fn.split('/')[-1].replace(glob_select,f'backsub_{glob_select}'))

            if config['inject_artpop']:
                logger.info('Injecting ArtPop model for ' + fn)
                if bandpass == 'R': artpop_fn = options['color1']['artpop_model_fn']
                else: artpop_fn = options['color2']['artpop_model_fn']
                image = misc.inject_model(fn, artpop_fn, config['ypos_inject'], config['xpos_inject'])
                io.write_pixels(fn, image, header=fits.open(fn)[config['ext']].header)

            if config['subtract_background']:
                logger.info('Subtracting background for ' + fn_updated)
                background_subtraction(filename = fn_updated, config = config, fn_id=back_fn_id)

            logger.info('Solving field for ' + fn_updated)
            solution = solve_field(fn_updated, index_path=index_path,
                                   tmp_path=tmp_path,
                                   target_radec=center,
                                   search_radius=0.5,
                                   identifier='OBJECT')
            fn_base = fn_updated.split('/')[-1].split('.proc.fits')[0]
            utils.mkdir_if_needed(os.path.join(frame_out,'catalogs'))
            cat_fn = os.path.join(frame_out,'catalogs',f'{fn_base}.cat')
            cat = extract_bright_stars(fn_updated,catalog_path=cat_fn)
            check = check_astrometric_solution(ref_cat,
                                               header=solution.header,
                                               cat=cat, max_sep=1)
            tweak = TweakWCS(solution.header, check.cat_match, check.ref_match)
            tweak.optimize()
            tweak.update_header(solution.header)
            tweak.update_header(solution.fitsio_header)

            if config['inject_artpop']:
                w = WCS(solution.header)
                sky = w.pixel_to_world(config['xpos_inject'],config['ypos_inject'])
                sky_pos.append(sky)

            astrom.append(solution)
            if config['make_plots']:
                fig_dir = os.path.join(frame_out, 'diagnostic-plots')
                utils.mkdir_if_needed(fig_dir)
                check = check_astrometric_solution(ref_cat,
                                                   header=solution.header,
                                                   cat=cat, max_sep=1,
                                                   make_plot=config['make_plots'],
                                                   xlim=[-0.3, 0.3],
                                                   ylim=[-0.3, 0.3])
                fig_fn = os.path.basename(fn_updated).replace('.fits',
                                                      '_check_wcs.png')
                check.fig.savefig(os.path.join(fig_dir, fig_fn), dpi=250)

            files_updated.append(fn_updated)

        logger.end_tqdm()

        # resample images
        logger.info(f'Registering and writing {ftype} frames')
        reg_files = []
        for fn, sol in zip(files_updated, astrom):
            resamp = improc.resample_image(
                fn, center[0], center[1], config['pixscale'], config['output_dimensions']['x'],
                config['output_dimensions']['y'], sol.fitsio_header)

            good_keys = ['EXTNAME','DATE_OBS','BLANK','GAIN','RDNOISE','SATURATE','EXPTIME','TEXPTIME','OBJECT','MJD_OBS','UTC_OBS','LST_OBS','PIXSIZE','FILTER','IMAGETYP','CCDSUM','LBCBIN','LBCSPRED','LBCBACK','CCDTEM','LBCOBFIL','FILTEOFF','SUBTRACT_OVERSCAN','TRIM_IMAGE','ORIGIN','FILENAME','OBS_ID','PROPID','OS_NUM','LBCOBID','LBCOBNAM','PARTNER','PI_NAME','AIRMASS','LBTLAT','LBTLONG','LBTELEV','PIXSCAL','DITHSEQ','DITHOFFX','DITHOFFY','TELESCOP','INSTRUME','LBCFWHM','LBTPRES','LBTRHUM','LBTTEMP','LBTWNDIR','LBTWNSPD','GUISTAT','DETECTOR','LBCNCHIP','INSPRE','DETSWV','MIRRORX','MIRRORY','MIRRORZ','MIRRORRX','MIRRORRY','MIRRORRZ','Z4','Z5','Z6','Z7','Z8','Z11','Z22','BACKSUB_TYPE','BACKSUB_VAL']
            header_old = fits.getheader(fn)
            header_new = resamp.wcs.to_header()
            for k, v in header_old.items():
                if k in good_keys:
                    header_new[k] = v

            out_fn = os.path.basename(fn).replace('.fits', '_reg.fits')
            out_fn = os.path.join(frame_out, out_fn)
            fits.writeto(out_fn, resamp.pixels, header_new, overwrite=True)
            reg_files.append(out_fn)

        # Create exposure map
        logger.info(f'Creating exposure map for {ftype} frames')

        stack_data = np.ndarray((num,config['output_dimensions']['y'],config['output_dimensions']['x']))
        exposure_map = np.zeros((config['output_dimensions']['y'],config['output_dimensions']['x']))
        im_num = 0
        logger.start_tqdm()

        for fn in tqdm(reg_files):
            hdul = fits.open(fn)
            exptime = hdul[0].header['EXPTIME']
            stack_data[im_num] = hdul[0].data
            ind = np.nonzero(hdul[0].data)
            exposure_map[ind] += exptime
            hdul.close()
            im_num+=1
        logger.end_tqdm()

        hdu_exp = fits.PrimaryHDU(exposure_map)
        hdu_exp.writeto(os.path.join(frame_out,f'lbc{bandpass.lower()}_M81blob_exposuremap.fits'),overwrite=True)

    if config['subtract_background'] and config['model_subtract_star'] : config['glob_select'] = f'backsub_{glob_select}'

    return config, sky_pos, src


def calibrate_images(config):

    if type(config) is not dict:
        with open(config, 'r') as config:
            config = yaml.load(config, Loader=yaml.FullLoader)

    flux_apers = config['flux_apers']
    num_apers = len(utils.list_of_strings(flux_apers))
    extra_params = f'ALPHA_J2000,DELTA_J2000,ISOAREA_IMAGE,FLUX_APER({num_apers})'
    chosen_aper = config['chosen_aper']

    sci_ims = glob.glob(os.path.join(config['image_dir'],'sci','lbc*reg.fits'))
    cali_ims = glob.glob(os.path.join(config['image_dir'],'cali','lbc*reg.fits'))
    cali_cats_r = []
    sci_cats_r = []
    cali_cats_b = []
    sci_cats_b = []


    # Run SE on images, store catalogs
    logger.info('Creating star catalogs...')
    for ftype, ims in zip(['cali','sci'],[cali_ims,sci_ims]):

        if ftype == 'cali' and not config['run_on_cali_files']: continue
        if ftype == 'sci' and not config['run_on_sci_files']: continue

        for im in ims:
            fn_base = im.split('/')[-1].split('.fits')[0]
            cat_fn = os.path.join(config['image_dir'],f'{ftype}/catalogs/{fn_base}.cat')
            cat = detection.sextractor.run(im,DETECT_MINAREA=3,DETECT_THRESH=5,PIXEL_SCALE=0.225,PHOT_APERTURES=flux_apers,
                                 catalog_path=cat_fn,extra_params=extra_params)

            star_query = 'FLAGS==0 and ISOAREA_IMAGE > 5 and \
                      FWHM_IMAGE > 1 and FWHM_IMAGE < 26'
            cat = cat[cat.to_pandas().query(star_query).index.values]

            cat['filename'] = fn_base+'.fits'

            if 'lbcr' in fn_base:
                if ftype == 'cali': cali_cats_r.append(cat)
                else: sci_cats_r.append(cat)

            else:
                if ftype == 'cali': cali_cats_b.append(cat)
                else: sci_cats_b.append(cat)

    # Sort catalogs, remove some known bad objects
    ## Make aperture plots and color term plots
    # http://spiff.rit.edu/richmond/sne/sn2011fe/color_terms.html
    panstarrs = Table.read(os.path.join(config['image_dir'],config['reference_catalog_withBESSEL']),format='ascii')
    panstarrs_coord = SkyCoord(panstarrs['ra'],panstarrs['dec'],unit='deg')

    size = (51,51)

    cat_all_cali_r = []
    cat_all_cali_b = []
    cat_all_cali_p = []
    cat_all_sci_r = []
    cat_all_sci_b = []
    cat_all_sci_p = []

    logger.info('Sorting catalogs...')
    for ftype, cats in zip(['cali','sci'],[cali_cats_r,sci_cats_r]):

        if ftype == 'cali' and not config['run_on_cali_files']: continue
        if ftype == 'sci' and not config['run_on_sci_files']: continue

        # Get catalogs to compare to
        if ftype == 'cali':
            cats_to_search = cali_cats_b
            texp = config['texp_cali']
        else:
            cats_to_search = sci_cats_b
            texp = config['texp_sci']


        for cat_r in cats:
            cat_b = find_closest_cat(cat_r,cats_to_search)

            r_coord = SkyCoord(cat_r['ALPHA_J2000'],cat_r['DELTA_J2000'],unit='deg')
            b_coord = SkyCoord(cat_b['ALPHA_J2000'],cat_b['DELTA_J2000'],unit='deg')

            ## Match and organize catalogs
            # Find r-band obj that appear in panstarrs
            idx_r_to_p, sep_r_to_p, _ = r_coord.match_to_catalog_sky(panstarrs_coord)
            mask_r_to_p = sep_r_to_p.arcsec < config['allowed_sep']
            obj_in_panstarrs_r = cat_r[mask_r_to_p]

            # Find b-band obj that appear in panstarrs
            idx_b_to_p, sep_b_to_p, _ = b_coord.match_to_catalog_sky(panstarrs_coord)
            mask_b_to_p = sep_b_to_p.arcsec < config['allowed_sep']
            obj_in_panstarrs_b = cat_b[mask_b_to_p]

            # Find smaller of the catalogs matched to panstarrs
            len_b = len(obj_in_panstarrs_b)
            len_r = len(obj_in_panstarrs_r)
            if len_b < len_r:
                obj_in_panstarrs_small = obj_in_panstarrs_b
                obj_in_panstarrs_large = obj_in_panstarrs_r
            else:
                obj_in_panstarrs_small = obj_in_panstarrs_r
                obj_in_panstarrs_large = obj_in_panstarrs_b

            # Match the smaller catalog to the larger catalog (both of which have been matched to panstarrs) and vice versa
            obj_in_panstarrs_small_coord = SkyCoord(obj_in_panstarrs_small['ALPHA_J2000'],obj_in_panstarrs_small['DELTA_J2000'],unit='deg')
            obj_in_panstarrs_large_coord = SkyCoord(obj_in_panstarrs_large['ALPHA_J2000'],obj_in_panstarrs_large['DELTA_J2000'],unit='deg')

            idx, sep, _ = obj_in_panstarrs_small_coord.match_to_catalog_sky(obj_in_panstarrs_large_coord)
            mask = sep.arcsec < config['allowed_sep']
            obj_in_panstarrs_small = obj_in_panstarrs_small[mask]
            obj_in_panstarrs_large = obj_in_panstarrs_large[idx[mask]]

            if len(obj_in_panstarrs_large) == 0: continue

            if 'lbcr' in obj_in_panstarrs_large[0]['filename']:
                obj_in_panstarrs_r = obj_in_panstarrs_large
                obj_in_panstarrs_b = obj_in_panstarrs_small
            else:
                obj_in_panstarrs_b = obj_in_panstarrs_large
                obj_in_panstarrs_r = obj_in_panstarrs_small


            # Now both obj_in_panstars catalogs should have all the same matching objects in the same order
            # Get objects in Pan-STARRS that match to the obj_in_panstarrs catalogs
            obj_in_panstarrs_coord = SkyCoord(obj_in_panstarrs_r['ALPHA_J2000'],obj_in_panstarrs_r['DELTA_J2000'],unit='deg')
            idx, sep, _ = obj_in_panstarrs_coord.match_to_catalog_sky(panstarrs_coord)
            panstarrs_matched = panstarrs[idx]

            # Phew! Catalog cross-matching complete. Now calculate mags and stuff!
            # Get instrumental magnitudes for b- and r-bands
            mag_inst_r = []
            mag_inst_b = []
            for obj_idx in range(len(obj_in_panstarrs_r)):
                file_r = obj_in_panstarrs_r[obj_idx]['filename']
                file_b = obj_in_panstarrs_b[obj_idx]['filename']
                im_hdul_r = fits.open(os.path.join(config['image_dir'],ftype,file_r),ignore_blank=True)
                im_hdul_b = fits.open(os.path.join(config['image_dir'],ftype,file_b),ignore_blank=True)

                position_r = (obj_in_panstarrs_r[obj_idx]['X_IMAGE']-1,obj_in_panstarrs_r[obj_idx]['Y_IMAGE']-1)
                position_b = (obj_in_panstarrs_b[obj_idx]['X_IMAGE']-1,obj_in_panstarrs_b[obj_idx]['Y_IMAGE']-1)

                cutout_r = Cutout2D(im_hdul_r[0].data, position_r, size)
                cutout_b = Cutout2D(im_hdul_b[0].data, position_b, size)

                x_corrected_r = position_r[0]-cutout_r.origin_original[0]
                y_corrected_r = position_r[1]-cutout_r.origin_original[1]
                x_corrected_b = position_b[0]-cutout_b.origin_original[0]
                y_corrected_b = position_b[1]-cutout_b.origin_original[1]
                position_corrected_r = (x_corrected_r,y_corrected_r)
                position_corrected_b = (x_corrected_b,y_corrected_b)

                aperture_r = CircularAperture(position_corrected_r, r=config['chosen_aper'])
                aperture_b = CircularAperture(position_corrected_b, r=config['chosen_aper'])

                phot_table_r = aperture_photometry(cutout_r.data, aperture_r)
                phot_table_b = aperture_photometry(cutout_b.data, aperture_b)

                phot_table_r['aperture_sum'].info.format = '%.8g'  # for consistent table output
                phot_table_b['aperture_sum'].info.format = '%.8g'  # for consistent table output

                aper_sum_r = phot_table_r[0]['aperture_sum']
                aper_sum_b = phot_table_r[0]['aperture_sum']

                mag_inst_r.append(-2.5 * np.log10(aper_sum_r*config['gain']/texp))
                mag_inst_b.append(-2.5 * np.log10(aper_sum_b*config['gain']/texp))
                im_hdul_r.close()
                im_hdul_b.close()
            obj_in_panstarrs_r['mag_inst'] = mag_inst_r
            obj_in_panstarrs_b['mag_inst'] = mag_inst_b

            # Find airmass
            lbt = EarthLocation.of_site('Large Binocular Telescope')
            glob_select = config['glob_select']
            chip = config['chip_num']
            file_b = os.path.join(config['image_dir'],ftype,obj_in_panstarrs_b[0]['filename'].split('-chip')[0]+f'-chip{chip}{glob_select}')
            file_r = os.path.join(config['image_dir'],ftype,obj_in_panstarrs_r[0]['filename'].split('-chip')[0]+f'-chip{chip}{glob_select}')
            hdul_b = fits.open(file_b,ignore_blank=True)
            hdul_r = fits.open(file_r,ignore_blank=True)
            utc_time_r = hdul_r[0].header['DATE_OBS']
            utc_time_b = hdul_b[0].header['DATE_OBS']

            coord_r = SkyCoord(obj_in_panstarrs_r['ALPHA_J2000'],obj_in_panstarrs_r['DELTA_J2000'],unit='deg')
            coord_b = SkyCoord(obj_in_panstarrs_b['ALPHA_J2000'],obj_in_panstarrs_b['DELTA_J2000'],unit='deg')

            obj_in_panstarrs_b['airmass'] = calculate_air_mass(coord_b, utc_time_b, lbt)
            obj_in_panstarrs_r['airmass'] = calculate_air_mass(coord_r, utc_time_r, lbt)
            hdul_b.close()
            hdul_r.close()

            # Re-determine which catalog is which
            if ftype == 'cali':
                if len(cat_all_cali_r) == 0:
                    cat_all_cali_r = obj_in_panstarrs_r
                    cat_all_cali_b = obj_in_panstarrs_b
                    cat_all_cali_p = panstarrs_matched
                else:
                    cat_all_cali_r = vstack([cat_all_cali_r,obj_in_panstarrs_r])
                    cat_all_cali_b = vstack([cat_all_cali_b,obj_in_panstarrs_b])
                    cat_all_cali_p = vstack([cat_all_cali_p,panstarrs_matched])
            else:
                if len(cat_all_sci_r) == 0:
                    cat_all_sci_r = obj_in_panstarrs_r
                    cat_all_sci_b = obj_in_panstarrs_b
                    cat_all_sci_p = panstarrs_matched
                else:
                    cat_all_sci_r = vstack([cat_all_sci_r,obj_in_panstarrs_r])
                    cat_all_sci_b = vstack([cat_all_sci_b,obj_in_panstarrs_b])
                    cat_all_sci_p = vstack([cat_all_sci_p,panstarrs_matched])


        if ftype == 'cali':
            cat_all_cali_r,_ = make_region_files(cat_all_cali_r, ftype, out_dir=os.path.join(config['out_dir'],ftype))
            cat_all_cali_b,_ = make_region_files(cat_all_cali_b, ftype, out_dir=os.path.join(config['out_dir'],ftype))
        else:
            cat_all_sci_r,_ = make_region_files(cat_all_sci_r, ftype, out_dir=os.path.join(config['out_dir'],ftype))
            cat_all_sci_b,_ = make_region_files(cat_all_sci_b, ftype, out_dir=os.path.join(config['out_dir'],ftype))

    # The catalogs should all the sorted now.
    logger.info('Removing bad stars from catalogs...')
    ra_cut = np.array(misc.list_of_floats(config['ra_cut']))
    dec_cut = np.array(misc.list_of_floats(config['dec_cut']))

    if len(ra_cut) != len(dec_cut):
        logger.warning('ra_cut and dec_cut are different lengths. Check cuts and try again.')
        sys.exit('Calibration stopped.')

    for ftype in ['cali','sci']:

        if ftype == 'cali' and not config['run_on_cali_files']: continue
        if ftype == 'sci' and not config['run_on_sci_files']: continue

        ## Remove bad detections by hand
        if ftype == 'cali':
            coord_b = SkyCoord(cat_all_cali_b['ALPHA_J2000'],cat_all_cali_b['DELTA_J2000'],unit='deg')
            coord_cut = SkyCoord(ra_cut,dec_cut,unit='deg')

            idx, sep, _ = coord_b.match_to_catalog_sky(coord_cut)
            mask = sep.arcsec < 1.
            cat_all_cali_b = cat_all_cali_b[~mask]
            cat_all_cali_r = cat_all_cali_r[~mask]
            cat_all_cali_p = cat_all_cali_p[~mask]

        else:
            coord_b = SkyCoord(cat_all_sci_b['ALPHA_J2000'],cat_all_sci_b['DELTA_J2000'],unit='deg')
            coord_cut = SkyCoord(ra_cut,dec_cut,unit='deg')

            idx, sep, _ = coord_b.match_to_catalog_sky(coord_cut)
            mask = sep.arcsec < 1.
            cat_all_sci_b = cat_all_sci_b[~mask]
            cat_all_sci_r = cat_all_sci_r[~mask]
            cat_all_sci_p = cat_all_sci_p[~mask]


    # Aperture correction
    # Here are some hand-picked stars for the correction
    ra_aper = np.array(misc.list_of_floats(config['calibration_stars_ra']))
    dec_aper = np.array(misc.list_of_floats(config['calibration_stars_dec']))

    # This is going to be a repeat of the calculations above, but for all exposures.
    stars_for_aper_corrections = Table()
    stars_for_aper_corrections['ra'] = ra_aper
    stars_for_aper_corrections['dec'] = dec_aper


    for ftype in ['cali','sci']:

        if ftype == 'cali' and not config['run_on_cali_files']: continue
        if ftype == 'sci' and not config['run_on_sci_files']: continue

        logger.info(f'Calculating aperture corrections for {ftype} images...')

        if ftype == 'cali':
            cat_all_r = cat_all_cali_r
            cat_all_b = cat_all_cali_b
            cat_all_p = cat_all_cali_p
            texp = config['texp_cali']
        else:
            cat_all_r = cat_all_sci_r
            cat_all_b = cat_all_sci_b
            cat_all_p = cat_all_sci_p
            texp = config['texp_sci']


        coord_r = SkyCoord(cat_all_r['ALPHA_J2000'],cat_all_r['DELTA_J2000'],unit='deg')
        coord_aper = SkyCoord(ra_aper,dec_aper,unit='deg')

        idx, sep, _ = coord_r.match_to_catalog_sky(coord_aper)
        mask = sep.arcsec < 1.
        cat_magcal_b = cat_all_b[mask]
        cat_magcal_r = cat_all_r[mask]
        cat_magcal_p = cat_all_p[mask]

        logger.info(f'Check catalog lengths (these should be the same; b,r,p): {len(cat_magcal_b)}, {len(cat_magcal_r)}, {len(cat_magcal_p)}')

        # Get aperture sums for all radii, in both r- and b-band
        size = (51,51)
        radii = np.arange(1,25)
        fig, axs = plt.subplots(int(2*(len(cat_magcal_p)+1)/3), 3, figsize=(15,260))
        mask_r = []
        mask_b = []

        row = 0
        col = 0

        for stars, band in zip([cat_magcal_r, cat_magcal_b],['r','b']):
            aper_sums = np.zeros((len(cat_magcal_p),len(radii)))
            star_count = 0

            if band == 'r': color = 'darkred'
            else: color = 'darkblue'
            for star in stars:
                saturation_warning = False
                file = star['filename']
                file_stub = file.split('-chip')[0]
                star_id = star['id']
                im_hdul = fits.open(os.path.join(config['image_dir'],ftype,file),ignore_blank=True)
                position = (star['X_IMAGE']-1,star['Y_IMAGE']-1)

                cutout = Cutout2D(im_hdul[0].data, position, size)

                if np.sum(cutout.data >= config['saturation_limit']) != 0:
                    if band == 'r': mask_r.append(False)
                    else: mask_b.append(False)
                    saturation_warning = True
                else:
                    if band == 'r': mask_r.append(True)
                    else: mask_b.append(True)

                x_corrected = position[0]-cutout.origin_original[0]
                y_corrected = position[1]-cutout.origin_original[1]
                position_corrected = (x_corrected,y_corrected)

                apertures = [CircularAperture(position_corrected, r=r) for r in radii]

                phot_table = aperture_photometry(cutout.data, apertures)
                for phot_col in phot_table.colnames:
                     phot_table[phot_col].info.format = '%.8g'  # for consistent table output
                aper_sums[star_count] = np.array([phot_table[0][f'aperture_sum_{i}'] for i in np.arange(len(radii))])

                norm=colors.LogNorm(vmin=10, vmax=12000)
                axs[row, col].imshow(cutout.data, origin='lower', norm=norm, cmap='gray')
                axs[row, col].set_title(f'{file_stub} | {star_id}')
                for aper in apertures:
                    aper.plot(axs[row, col],color=color)
                if saturation_warning:
                    axs[row, col].text(0,0,'SATURATED',fontsize=30,color='gold')

                col +=1
                if col==3:
                    row+=1
                    col=0

                star_count+=1

                im_hdul.close()
            if band == 'r':
                aper_sums_all_r = aper_sums
            else: aper_sums_all_b = aper_sums


        plt.savefig(os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'aperture_plots_all_exposures.png'))

        # Remove saturated stars
        full_mask = ~(~np.asarray(mask_r)+~np.asarray(mask_b))
        aper_sums_all_r = aper_sums_all_r[full_mask]
        aper_sums_all_b = aper_sums_all_b[full_mask]
        cat_magcal_r = cat_magcal_r[full_mask]
        cat_magcal_b = cat_magcal_b[full_mask]
        cat_magcal_p = cat_magcal_p[full_mask]


        # Remove bad stars from each catalog
        if config['check_aperture_plots']:
            plot_fn = os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'aperture_plots_all_exposures.png')
            request = f'\nCheck aperture plot ({plot_fn}) and input star numbers that you wish to exlude from the calibration measurement.\nNumbers should be comma-deliminated with no spaces.\n\nInput numbers for r-band {ftype} stars in all exposures: '
            cut_stars_r = interactive.get_input(request, anything_acceptable=True, exit_response='stop', full_stop = False)
            request = f'Input numbers for b-band {ftype} stars: '
            cut_stars_b = interactive.get_input(request, anything_acceptable=True, exit_response='stop', full_stop = False)
            cut_stars_r = misc.list_of_floats(cut_stars_r)
            cut_stars_b = misc.list_of_floats(cut_stars_b)
        else:
            cut_stars_r = misc.list_of_floats(config['all_exposures'][ftype]['cut_star_nums_r'])
            cut_stars_b = misc.list_of_floats(config['all_exposures'][ftype]['cut_star_nums_b'])

        mask=[]

        for idx in range(len(cat_magcal_p)):
            keep_star=True
            if cat_magcal_r[idx]['id'] in cut_stars_r:
                keep_star = False
            elif cat_magcal_b[idx]['id'] in cut_stars_b:
                keep_star=False
            mask.append(keep_star)

        cat_magcal_r = cat_magcal_r[mask]
        cat_magcal_b = cat_magcal_b[mask]
        cat_magcal_p = cat_magcal_p[mask]
        aper_sums_all_r = aper_sums_all_r[mask]
        aper_sums_all_b = aper_sums_all_b[mask]

        logger.info(f'Plotting curves of growth for {ftype} images...')
        # Plot curves of growth
        num_exposures = len(np.unique(cat_magcal_r['filename']))
        if num_exposures == 1: num_exposures+=1
        bin_size = 0.02

        fig, axs = plt.subplots(num_exposures,3,figsize=(20,num_exposures*6))

        flux_fracs_all_r = []
        flux_fracs_all_b = []
        flux_levels_all_r = []
        flux_levels_all_b = []
        flux_medians_r = []
        flux_medians_b = []

        files_r = [cat_magcal_r[0]['filename']]
        files_b = [cat_magcal_b[0]['filename']]
        prev_exp = cat_magcal_r[0]['filename']

        row = 0

        for idx in range(len(cat_magcal_p)):

            if cat_magcal_r[idx]['filename'] != prev_exp:
                flux_mean_all_r = np.mean(flux_fracs_all_r)
                flux_median_all_r = np.median(flux_fracs_all_r)
                flux_mean_all_b = np.mean(flux_fracs_all_b)
                flux_median_all_b = np.median(flux_fracs_all_b)

                flux_medians_r.append(flux_median_all_r)
                flux_medians_b.append(flux_median_all_b)
                files_r.append(cat_magcal_r[idx]['filename'])
                files_b.append(cat_magcal_b[idx]['filename'])

                # Find best bins
                min_bin = min(flux_fracs_all_r) - 0.01
                max_bin = max(flux_fracs_all_b) + 0.01

                axs[row][2].set_xlabel(f'Fraction of flux outside a {chosen_aper} pixel radius aperture',fontsize=20)
                axs[row][2].set_ylabel('Number of stars',fontsize=20)
                axs[row][2].hist(np.asarray(flux_fracs_all_r),color='red',bins=np.arange(min_bin,max_bin,bin_size),alpha=0.5,label='r-band');
                axs[row][2].hist(np.asarray(flux_fracs_all_b),color='blue',bins=np.arange(min_bin,max_bin,bin_size),alpha=0.5,label='b-band');
                axs[row][2].vlines(flux_mean_all_r,0,18,color='darkred',lw=2)
                axs[row][2].vlines(flux_mean_all_b,0,18,color='darkblue',lw=2)
                axs[row][2].vlines(flux_median_all_r,0,18,color='darkred',lw=2,linestyle='dashed')
                axs[row][2].vlines(flux_median_all_b,0,18,color='darkblue',lw=2,linestyle='dashed')
                axs[row][2].set_ylim(0,5.25)
                axs[row][2].set_xlim(0.1,0.8)

                flux_fracs_all_r = []
                flux_fracs_all_b = []
                flux_levels_all_r = []
                flux_levels_all_b = []

                row+=1


            file_stub_r = cat_magcal_r[idx]['filename'].split('-chip')[0]
            file_stub_b = cat_magcal_b[idx]['filename'].split('-chip')[0]
            aper_sum_r = aper_sums_all_r[idx]
            aper_sum_b = aper_sums_all_b[idx]

            axs[row][0].set_xlabel('Aperture radius (pixels)',fontsize=25)
            axs[row][0].set_ylabel('Normalized Flux',fontsize=25)
            axs[row][1].set_xlabel('Aperture radius (pixels)',fontsize=25)
            axs[row][1].set_ylabel('Flux (counts)',fontsize=25)

            der_r = np.diff(aper_sum_r) / np.diff(radii)
            der_b = np.diff(aper_sum_b) / np.diff(radii)

            # Find the aperture radius where the derivative crosses zero (or meets other threshold)
            # Estimate the flux at that radius
            if len(np.where(np.diff(np.sign(der_r)))[0]) == 0:
                flux_level_r = max(aper_sum_r)
            else:
                zero_crossing_r = np.where(np.diff(np.sign(der_r)))[0][0]
                flux_level_r = aper_sum_r[zero_crossing_r+1]

            if len(np.where(np.diff(np.sign(der_b)))[0]) == 0:
                flux_level_b = max(aper_sum_b)
            else:
                zero_crossing_b = np.where(np.diff(np.sign(der_b)))[0][0]
                flux_level_b = aper_sum_b[zero_crossing_b+1]


            flux_fraction_b = 1 - (aper_sum_b[np.where(radii==config['chosen_aper'])] / flux_level_b)
            flux_fraction_r = 1 - (aper_sum_r[np.where(radii==config['chosen_aper'])] / flux_level_r)

            flux_fracs_all_r.append(flux_fraction_r)
            flux_levels_all_r.append(flux_level_r)
            flux_fracs_all_b.append(flux_fraction_b)
            flux_levels_all_b.append(flux_level_b)

            axs[row][0].plot(radii,aper_sum_r*config['gain']/(flux_level_r*config['gain']),color='darkred')
            axs[row][0].plot(radii,aper_sum_b*config['gain']/(flux_level_b*config['gain']),color='darkblue')
            axs[row][1].plot(radii,aper_sum_r*config['gain'],color='darkred')
            axs[row][1].plot(radii,aper_sum_b*config['gain'],color='darkblue')
            rad = (radii[:-1] +radii[1:]) / 2
            axs[row][0].plot(rad,der_r/max(der_r),color='darkred')
            axs[row][0].plot(rad,der_b/max(der_b),color='darkblue')

            axs[row][0].axvline(x=config['chosen_aper'],c='black',linestyle='dashed')
            axs[row][1].axvline(x=config['chosen_aper'],c='black',linestyle='dashed')
            axs[row][0].plot(np.arange(25),np.zeros(25),color='black')
            #axs[row][0].legend(fontsize=15,loc='upper right')
            #axs[row][1].legend(fontsize=15,loc='upper right')
            axs[row][0].set_ylim(-0.5,1.25)
            #axs[row][0].set_xlim(0,40)
            #axs[row][1].set_xlim(0,40)
            axs[row][2].set_xlim(0.1,0.8)
            #axs[row][1].set_ylim(0,1e7)
            axs[row][0].text(7,-0.35,f'{file_stub_r}\n{file_stub_b}',fontsize=20)


            prev_exp = cat_magcal_r[idx]['filename']

        flux_mean_all_r = np.mean(flux_fracs_all_r)
        flux_median_all_r = np.median(flux_fracs_all_r)
        flux_mean_all_b = np.mean(flux_fracs_all_b)
        flux_median_all_b = np.median(flux_fracs_all_b)
        flux_medians_r.append(flux_median_all_r)
        flux_medians_b.append(flux_median_all_b)

        # Find best bins
        min_bin = min(flux_fracs_all_r) - 0.01
        max_bin = max(flux_fracs_all_b) + 0.01

        axs[row][2].set_xlabel(f'Fraction of flux outside a {chosen_aper} pixel radius aperture',fontsize=20)
        axs[row][2].set_ylabel('Number of stars',fontsize=20)
        axs[row][2].hist(np.asarray(flux_fracs_all_r),color='red',bins=np.arange(min_bin,max_bin,bin_size),alpha=0.5,label='r-band');
        axs[row][2].hist(np.asarray(flux_fracs_all_b),color='blue',bins=np.arange(min_bin,max_bin,bin_size),alpha=0.5,label='b-band');
        axs[row][2].vlines(flux_mean_all_r,0,18,color='darkred',lw=2)
        axs[row][2].vlines(flux_mean_all_b,0,18,color='darkblue',lw=2)
        axs[row][2].vlines(flux_median_all_r,0,18,color='darkred',lw=2,linestyle='dashed')
        axs[row][2].vlines(flux_median_all_b,0,18,color='darkblue',lw=2,linestyle='dashed')
        axs[row][2].set_ylim(0,5.25);

        plt.savefig(os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'curveofgrowth_separate_exposures_{chosen_aper}pix.png'))

        # Compare median flux results for different exposures
        bin_min = min([min(flux_medians_r),min(flux_medians_b)])-bin_size
        bin_max = max([max(flux_medians_r),max(flux_medians_b)])+bin_size

        plt.hist(flux_medians_r,color='red',alpha=0.5,bins=np.arange(bin_min,bin_max,bin_size));
        plt.hist(flux_medians_b,color='blue',alpha=0.5,bins=np.arange(bin_min,bin_max,bin_size));
        plt.vlines(np.median(flux_medians_r),0,4,color='darkred',lw=2,linestyle='dashed')
        plt.vlines(np.median(flux_medians_b),0,4,color='darkblue',lw=2,linestyle='dashed')
        plt.xlabel(f'Median fraction of flux outside {chosen_aper}-pix aperture for all exposures',fontsize=20)
        plt.text(.44,3.5,f'median median r-band: {round(np.median(flux_medians_r),4)}\nmedian median b-band: {round(np.median(flux_medians_b),4)}',fontsize=20)

        plt.savefig(os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'flux_fraction_median_separate_exposures_{chosen_aper}pix_withoutoutlier.png'))

        print('Median median r-band: ',np.median(flux_medians_r),'\nMedian median b-band: ',np.median(flux_medians_b))

        cat_all_r,cat_all_b,cat_all_p,zp_r,zp_b,color_term_r,color_term_b = apply_aper_color_corrections(os.path.join(config['image_dir'],ftype),cat_all_r,cat_all_b,cat_all_p,'mag_inst','mag_calibrated',chosen_aper,flux_medians_r,flux_medians_b,files_r,files_b,config['gain'], texp, config['saturation_limit'])
        plot_residuals(cat_all_r,cat_all_b,cat_all_p,'mag_calibrated',save_figs=True,fig_name=os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'residuals_after_aper_color_corrections_{chosen_aper}pix_separate_exposures_mag_inst_fixed.png'))


        weight_exp = True
        mask_std = 2.5
        files_r = np.asarray(np.unique(cat_all_r['filename']))
        files_b = np.asarray(np.unique(cat_all_b['filename']))


        ims_r = np.ndarray((len(files_r),config['output_dimensions']['y'],config['output_dimensions']['x']))
        ims_b = np.ndarray((len(files_b),config['output_dimensions']['y'],config['output_dimensions']['x']))
        noise_r = ims_r*0.
        noise_b = ims_b*0.
        weight_r = np.zeros(len(files_r))
        weight_b = np.zeros(len(files_b))
        x,y = np.mgrid[:2048, :4608]

        # Get exposure map
        logger.info('Creating exposure map...')
        exp_map_r = fits.open(os.path.join(config['image_dir'],ftype,'lbcr_M81blob_exposuremap.fits'),ignore_blank=True)[0].data
        exp_map_b = fits.open(os.path.join(config['image_dir'],ftype,'lbcb_M81blob_exposuremap.fits'),ignore_blank=True)[0].data

        overlap_r = exp_map_r == np.ndarray.max(exp_map_r)
        overlap_b = exp_map_b == np.ndarray.max(exp_map_b)

        # Get image data, background levels
        i=0
        avg_zp_r = np.mean(zp_r['zp'][0:-1])
        avg_zp_b = np.mean(zp_b['zp'][0:-1])
        back_vals_r = []
        back_vals_b = []

        for im_r,im_b in zip(files_r,files_b):
            hdu_r = fits.open(os.path.join(config['image_dir'],ftype,im_r),ignore_blank=True)
            hdu_b = fits.open(os.path.join(config['image_dir'],ftype,im_b),ignore_blank=True)
            zp_im_r = zp_r['zp'][zp_r['filename']==im_r][0]
            zp_im_b = zp_b['zp'][zp_b['filename']==im_b][0]
            data_r = hdu_r[config['ext']].data*overlap_r
            data_b = hdu_b[config['ext']].data*overlap_b
            hdu_r[config['ext']].header['ZEROPOINT'] = -zp_im_r
            hdu_b[config['ext']].header['ZEROPOINT'] = -zp_im_b
            hdu_r.writeto(os.path.join(config['image_dir'],ftype,im_r),overwrite=True)
            hdu_b.writeto(os.path.join(config['image_dir'],ftype,im_b),overwrite=True)

            # Normalize data to common photometric zeropoint
            ims_r[i] = data_r * 10**(-(avg_zp_r-zp_im_r)/2.5)
            ims_b[i] = data_b * 10**(-(avg_zp_b-zp_im_b)/2.5)

            # Get background levels
            if weight_exp:
                '''
                c0_0_r = hdu_r.header['BACKSUB_PARAMVALS_0']
                c0_0_b = hdu_b.header['BACKSUB_PARAMVALS_0']
                c1_0_r = hdu_r.header['BACKSUB_PARAMVALS_1']
                c1_0_b = hdu_b.header['BACKSUB_PARAMVALS_1']
                c0_1_r = hdu_r.header['BACKSUB_PARAMVALS_2']
                c0_1_b = hdu_b.header['BACKSUB_PARAMVALS_2']
                back_medians_r.append(hdu_r.header['BACKSUB_MEDIAN'])
                back_medians_b.append(hdu_b.header['BACKSUB_MEDIAN'])

                weight_r[i] = 1.0/np.mean(c0_0_r + (c1_0_r*x) + (c0_1_r*y))
                weight_b[i] = 1.0/np.mean(c0_0_b + (c1_0_b*x) + (c0_1_b*y))
                '''

                back_vals_r.append(hdu_r[config['ext']].header['BACKSUB_VAL'])
                back_vals_b.append(hdu_b[config['ext']].header['BACKSUB_VAL'])
                weight_r[i] = 1.0/back_vals_r[i]
                weight_b[i] = 1.0/back_vals_b[i]

            else:
                weight_r = None
                weight_b = None

            i+=1

        # Create median stack
        if config['stack_images']:
            logger.info(f'Stacking {ftype} images...')
            median_stack_r = np.median(ims_r,axis=0)
            median_stack_b = np.median(ims_b,axis=0)
            std_stack_r = np.std(ims_r,axis=0)
            std_stack_b = np.std(ims_b,axis=0)
            madstd_stack_r = mad_std(ims_r,axis=0)
            madstd_stack_b = mad_std(ims_b,axis=0)

            # Calculate noise
            noise_r = np.sqrt((median_stack_r*config['gain'])+(np.mean(back_vals_r)*config['gain']+(config['readnoise']**2)))/config['gain']
            noise_b = np.sqrt((median_stack_b*config['gain'])+(np.mean(back_vals_b)*config['gain']+(config['readnoise']**2)))/config['gain']

            hdu = fits.PrimaryHDU(noise_r)
            hdu.writeto(os.path.join(config['image_dir'],ftype,f'noise_{ftype}_r.fits'), overwrite=True)
            hdu = fits.PrimaryHDU(noise_b)
            hdu.writeto(os.path.join(config['image_dir'],ftype,f'noise_{ftype}_b.fits'), overwrite=True)

            # Mask images
            ims_r = np.ma.asarray(ims_r)
            ims_b = np.ma.asarray(ims_b)

            for i in range(len(files_r)):

                stddiff = abs(median_stack_r - ims_r[i]) / noise_r
                mask = np.array(stddiff > 3.0)
                ims_r[i,(mask)] = np.ma.masked

            for i in range(len(files_b)):
                stddiff = abs(median_stack_b - ims_b[i]) / noise_b
                mask = np.array(stddiff > 3.0)
                ims_b[i][(mask)] = np.ma.masked

            # Create mean stack
            mean_stack_r = np.ma.average(ims_r,axis=0,weights=weight_r)
            mean_stack_b = np.ma.average(ims_b,axis=0,weights=weight_b)

            mean_stack_r = np.asarray(mean_stack_r.data)
            mean_stack_b = np.asarray(mean_stack_b.data)

            hdu_r[config['ext']].header['SKY_MEAN'] = np.mean(back_vals_r)
            hdu_b[config['ext']].header['SKY_MEAN'] = np.mean(back_vals_b)

            '''
            hdu_r.header.pop('BACKSUB_PARAMS_0')
            hdu_b.header.pop('BACKSUB_PARAMS_0')
            hdu_r.header.pop('BACKSUB_PARAMS_1')
            hdu_b.header.pop('BACKSUB_PARAMS_1')
            hdu_r.header.pop('BACKSUB_PARAMS_2')
            hdu_b.header.pop('BACKSUB_PARAMS_2')
            hdu_r.header.pop('BACKSUB_PARAMVALS_0')
            hdu_b.header.pop('BACKSUB_PARAMVALS_0')
            hdu_r.header.pop('BACKSUB_PARAMVALS_1')
            hdu_b.header.pop('BACKSUB_PARAMVALS_1')
            hdu_r.header.pop('BACKSUB_PARAMVALS_2')
            hdu_b.header.pop('BACKSUB_PARAMVALS_2')
            hdu_r.header.pop('BACKSUB_MEDIAN')
            hdu_b.header.pop('BACKSUB_MEDIAN')
            '''

            mean_stack_r,num_dead_pix_r = replace_dead_pixels(mean_stack_r, mask=~overlap_r)
            mean_stack_b,num_dead_pix_b = replace_dead_pixels(mean_stack_b, mask=~overlap_b)

            hdu = fits.PrimaryHDU(mean_stack_r,header=hdu_r[config['ext']].header)
            hdu.writeto(os.path.join(config['image_dir'],ftype,config['final_stack_fn_r']),overwrite=True)
            hdu = fits.PrimaryHDU(mean_stack_b,header=hdu_b[config['ext']].header)
            hdu.writeto(os.path.join(config['image_dir'],ftype,config['final_stack_fn_b']),overwrite=True)

        if ftype == 'cali':
            cat_all_cali_r = cat_all_r
            cat_all_cali_b = cat_all_b
            cat_all_cali_p = cat_all_p
        else:
            cat_all_sci_r = cat_all_r
            cat_all_sci_b = cat_all_b
            cat_all_sci_p = cat_all_p

    if config['run_on_sci_files'] and config['construct_psf']:
        logger.info('Calculating PSF...')
        # Get star for PSF
        psf_star = SkyCoord(ra=[config['psf_star']['ra']], dec=[config['psf_star']['dec']],unit='deg')
        catalog = SkyCoord(ra=cat_all_sci_r['ALPHA_J2000'], dec=cat_all_sci_r['DELTA_J2000'],unit='deg')
        idx, _, _ = psf_star.match_to_catalog_sky(catalog)
        psf_star = catalog[idx]
        #print('psf_star: ', psf_star)
        # Create cutout
        #mask = cat_all_sci_b['id'] == psf_star['id']

        #print('mask: ', mask)
        star_r = cat_all_sci_r[idx]#[mask]
        star_b = cat_all_sci_b[idx]#[mask]

        # Calculate median PSF
        # Get median PSF of that star from all exposures
        star_coord = SkyCoord(star_b['ALPHA_J2000'],star_b['DELTA_J2000'],unit='deg')

        files_r = np.unique(np.asarray(cat_all_sci_r['filename']))
        files_b = np.unique(np.asarray(cat_all_sci_b['filename']))
        psfs_r = []
        psfs_b = []

        for file_r, file_b in zip(files_r,files_b):

            hdul_r = fits.open(os.path.join(config['image_dir'],'sci',file_r),ignore_blank=True)
            hdul_b = fits.open(os.path.join(config['image_dir'],'sci',file_b),ignore_blank=True)
            wcs_r = WCS(hdul_r[0].header)
            wcs_b = WCS(hdul_b[0].header)
            cutout_r = Cutout2D(hdul_r[0].data, star_coord, size, wcs=wcs_r)
            cutout_b = Cutout2D(hdul_b[0].data, star_coord, size, wcs=wcs_b)

            psf_r = cutout_r.data/np.sum(cutout_r.data)
            psf_b = cutout_b.data/np.sum(cutout_b.data)

            psfs_r.append(psf_r)
            psfs_b.append(psf_b)

        psf_median_r = np.median(psfs_r,axis=0)
        psf_median_b = np.median(psfs_b,axis=0)

        hdu_r = fits.PrimaryHDU(psf_median_r)
        hdu_r.writeto(os.path.join(config['image_dir'],'sci',config['psf_fn_r']),overwrite=True)
        hdu_b = fits.PrimaryHDU(psf_median_b)
        hdu_b.writeto(os.path.join(config['image_dir'],'sci',config['psf_fn_b']),overwrite=True)
        '''
        fig,ax = plt.subplots(1,2,figsize=(12,6))
        ax[0].imshow(psf_median_r,origin='lower',cmap='gray')
        ax[0].set_title('Median r-band PSF')
        ax[1].imshow(psf_median_b,origin='lower',cmap='gray')
        ax[1].set_title('Median b-band PSF')
        '''

    # Do everthything again but for stacked images...
    logger.info('Repeating calibration for stacked images...')
    for ftype in ['sci','cali']:

        if ftype == 'cali':
            texp = config['texp_cali']
            if not config['run_on_cali_files']: continue
        if ftype == 'sci':
            texp = config['texp_sci']
            if not config['run_on_sci_files']: continue

        # Run SE on stacks
        extra_params = f'ALPHA_J2000,DELTA_J2000,ISOAREA_IMAGE'
        stack_r = os.path.join(config['image_dir'],ftype,config['final_stack_fn_r'])
        stack_b = os.path.join(config['image_dir'],ftype,config['final_stack_fn_b'])

        # Run SE on images, store catalogs
        for im in [stack_r,stack_b]:

            fn_base = im.split('/')[-1].split('.fits')[0]
            cat_fn = os.path.join(config['image_dir'],ftype,'catalogs',f'{fn_base}.cat')

            cat = detection.sextractor.run(im,DETECT_MINAREA=3,DETECT_THRESH=5,PIXEL_SCALE=0.225,
                                 catalog_path=cat_fn,extra_params=extra_params)

            star_query = 'FLAGS==0 and ISOAREA_IMAGE > 5 and \
                      FWHM_IMAGE > 1 and FWHM_IMAGE < 26'
            cat = cat[cat.to_pandas().query(star_query).index.values]

            cat['filename'] = fn_base+'.fits'

            if config['stack_color_id_r'] in fn_base:
                stack_cat_r = cat
            else:
                stack_cat_b = cat
            print('BASE: ',fn_base)
        # Clean catalogs, match to Pan-STARRS

        # Sort catalogs, remove some known bad objects
        panstarrs = Table.read(config['reference_catalog_withBESSEL'],format='ascii')
        panstarrs_coord = SkyCoord(panstarrs['ra'],panstarrs['dec'],unit='deg')

        r_coord = SkyCoord(stack_cat_r['ALPHA_J2000'],stack_cat_r['DELTA_J2000'],unit='deg')
        b_coord = SkyCoord(stack_cat_b['ALPHA_J2000'],stack_cat_b['DELTA_J2000'],unit='deg')

        ## Match and organize catalogs
        # Find r-band obj that appear in panstarrs
        idx_r_to_p, sep_r_to_p, _ = r_coord.match_to_catalog_sky(panstarrs_coord)
        mask_r_to_p = sep_r_to_p.arcsec < config['allowed_sep']
        obj_in_panstarrs_r = stack_cat_r[mask_r_to_p]
        print('Matched panstarrs catalog length: ',len(obj_in_panstarrs_r))
        # Find b-band obj that appear in panstarrs
        idx_b_to_p, sep_b_to_p, _ = b_coord.match_to_catalog_sky(panstarrs_coord)
        mask_b_to_p = sep_b_to_p.arcsec < config['allowed_sep']
        obj_in_panstarrs_b = stack_cat_b[mask_b_to_p]

        # Find smaller of the catalogs matched to panstarrs
        len_b = len(obj_in_panstarrs_b)
        len_r = len(obj_in_panstarrs_r)
        if len_b < len_r:
            obj_in_panstarrs_small = obj_in_panstarrs_b
            obj_in_panstarrs_large = obj_in_panstarrs_r
        else:
            obj_in_panstarrs_small = obj_in_panstarrs_r
            obj_in_panstarrs_large = obj_in_panstarrs_b

        # Match the smaller catalog to the larger catalog (both of which have been matched to panstarrs) and vice versa
        obj_in_panstarrs_small_coord = SkyCoord(obj_in_panstarrs_small['ALPHA_J2000'],obj_in_panstarrs_small['DELTA_J2000'],unit='deg')
        obj_in_panstarrs_large_coord = SkyCoord(obj_in_panstarrs_large['ALPHA_J2000'],obj_in_panstarrs_large['DELTA_J2000'],unit='deg')

        idx, sep, _ = obj_in_panstarrs_small_coord.match_to_catalog_sky(obj_in_panstarrs_large_coord)
        mask = sep.arcsec < config['allowed_sep']
        obj_in_panstarrs_small = obj_in_panstarrs_small[mask]
        obj_in_panstarrs_large = obj_in_panstarrs_large[idx[mask]]
        print('Matched panstarrs catalog length: ',len(obj_in_panstarrs_large))
        if '_r' in obj_in_panstarrs_large[0]['filename']:
            obj_in_panstarrs_r = obj_in_panstarrs_large
            obj_in_panstarrs_b = obj_in_panstarrs_small
        else:
            obj_in_panstarrs_b = obj_in_panstarrs_large
            obj_in_panstarrs_r = obj_in_panstarrs_small

        # Now both obj_in_panstars catalogs should have all the same matching objects in the same order
        # Get objects in Pan-STARRS that match to the obj_in_panstarrs catalogs
        obj_in_panstarrs_coord = SkyCoord(obj_in_panstarrs_r['ALPHA_J2000'],obj_in_panstarrs_r['DELTA_J2000'],unit='deg')
        idx, sep, _ = obj_in_panstarrs_coord.match_to_catalog_sky(panstarrs_coord)
        panstarrs_matched = panstarrs[idx]

        stack_cat_r = obj_in_panstarrs_r
        stack_cat_b = obj_in_panstarrs_b
        stack_cat_p = panstarrs_matched

        # Phew! Catalog cross-matching complete.

        # Get instrumental magnitudes for b- and r-bands
        mag_inst_r = []
        mag_inst_b = []
        im_hdul_r = fits.open(os.path.join(config['image_dir'],ftype,config['final_stack_fn_r']))
        im_hdul_b = fits.open(os.path.join(config['image_dir'],ftype,config['final_stack_fn_b']))
        for obj_idx in range(len(stack_cat_r)):

            position_r = (stack_cat_r[obj_idx]['X_IMAGE']-1,stack_cat_r[obj_idx]['Y_IMAGE']-1)
            position_b = (stack_cat_b[obj_idx]['X_IMAGE']-1,stack_cat_b[obj_idx]['Y_IMAGE']-1)

            cutout_r = Cutout2D(im_hdul_r[0].data, position_r, size)
            cutout_b = Cutout2D(im_hdul_b[0].data, position_b, size)

            x_corrected_r = position_r[0]-cutout_r.origin_original[0]
            y_corrected_r = position_r[1]-cutout_r.origin_original[1]
            x_corrected_b = position_b[0]-cutout_b.origin_original[0]
            y_corrected_b = position_b[1]-cutout_b.origin_original[1]
            position_corrected_r = (x_corrected_r,y_corrected_r)
            position_corrected_b = (x_corrected_b,y_corrected_b)

            aperture_r = CircularAperture(position_corrected_r, r=chosen_aper)
            aperture_b = CircularAperture(position_corrected_b, r=chosen_aper)

            phot_table_r = aperture_photometry(cutout_r.data, aperture_r)
            phot_table_b = aperture_photometry(cutout_b.data, aperture_b)

            phot_table_r['aperture_sum'].info.format = '%.8g'  # for consistent table output
            phot_table_b['aperture_sum'].info.format = '%.8g'  # for consistent table output

            aper_sum_r = phot_table_r[0]['aperture_sum']
            aper_sum_b = phot_table_r[0]['aperture_sum']

            mag_inst_r.append(-2.5 * np.log10(aper_sum_r*config['gain']/texp))
            mag_inst_b.append(-2.5 * np.log10(aper_sum_b*config['gain']/texp))
        im_hdul_r.close()
        im_hdul_b.close()
        stack_cat_r['mag_inst'] = mag_inst_r
        stack_cat_b['mag_inst'] = mag_inst_b

        ## Remove bad detections by hand
        ra_cut = np.array(misc.list_of_floats(config['ra_cut']))
        ra_cut = np.array(misc.list_of_floats(config['dec_cut']))
        coord_b = SkyCoord(stack_cat_b['ALPHA_J2000'],stack_cat_b['DELTA_J2000'],unit='deg')
        coord_cut = SkyCoord(ra_cut,dec_cut,unit='deg')

        idx, sep, _ = coord_b.match_to_catalog_sky(coord_cut)
        mask = sep.arcsec < config['allowed_sep']
        stack_cat_b = stack_cat_b[~mask]
        stack_cat_r = stack_cat_r[~mask]
        stack_cat_p = stack_cat_p[~mask]

        # Do aperture corrections

        # Here are some hand-picked stars for the correction
        ra_aper = np.array(misc.list_of_floats(config['calibration_stars_ra']))
        dec_aper = np.array(misc.list_of_floats(config['calibration_stars_dec']))

        # This is going to be a repeat of the calculations above, but for all exposures.
        stars_for_aper_corrections = Table()
        stars_for_aper_corrections['ra'] = ra_aper
        stars_for_aper_corrections['dec'] = dec_aper

        coord_r = SkyCoord(stack_cat_r['ALPHA_J2000'],stack_cat_r['DELTA_J2000'],unit='deg')
        coord_aper = SkyCoord(ra_aper,dec_aper,unit='deg')

        idx, sep, _ = coord_r.match_to_catalog_sky(coord_aper)
        mask = sep.arcsec < config['allowed_sep']
        cat_magcal_stack_b = stack_cat_b[mask]
        cat_magcal_stack_r = stack_cat_r[mask]
        cat_magcal_stack_p = stack_cat_p[mask]

        logger.info(f'Check catalog lengths (these should be the same; b,r,p): {len(cat_magcal_stack_b)}, {len(cat_magcal_stack_r)}, {len(cat_magcal_stack_p)}')

        # Get aperture sums for all radii, in both r- and b-band
        fig, axs = plt.subplots(int(2*len(cat_magcal_stack_p)/3)+1, 3, figsize=(15,25))
        mask_r = []
        mask_b = []

        row = 0
        col = 0
        for stars, band in zip([cat_magcal_stack_r, cat_magcal_stack_b],['r','b']):
            aper_sums = np.zeros((len(cat_magcal_stack_p),len(radii)))
            star_count = 0
            if band == 'r':
                color = 'darkred'
                im_hdul = fits.open(stack_r)
            else:
                color = 'darkblue'
                im_hdul = fits.open(stack_b)
            for star in stars:
                saturation_warning = False
                position = (star['X_IMAGE']-1,star['Y_IMAGE']-1)

                cutout = Cutout2D(im_hdul[0].data, position, size)

                if np.sum(cutout.data >= config['saturation_limit']) != 0:
                    if band == 'r': mask_r.append(False)
                    else: mask_b.append(False)
                    saturation_warning = True
                else:
                    if band == 'r': mask_r.append(True)
                    else: mask_b.append(True)

                x_corrected = position[0]-cutout.origin_original[0]
                y_corrected = position[1]-cutout.origin_original[1]
                position_corrected = (x_corrected,y_corrected)

                apertures = [CircularAperture(position_corrected, r=r) for r in radii]

                phot_table = aperture_photometry(cutout.data, apertures)
                for phot_col in phot_table.colnames:
                     phot_table[phot_col].info.format = '%.8g'  # for consistent table output
                aper_sums[star_count] = np.array([phot_table[0][f'aperture_sum_{i}'] for i in np.arange(len(radii))])

                norm=colors.LogNorm(vmin=10, vmax=12000)
                axs[row, col].imshow(cutout.data, origin='lower', norm=norm, cmap='gray')
                axs[row, col].set_title(f'{star_count}')
                for aper in apertures:
                    aper.plot(axs[row, col],color=color)
                if saturation_warning:
                    axs[row, col].text(0,0,'SATURATED',fontsize=30,color='gold')

                col +=1
                if col==3:
                    row+=1
                    col=0

                star_count+=1

            if band == 'r':
                aper_sums_stack_r = aper_sums
            else: aper_sums_stack_b = aper_sums


        plt.savefig(os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'aperture_plots_stacks.png'))

        # Remove saturated stars
        full_mask = ~(~np.asarray(mask_r)+~np.asarray(mask_b))
        aper_sums_stack_r = aper_sums_stack_r[full_mask]
        aper_sums_stack_b = aper_sums_stack_b[full_mask]
        cat_magcal_stack_r = cat_magcal_stack_r[full_mask]
        cat_magcal_stack_b = cat_magcal_stack_b[full_mask]
        cat_magcal_stack_p = cat_magcal_stack_p[full_mask]

        # Remove bad stars from each catalog
        if config['check_aperture_plots']:
            plot_fn = os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'aperture_plots_stacks.png')
            request = f'\nCheck aperture plot ({plot_fn}) and input star numbers that you wish to exlude from the calibration measurement.\nNumbers should be comma-deliminated with no spaces.\n\nInput numbers for r-band {ftype} stars in the stacked image: '
            cut_stars_r = interactive.get_input(request, anything_acceptable=True, exit_response='stop', full_stop = False)
            request = f'Input numbers for b-band {ftype} stars: '
            cut_stars_b = interactive.get_input(request, anything_acceptable=True, exit_response='stop', full_stop = False)
            cut_stars_r = misc.list_of_floats(cut_stars_r)
            cut_stars_b = misc.list_of_floats(cut_stars_b)
        else:
            cut_stars_r = misc.list_of_floats(config['stack'][ftype]['cut_star_nums_r'])
            cut_stars_b = misc.list_of_floats(config['stack'][ftype]['cut_star_nums_b'])

        mask=[]

        for idx in range(len(cat_magcal_stack_p)):
            keep_star=True
            if idx in cut_stars_r:
                keep_star=False
            elif idx in cut_stars_b:
                keep_star = False
            mask.append(keep_star)

        cat_magcal_stack_r = cat_magcal_stack_r[mask]
        cat_magcal_stack_b = cat_magcal_stack_b[mask]
        cat_magcal_stack_p = cat_magcal_stack_p[mask]
        aper_sums_stack_r = aper_sums_stack_r[mask]
        aper_sums_stack_b = aper_sums_stack_b[mask]

        # Plot curves of growth
        bin_size = 0.02

        fig, axs = plt.subplots(1,3,figsize=(20,6))

        flux_fracs_stack_r = []
        flux_fracs_stack_b = []
        flux_levels_stack_r = []
        flux_levels_stack_b = []


        for idx in range(len(cat_magcal_stack_p)):

            aper_sum_r = aper_sums_stack_r[idx]
            aper_sum_b = aper_sums_stack_b[idx]

            axs[0].set_xlabel('Aperture radius (pixels)',fontsize=25)
            axs[0].set_ylabel('Normalized Flux',fontsize=25)
            axs[1].set_xlabel('Aperture radius (pixels)',fontsize=25)
            axs[1].set_ylabel('Flux (counts)',fontsize=25)

            der_r = np.diff(aper_sum_r) / np.diff(radii)
            der_b = np.diff(aper_sum_b) / np.diff(radii)

            # Find the aperture radius where the derivative crosses zero (or meets other threshold)
            # Estimate the flux at that radius
            if len(np.where(np.diff(np.sign(der_r)))[0]) == 0:
                flux_level_r = max(aper_sum_r)
            else:
                zero_crossing_r = np.where(np.diff(np.sign(der_r)))[0][0]
                flux_level_r = aper_sum_r[zero_crossing_r+1]

            if len(np.where(np.diff(np.sign(der_b)))[0]) == 0:
                flux_level_b = max(aper_sum_b)
            else:
                zero_crossing_b = np.where(np.diff(np.sign(der_b)))[0][0]
                flux_level_b = aper_sum_b[zero_crossing_b+1]


            flux_fraction_b = 1 - (aper_sum_b[np.where(radii==chosen_aper)] / flux_level_b)
            flux_fraction_r = 1 - (aper_sum_r[np.where(radii==chosen_aper)] / flux_level_r)

            flux_fracs_stack_r.append(flux_fraction_r)
            flux_levels_stack_r.append(flux_level_r)
            flux_fracs_stack_b.append(flux_fraction_b)
            flux_levels_stack_b.append(flux_level_b)

            axs[0].plot(radii,aper_sum_r*config['gain']/(flux_level_r*config['gain']),color='darkred')
            axs[0].plot(radii,aper_sum_b*config['gain']/(flux_level_b*config['gain']),color='darkblue')
            axs[1].plot(radii,aper_sum_r*config['gain'],color='darkred')
            axs[1].plot(radii,aper_sum_b*config['gain'],color='darkblue')
            rad = (radii[:-1] +radii[1:]) / 2
            axs[0].plot(rad,der_r/max(der_r),color='darkred')
            axs[0].plot(rad,der_b/max(der_b),color='darkblue')
            axs[0].axvline(x=chosen_aper,c='black',linestyle='dashed')
            axs[1].axvline(x=chosen_aper,c='black',linestyle='dashed')
            axs[0].plot(np.arange(25),np.zeros(25),color='black')
            axs[0].set_ylim(-0.5,1.25)
            axs[2].set_xlim(0.1,0.8)
            #axs[1].set_ylim(0,1e7)

        flux_mean_stack_r = np.mean(flux_fracs_stack_r)
        flux_median_stack_r = np.median(flux_fracs_stack_r)
        flux_mean_stack_b = np.mean(flux_fracs_stack_b)
        flux_median_stack_b = np.median(flux_fracs_stack_b)
        flux_medians_stack = [flux_median_stack_r,flux_median_stack_b]


        # Find best bins
        min_bin = min(flux_fracs_stack_r) - 0.01
        max_bin = max(flux_fracs_stack_b) + 0.01

        axs[2].set_xlabel(f'Fraction of flux outside a {chosen_aper} pixel radius aperture',fontsize=20)
        axs[2].set_ylabel('Number of stars',fontsize=20)
        axs[2].hist(np.asarray(flux_fracs_stack_r),color='red',bins=np.arange(min_bin,max_bin,bin_size),alpha=0.5,label='r-band');
        axs[2].hist(np.asarray(flux_fracs_stack_b),color='blue',bins=np.arange(min_bin,max_bin,bin_size),alpha=0.5,label='b-band');
        axs[2].vlines(flux_mean_stack_r,0,18,color='darkred',lw=2)
        axs[2].vlines(flux_mean_stack_b,0,18,color='darkblue',lw=2)
        axs[2].vlines(flux_median_stack_r,0,18,color='darkred',lw=2,linestyle='dashed')
        axs[2].vlines(flux_median_stack_b,0,18,color='darkblue',lw=2,linestyle='dashed')
        axs[2].text(0.4,4,f'r-band median: {round(flux_median_stack_r,2)}\nb-band median: {round(flux_median_stack_b,2)}',fontsize=20)
        axs[2].set_ylim(0,5.25);

        plt.savefig(os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'curveofgrowth_stack_{chosen_aper}pix.png'))

        flux_medians_r = [flux_median_stack_r]
        flux_medians_b = [flux_median_stack_b]
        files_r = [config['final_stack_fn_r']]
        files_r = [config['final_stack_fn_b']]

        stack_cat_r,stack_cat_b,stack_cat_p,zp_stack_r,zp_stack_b,color_term_stack_r, color_term_stack_b = apply_aper_color_corrections(os.path.join(config['image_dir'],ftype),stack_cat_r,stack_cat_b,stack_cat_p,'mag_inst','mag_calibrated',chosen_aper,flux_medians_r,flux_medians_b,files_r,files_b,config['gain'], texp, config['saturation_limit'])
        plot_residuals(stack_cat_r,stack_cat_b,stack_cat_p,'mag_calibrated',save_figs=True,fig_name=os.path.join(config['image_dir'],ftype,'diagnostic-plots',f'residuals_after_aper_color_corrections_{chosen_aper}pix_stack.png'))

        with fits.open(os.path.join(config['out_dir'],ftype,zp_stack_r['filename'][0]), 'update') as hdu:
            hdu[config['ext']].header['STACK_ZPT'] = -zp_stack_r['zp'][0]
            hdu[config['ext']].header['STACK_COLORTERM_R'] = color_term_stack_r
            hdu[config['ext']].header['STACK_COLORTERM_B'] = color_term_stack_b
        with fits.open(os.path.join(config['out_dir'],ftype,zp_stack_b['filename'][0]), 'update') as hdu:
            hdu[config['ext']].header['STACK_ZPT'] = -zp_stack_b['zp'][0]
            hdu[config['ext']].header['STACK_COLORTERM_R'] = color_term_stack_r
            hdu[config['ext']].header['STACK_COLORTERM_B'] = color_term_stack_b

        logger.info('Calibration Results: \n')
        print('Stack zeropoints:')
        print(zp_stack_r['filename'][0], -zp_stack_r['zp'][0])
        print(zp_stack_b['filename'][0],-zp_stack_b['zp'][0])
        print('\nColor terms:')
        print('r-band: ',color_term_stack_r)
        print('b-band: ',color_term_stack_b)

    return
