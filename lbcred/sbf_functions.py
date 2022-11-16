import numpy as np
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy import optimize
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft
from collections import namedtuple
from astropy import units as u
from photutils import EllipticalAperture
from astropy.io import fits
import pymfit, os, sep


SBFResults = namedtuple(
    'SBFResults', 'ps_image ps_image_err ps_psf npix k p cov fit_func'
)


def _measure_power_spectrum(data):
    fft = np.fft.fftshift(np.fft.fft2(data))
    ps = np.abs(fft)**2
    return ps


def _get_wavenumbers(window_length, num_radial_bins=45):
    k = np.fft.fftshift(np.fft.fftfreq(window_length))
    kx, ky = np.meshgrid(k, k)
    k = np.sqrt(kx**2 + ky**2)
    k, _ = azimuthal_average(k, num_radial_bins)

    return k


def azimuthal_average(image, num_radial_bins=50):
    """
    Compute radial profile of image.

    Parameters
    ----------
    image : ndarray
        Input image.
    num_radial_bins : int
        Number of radial bins in profile.

    Returns
    -------
    radial_mean : ndarray
        Mean intensity within each annulus.
    radial_err : ndarray
        Standard error on the mean: sigma / sqrt(N).
    """
    ny, nx = image.shape
    yy, xx = np.mgrid[:ny, :nx]
    center = np.array(image.shape) / 2

    r = np.hypot(xx - center[1], yy - center[0])
    rbin = (num_radial_bins * r/r.max()).astype(np.int)

    radial_mean = ndimage.mean(
        image, labels=rbin, index=np.arange(1, rbin.max() + 1))

    radial_stddev = ndimage.standard_deviation(
        image, labels=rbin, index=np.arange(1, rbin.max() + 1))

    npix = ndimage.sum(np.ones_like(image), labels=rbin,
                       index=np.arange(1, rbin.max() + 1))

    radial_err = radial_stddev / np.sqrt(npix)
    return radial_mean, radial_err


sbf_xlabel = r'Wave Number $\left[\mathrm{pixel}^{-1}\right]$'
sbf_ylabel = r'Power'
def sbf_results(sbf, residual_image, subplots=None, xlabel=sbf_xlabel,
                ylabel=sbf_ylabel, xscale='linear', percentiles=[0.1, 99.9],
                yscale='log', plot_errors=True, ylim_factors=[0.5, 1.1],
                cmap='gray_r',save_fn=None, normalize_ps = False, plot_blank_fields=False, blank_results = {}): 

    if subplots is None:
        fig, ax = plt.subplots(1, 2, figsize=(15, 6.5))
        fig.subplots_adjust(wspace=0.2)
    else:
        fig, ax = subplots

    vmin, vmax = np.nanpercentile(residual_image, percentiles)

    ax[0].imshow(residual_image, origin='lower', cmap=cmap,
                 vmin=vmin, vmax=vmax)

    ax[0].set(xticks=[], yticks=[])
    
    if normalize_ps:
        normalization = sbf.p[1] / sbf.npix

    else:
        normalization = 1.
        
    ax[1].axhline(y=sbf.p[1] / normalization / sbf.npix, ls='--', c='gray', lw=3)   # White noise floor
    norm = sbf.fit_func(sbf.k, *sbf.p / normalization).max() / sbf.ps_psf.max()
    ax[1].plot(sbf.k, sbf.ps_psf * norm / sbf.npix, ls='--',
                   c='gray', lw=3)                                  # PSF contribution to power spectrum

    ax[1].plot(sbf.k, sbf.fit_func(sbf.k, *sbf.p / normalization) / sbf.npix,
                   c='slateblue', lw=3)                             # Power spectrum fit
    if plot_errors:
        ax[1].errorbar(sbf.k, sbf.ps_image / normalization / sbf.npix,
                       yerr=sbf.ps_image_err / sbf.npix,
                       fmt='k.', lw=2, capsize=3, capthick=1.5, zorder=10)  # Measured power spectrum - points with errorbars
    else:
        ax[1].plot(sbf.k, sbf.ps_image / normalization / sbf.npix, 'k-', lw=2)              # Measured power spectrum - line
        
    if plot_blank_fields:
        for key in blank_results:
            ax[1].plot(blank_results[key]['results'].k, blank_results[key]['results'].ps_image / normalization / blank_results[key]['results'].npix, color='tan', lw=2, zorder=0, alpha=0.5)

    ymax = (((sbf.ps_image / normalization) + sbf.ps_image_err) / sbf.npix).max()

    ax[1].set_xscale(xscale)
    ax[1].set_yscale(yscale)
    ax[1].set_ylim(sbf.p[1] / normalization * ylim_factors[0] / sbf.npix,
                   ymax * ylim_factors[1])
    ax[1].set_xlabel(xlabel, fontsize=25)
    ax[1].set_ylabel(ylabel, fontsize=25)
    ax[1].tick_params(labelsize=20)
    ax[1].minorticks_on()
    if save_fn is not None:
        fig.savefig(save_fn, bbox_inches='tight', dpi=500)

    return fig, ax


def measure_sbf(normed_res_image, psf, mask=None, k_range=[0.01, 0.4],
                fit_param_guess=[100, 50], num_radial_bins=45,
                use_sigma=False, **kwargs):

    res_image = normed_res_image.copy()
    '''
    if res_image.shape[0] > psf.shape[0]:
        psf_padded = np.pad(
            psf, (res_image.shape[0] - psf.shape[0])//2, 'constant')
        psf_padded /= psf_padded.sum()
    else:
        shape = psf.shape[0]
        if shape%2 == 0:
            x =
            shape-=1
        psf_padded = Cutout2D(psf, (x,y), (shape,shape)).data
    '''

    psf_padded = np.pad(
        psf, (res_image.shape[0] - psf.shape[0])//2, 'constant')
    psf_padded /= psf_padded.sum()
    ps_psf = _measure_power_spectrum(psf_padded)

    npix = np.product(res_image.shape)
    if mask is not None:
        res_image[mask.astype(bool)] = 0.0
        mask = (~mask.astype(bool)).astype(float)
        npix = mask.sum()
        ps_mask = _measure_power_spectrum(mask)
        ps_psf = convolve_fft(ps_psf, ps_mask, boundary='fill',
                              normalize_kernel=True)

    ps_image = _measure_power_spectrum(res_image)
    ps_image, ps_image_err = azimuthal_average(ps_image, num_radial_bins)
    ps_psf, _ = azimuthal_average(ps_psf, num_radial_bins)
    wavenumbers = _get_wavenumbers(res_image.shape[0], num_radial_bins)

    # apply cut on wave number
    k_cut = (wavenumbers >= k_range[0]) & (wavenumbers <= k_range[1])
    ps_image = ps_image[k_cut]
    ps_psf = ps_psf[k_cut]
    ps_image_err = ps_image_err[k_cut]
    wavenumbers = wavenumbers[k_cut]

    # define fitting function: psf(k)*p0 + p1
    psf_k = interp1d(wavenumbers, ps_psf)
    fit_func = lambda k, p0, p1: psf_k(k) * p0 + p1

    # perform fit
    sigma = ps_image_err if use_sigma else None
    popt, pcov = optimize.curve_fit(
        fit_func, wavenumbers, ps_image, p0=fit_param_guess,
        sigma=sigma, **kwargs)

    # consolidate results
    results = SBFResults(ps_image=ps_image,
                         ps_image_err=ps_image_err,
                         ps_psf=ps_psf,
                         npix=npix,
                         k=wavenumbers,
                         p=popt,
                         cov=pcov,
                         fit_func=fit_func)

    return results


def get_sbf_distance(results, zeropoint, color, gain, texp, colorterm = None, extinction_correction=None, blank_field_results=None):

    if blank_field_results != None:
        sbf_mag = zeropoint - 2.5*np.log10((results.p[0]-blank_field_results.p[0])*gain/texp/results.npix)

    else:
        sbf_mag = zeropoint - 2.5*np.log10(results.p[0]*gain/texp/results.npix)

    if colorterm is not None:
        sbf_mag -= colorterm*color

    if extinction_correction is not None:
        sbf_mag -= extinction_correction

    M_r_a = (6.09 * color) - 8.81            # Eqn. 2 from Jerjen 2000
    M_r_b = 1.89*((color-0.77)**2) - 1.26    # Eqn. 3 from Jerjen 2000

    dist_mod_a = sbf_mag - M_r_a
    dist_mod_b = sbf_mag - M_r_b

    d_a = 10**((dist_mod_a + 5)/5) # parsecs
    d_b = 10**((dist_mod_b + 5)/5) # parsecs

    return sbf_mag, d_a, d_b

def elliptical_mask(shape, pars, scale=2):
    ell = EllipticalAperture(
        [pars['X0'], pars['Y0']], scale * pars['r_e'],
        scale * pars['r_e'] * (1 - pars['ell']),
        np.deg2rad(pars['PA'] - 90)
    )
    ell_mask = ~ell.to_mask().to_image(shape).astype(bool)
    return ell_mask

def get_sbf_mask_resid(model_fn, resid_fn, sersic_params, grow_obj, scale, config, blank_field=False):

    resid = fits.open(resid_fn)
    model = fits.open(model_fn)
    color = config['color_name']
    resid_no_norm = resid[config['ext']].data
    resid[config['ext']].data = resid[config['ext']].data/np.sqrt(model[config['ext']].data)
    # Make aper mask
    # Ellipticity = 1 - b/a, b = semiminor axis, a = semimajor axis
    aper_mask = elliptical_mask(resid[config['ext']].data.shape, sersic_params, scale=scale)

    sbf_resid_fn = model_fn.replace(f'.fits',f'_sbf_resid.fits')
    resid.writeto(sbf_resid_fn,overwrite=config['overwrite'])

    # Make pymfit masks
    if blank_field:
        mask_fn = resid_fn.replace(f'.fits','_sbf_mask.fits')
    else:
        mask_fn = model_fn.replace(f'.fits','_sbf_mask.fits')
    dist_mod = 5*np.log10(config['assumed_distance']*10**6)-5
    app_mag = dist_mod+config['given_abs_mag_to_mask']
    resid_no_norm = np.ascontiguousarray(resid_no_norm)
    obj, seg, bkg, img = pymfit.masking.detect_sources(resid_no_norm, config['masking_sbf']['thresh'], config['masking_sbf']['backsize'], mask=aper_mask, return_all=True, kern_sig=config['masking_sbf']['kern_sig'])
    flux=obj['flux']

    obj = obj[(-2.5*np.log10(flux*config['gain']/config['exposure_time'])+config['zpt']) < app_mag]
    seg_mask = pymfit.masking.make_seg_mask(seg, grow_obj, config['masking_sbf']['thresh'])
    obj_mask = pymfit.masking.make_obj_mask(obj, resid[config['ext']].data.shape, grow_obj)
    mask = (seg_mask | obj_mask).astype(int)
    fits.writeto(mask_fn, mask, overwrite=config['overwrite'])

    # Combine aperture and pymfit masks, save masks to file
    sbf_mask = ((aper_mask.astype(bool) | mask.astype(bool))).astype(float)
    mask_file = fits.open(mask_fn)
    mask_file[0].data = sbf_mask
    mask_file.writeto(mask_fn,overwrite=True)

    vmin, vmax = np.nanpercentile(resid[config['ext']].data, [5e-1, 92.])
    plt.clf()
    plt.imshow(resid[config['ext']].data,vmin=vmin,vmax=vmax,cmap='gray_r')
    plt.imshow(sbf_mask,cmap='BuPu',alpha=0.3)
    plt.title(f'SFB residual : {color}-band')
    plt.savefig(os.path.join(config['out_dir'], f'sbf_resid_mask_{color}.png'))

    return resid, sbf_mask
