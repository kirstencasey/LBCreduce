import os
import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude, match_coordinates_sky, SkyCoord
from astropy.table import Table, vstack
from ..log import logger
from .. import package_dir, utils


__all__ = ['match_sky_coordinates']


def match_sky_coordinates(cat, ref_cat, cat_cols='ra,dec',
                          ref_cols='ra,dec', sep_max=2.0):
    """
    Match two catalogs based on their sky coordinates.

    Parameters
    ----------
    cat : astropy.table.Table
        Catalog you want to match.
    ref_cat : astropy.table.Table
        Reference catalog to match 'cat' with.  
    cat_cols : list or str (optional)
        Column names of 'cat' corresponding to RA and Dec. Will assume 
        SExtractor names by default.
    ref_cols : list or str (optional)
        Column names of 'ref_cat' corresponding to RA and Dec.
    sep_max : float (optional)
        Maximum allowable star separation between image and
        reference catalogue (arcsec). Default is 2.0 arcsec.

    Returns
    -------
    cat_match : astropy.table.Table
        Matched image catalog.
    ref_match : astropy.table.Table
        Matched reference catalog.
    sep : astropy.units.Quantity
        Angular separation between the matched sources.

    Notes
    -----
    The returned catalogs will be the same length, with each row corresponing 
    to the (hopefully same) matched object.
    """
    cat_cols = utils.list_of_strings(cat_cols)
    ref_cols = utils.list_of_strings(ref_cols)
    cat_sc = SkyCoord(cat[cat_cols[0]], cat[cat_cols[1]], unit='deg')
    ref_sc = SkyCoord(ref_cat[ref_cols[0]], ref_cat[ref_cols[1]], unit='deg')
    idx, sep, _ = match_coordinates_sky(cat_sc, ref_sc)
    match = sep.arcsec < sep_max
    cat_match = cat[match].copy()
    ref_match = ref_cat[idx[match]].copy()
    sep = sep[match]
    return cat_match, ref_match, sep
