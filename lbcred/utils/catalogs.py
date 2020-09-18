import os
import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude, match_coordinates_sky, SkyCoord
from astropy.table import Table, vstack
from ..log import logger
from .. import package_dir, utils


__all__ = ['match_sky_coordinates','sextractor_cat_to_ds9reg']


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


def sextractor_cat_to_ds9reg(cat, color='green', tag='all', winparams=False,
                             outfile='sextractor.reg', drawmode='ellipse', textparam=None):
    """
    Write a ds9 region file from SExtractor output.
    Parameters
    ----------
    cat : structured ndarray
        Output cat from a sextractor run.
    color : string, optional
      Region color (cyan, blue, magenta, red, green,
      yellow, white, or black)
    tag : string, optional
      ds9 tag for all the regions
    winparams : bool, optional
        If True, use sextractor's windowed parameters.
    outfile : string, optional
        Output reg file name.
    drawmode : string, optional
        Draw an 'ellipse' or 'point' for every object
    textparam : string, optional
        If not None, write this sextractor output parameter
        next to each object in the catalog.
    Notes
    -----
     i) Adapted from https://github.com/nhmc/Barak.git.
    ii) The sextractor output file must contain X_IMAGE and
        Y_IMAGE for drawmode=point, and for drawmode=ellipse,
        it must also contain A_IMAGE, B_IMAGE, and THETA_IMAGE.
        The corresponding 'WIN' parameters are acceptable with
        winparams set to True.
    """
    assert (drawmode=='ellipse') or (drawmode=='point')
    regions = ['global font="helvetica 10 normal" select=1 highlite=1 '
               'edit=0 move=1 delete=1 include=1 fixed=0 source']
    regions.append('image')
    fields = ['X_IMAGE', 'Y_IMAGE', 'A_IMAGE','B_IMAGE','THETA_IMAGE']
    if drawmode=='point':
        fields = fields[:2]
    if winparams:
        fields = [f.split('_')[0]+'WIN'+'_'+f.split('_')[1] for f in fields]
    if textparam is not None:
        textfmt = 'text={%s}'
        fields.append(textparam)
    else:
        textfmt = ''
    fmt = {'ellipse':'ellipse(%s %s %s %s %s) # '+textfmt+' color=%s tag={%s}',
           'point':'point(%s %s) # point=circle '+textfmt+' color=%s tag={%s}'}[drawmode]
    for row  in cat[fields]:
        vals = list(row)
        vals.extend([color, tag])
        regions.append(fmt % tuple(vals))
    print('writing to region file to', outfile)
    fh = open(outfile,'w')
    fh.write('\n'.join(regions))
    fh.close()
