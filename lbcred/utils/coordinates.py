"""
Functions for calculating the time and location of things.
"""
import datetime
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import get_sun, get_moon, SkyCoord, EarthLocation
from astropy.coordinates import Longitude, Latitude, AltAz

from .. import utils
from ..log import logger
from .. import ResultStruct


__all__ = ['get_today', 'get_moon_status', 'get_image_corners',
           'has_wcs', 'nms', 'to_skycoord']



nms_lon = Longitude('-105:28:41', unit=u.deg)
nms_lat = Latitude('32:53:23', unit=u.deg)
nms_elevation = 7300 * u.imperial.ft
nms = EarthLocation.from_geodetic(nms_lon, nms_lat, nms_elevation)


def get_today():
    """
    Return today's date.
    """
    return datetime.date.today().strftime('%m/%d/%Y')


def get_moon_status(skycoord, utc_time):
    assert type(skycoord) == SkyCoord, 'Coords are not a SkyCoord object!'
    try:
        utc_time = Time(utc_time)
        moon = get_moon(utc_time, nms)

        nms_altaz = AltAz(obstime=utc_time, location=nms)
        moon_nms = moon.transform_to(nms_altaz)
        target_nms = skycoord.transform_to(nms_altaz)
        sep_from_target = moon.separation(skycoord)

        sun = get_sun(utc_time)
        elongation = sun.separation(moon)
        phase = np.arctan2(sun.distance * np.sin(elongation),
                           moon.distance - sun.distance * np.cos(elongation))
        illumination = (1 + np.cos(phase.value)) / 2.0

        results = ResultStruct(
            alt=moon_nms.alt,
            az=moon_nms.az,
            phase=phase,
            illumination=illumination,
            sep_from_target=sep_from_target,
            target_alt=target_nms.alt,
            target_az=target_nms.az,
            success=True
        )
    except Exception as e:
        logger.warning('Moon data calculation failed: {}'.format(e))
        results = ResultStruct(success=False)
    return results


def get_image_corners(path_or_header):
    """
    Get the sky coordinates of the corners of an image taking SIP
    corrections into account.
    """
    header = utils.load_path_or_header(path_or_header)
    wcs = WCS(header)
    corners = wcs.calc_footprint()
    return corners


def has_wcs(path_or_header):
    """
    Check if header has a WCS.
    """
    header = utils.load_path_or_header(path_or_header)
    return 'WCSAXES' in header.keys()


def to_skycoord(coord):
    if type(coord) == SkyCoord:
        sc = coord
    elif type(coord[0]) == str:
        assert type(coord[1]) == str, 'ra & dec must be in the same format!'
        ra, dec = coord
        ra = ra.split(':')
        ra = f'{ra[0]}h{ra[1]}m{ra[2]}s'
        dec = dec.replace(':', 'm') + 's'
        sc =  SkyCoord(ra, dec)
    else:
        ra, dec = coord
        sc =  SkyCoord(ra, dec, unit='deg')
    return sc
