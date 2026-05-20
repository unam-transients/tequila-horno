import numpy as np
import astropy.io.fits
import os
from datetime import datetime

import horno.log


def _ihdu(fitspath):
    """
    Return the HDU index for the actual data in a FITS file. This is 0
    for an uncompressed FITS file and 1 for a fpack-compressed FITS
    file. The suffix of the fitspath is used to determine if the file is
    compressed or not; if the suffix is ".fz", the file is assumed to be
    compressed.
    """
    if fitspath[-3:] == ".fz":
        return 1
    else:
        return 0


def readraw(fitspath, name=None, logger=None):
    logger = horno.log.getlogger(logger)
    logger(name, "reading header and data from FITS file %s."
        % (os.path.basename(fitspath))
    )
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    header = hdu[ihdu].header
    data = np.array(hdu[ihdu].data, dtype=np.float32)
    hdu.close()
    return header, data


def readrawheader(fitspath, name=None, logger=None):
    logger = horno.log.getlogger(logger)
    logger(name, "reading header from raw FITS file %s." % (os.path.basename(fitspath))
    )
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    header = hdu[ihdu].header
    hdu.close()
    return header


def readrawdata(fitspath, name=None, logger=None):
    logger = horno.log.getlogger(logger)
    logger(name, "reading data from raw file %s." % (os.path.basename(fitspath)))
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    data = np.array(hdu[ihdu].data, dtype=np.float32)
    hdu.close()
    return data


def readproduct(fitspath, name=None, logger=None):
    logger = horno.log.getlogger(logger)
    logger(
        "reading header from product file %s." % (os.path.basename(fitspath))
    )
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    header = hdu[ihdu].header
    data = np.array(hdu[ihdu].data, dtype=np.float32)
    hdu.close()
    return header, data


def readproductheader(fitspath, name=None, logger=None):
    logger = horno.log.getlogger(logger)
    logger(name, 
        "reading header from product file %s." % (os.path.basename(fitspath))
    )
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    header = hdu[ihdu].header
    hdu.close()
    return header


def readproductdata(fitspath, name=None, logger=None):
    logger = horno.log.getlogger(logger)
    logger(name, 
        "reading data from product file %s." % (os.path.basename(fitspath))
    )
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    data = np.array(hdu[ihdu].data, dtype=np.float32)
    hdu.close()
    return data


def writeproduct(
    fitspath,
    data,
    filter=None,
    starttimestamp=None,
    endtimestamp=None,
    exposuretime=None,
    gain=None,
    name=None, logger=None,
):
    logger = horno.log.getlogger(logger)
    logger(name, "writing product file %s." % (os.path.basename(fitspath)))
    header = astropy.io.fits.Header()
    if filter is not None:
        header.append(("FILTER", filter))
    if starttimestamp is not None:
        header.append(
            (
                "DATE-OBS",
                datetime.utcfromtimestamp(starttimestamp).isoformat(
                    "T", "milliseconds"
                ),
            )
        )
    if endtimestamp is not None:
        header.append(
            (
                "DATE-END",
                datetime.utcfromtimestamp(endtimestamp).isoformat("T", "milliseconds"),
            )
        )
    if exposuretime is not None:
        header.append(("EXPTIME", exposuretime))
    if gain is not None:
        header.append(("GAIN", gain))
    astropy.io.fits.writeto(fitspath, data, header, overwrite=True)
    return
