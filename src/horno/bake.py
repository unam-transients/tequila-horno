import math
import os.path

import numpy as np

import horno.fits
import horno.image
import horno.instrument
import horno.path
import horno.log

_darkdata = None
_flatdata = None


def readdark(
    exposuretime, path="dark-{exposuretime:.0f}.fits", name="readdark", logger=None
):
    logger = horno.log.getlogger(logger)
    global _darkdata
    path = path.format(exposuretime=exposuretime)
    if os.path.exists(path):
        logger(name, "reading %s." % (path))
        _darkdata = horno.fits.readproductdata(path, name=name, logger=logger)
    else:
        raise RuntimeError("no dark found.")
    return _darkdata


def readflat(path="flat.fits", name="readflat", logger=None):
    logger = horno.log.getlogger(logger)
    global _flatdata
    path = path.format()
    if os.path.exists(path):
        logger(name, "reading %s." % (path))
        _flatdata = horno.fits.readproductdata(path, name=name, logger=logger)
    else:
        raise RuntimeError("no flat found.")
    return _flatdata


def writedark(
    path="dark-{exposuretime:.0f}.fits",
    exposuretime=None,
    name="writebias",
    logger=None,
):
    logger = horno.log.getlogger(logger)
    path = path.format(exposuretime=exposuretime)
    logger(name, "writing %s." % (path))
    horno.fits.writeproduct(
        path, _darkdata, exposuretime=exposuretime, name=name, logger=logger
    )
    return


def writeflat(path="flat.fits", name="writeflat", logger=None):
    logger = horno.log.getlogger(logger)
    path = path.format()
    logger(name, "writing %s." % (path))
    horno.fits.writeproduct(path, _flatdata, name=name, logger=logger)
    return


def bake(
    fitspath,
    name="bake",
    dooffset=False,
    dotrim=False,
    dodark=False,
    doflat=False,
    dosky=False,
    dowindow=False,
    dorotate=False,
    nwindow=None,
    nmargin=0,
    logger=None,
    verbose=True,
):

    logger = horno.log.getlogger(logger)

    if verbose:
        logger(name, "reading %s." % (os.path.basename(fitspath)))
    header, data = horno.fits.readraw(fitspath, name=name, logger=logger, verbose=verbose)

    # Set invalid pixels to nan.
    data[np.where(data == horno.instrument.datamax(header))] = np.nan

    if (
        dooffset
        and horno.instrument.offsetyslice(header) is not None
        and horno.instrument.offsetxslice(header) is not None
    ):
        if verbose:
            logger(name, "removing offset.")
        offset = np.nanmedian(data[horno.instrument.offsetyslice(header), horno.instrument.offsetxslice(header)])
        if verbose:
            logger(name, "offset is %.2f." % offset)
        data -= offset

    if (
        dotrim
        and horno.instrument.trimyslice(header) is not None
        and horno.instrument.trimxslice(header) is not None
    ):
        if verbose:
            logger(name, "trimming.")
        data = data[
            horno.instrument.trimyslice(header), horno.instrument.trimxslice(header)
        ]

    if dodark and _darkdata is not None:
        if verbose:
            logger(name, "subtracting dark.")
        data -= _darkdata

    if doflat and _flatdata is not None:
        if verbose:
            logger(name, "dividing by flat.")
        data /= _flatdata

    if dosky:
        median = np.nanmedian(data)
        if verbose:
            logger(name, "subtracting median sky of %.1f DN." % (median))
        # data -= np.nanmedian(data, axis=0, keepdims=True)
        # data -= np.nanmedian(data, axis=1, keepdims=True)
        data -= np.nanmedian(data, keepdims=True)

    if dorotate:
        if verbose:
            logger(name, "rotating to standard orientation.")
        data = horno.instrument.dorotate(header, data, name=name, logger=logger)

    if nwindow is not None:

        if verbose:
            logger(name, "windowing to %d by %d." % (nwindow, nwindow))

        assert nwindow <= data.shape[0]
        assert nwindow <= data.shape[1]
        ylo = int((data.shape[0] - nwindow) / 2)
        yhi = ylo + nwindow
        xlo = int((data.shape[1] - nwindow) / 2)
        xhi = xlo + nwindow
        data = data[ylo:yhi, xlo:xhi].copy()

    return header, data


def usefakebias(name=None, logger=None):
    global _biasdata
    _biasdata = None
    return _biasdata


def usefakedark(name=None, logger=None):
    logger = horno.log.getlogger(logger)
    global _darkdata
    _darkdata = None
    return _darkdata


def usefakeflat(name=None, logger=None):
    logger = horno.log.getlogger(logger)
    global _flatdata
    _flatdata = None
    return _flatdata


def makedark(
    fitspaths,
    exposuretime,
    detectortemperature=None,
    darkpath="dark-{exposuretime}.fits",
    fitspathsslice=None,
    name="makedark",
    logger=None,
):

    logger = horno.log.getlogger(logger)

    logger(name, "making %.0f second dark from %s." % (exposuretime, fitspaths))

    fitspathlist = horno.path.getrawfitspaths(
        fitspaths,
        exposuretime=exposuretime,
        detectortemperature=detectortemperature,
        fitspathsslice=fitspathsslice,
        name=name,
        logger=logger,
    )

    if len(fitspathlist) == 0:
        logger(name, "ERROR: no dark files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        fitsbasename = os.path.basename(fitspath)
        header = horno.fits.readrawheader(fitspath, name="makedark", logger=logger)
        # Reject darks taken in daylight, civil twilight, or evening nautical
        # twilight, as these have noticeable light leaks.
        startskystate = header["SSNSK"]
        endskystate = header["ESNSK"]
        startsunha = float(header["SSNHA"])
        if startskystate == "daylight" or endskystate == "daylight":
            logger(name, "rejecting %s as it was taken in daylight." % fitsbasename)
            continue
        if startskystate == "civiltwilight" or endskystate == "civiltwilight":
            logger(
                name, "rejecting %s as it was taken in civil twilight." % fitsbasename
            )
            continue
        if (
            startskystate == "nauticaltwilight" or endskystate == "nauticaltwilight"
        ) and startsunha > 0:
            logger(
                name,
                "rejecting %s as it was taken in evening nautical twilight."
                % fitsbasename,
            )
            continue
        logger(
            name,
            "accepting %s." % fitsbasename,
        )
        header, data = bake(fitspath, dooffset=True, dotrim=True, name="makedark", logger=logger)
        headerlist.append(header)
        datalist.append(data)

    logger(name, "averaging %d darks with rejection." % len(datalist))
    global _darkdata
    _darkdata, darksigma = horno.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    mean, sigma = horno.image.clippedmeanandsigma(_darkdata, sigma=5)
    logger(name, "dark is %.2f ± %.2f DN." % (mean, sigma))

    sigma = horno.image.clippedmean(darksigma, sigma=5) / math.sqrt(len(datalist))
    logger(name, "estimated noise in dark is %.2f DN." % sigma)

    horno.image.show(_darkdata, zscale=True)

    writedark(darkpath, exposuretime=exposuretime, name="makedark")

    logger(name, "finished.")

    return


def makeflat(
    fitspaths,
    flatpath="flat.fits",
    fitspathsslice=None,
    pmax=0.2,
    name="makeflat",
    logger=None,
):

    logger = horno.log.getlogger(logger)

    ############################################################################

    logger(name, "making flat %s." % fitspaths)

    ############################################################################

    logger(name, "making flat without mask.")

    fitspathlist = horno.path.getrawfitspaths(
        fitspaths, fitspathsslice=fitspathsslice, name=name, logger=logger
    )

    if len(fitspathlist) == 0:
        logger(name, "ERROR: no flat files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = bake(
            fitspath,
            name="makeflat",
            dooffset=True,
            dotrim=True,
            dodark=True,
        )
        centeryslice = slice(int(data.shape[0] * 1 / 4), int(data.shape[0] * 3 / 4))
        centerxslice = slice(int(data.shape[1] * 1 / 4), int(data.shape[1] * 3 / 4))
        if np.isnan(data[centeryslice, centerxslice]).all():
            logger(
                name,
                "rejected %s: no valid data in center." % os.path.basename(fitspath),
            )
            continue
        median = np.nanmedian(data[centeryslice, centerxslice])
        logger(name, "median in center is %.2f DN." % median)
        if median > horno.instrument.flatmax(header):
            logger(name, "rejecting image: median in center is too high.")
            continue
        logger(name, "accepted %s." % os.path.basename(fitspath))

        centeryslice = slice(
            int(data.shape[0] / 2 * 1 / 4), int(data.shape[0] / 2 * 3 / 4)
        )
        centerxslice = slice(
            int(data.shape[1] / 2 * 1 / 4), int(data.shape[1] / 2 * 3 / 4)
        )
        median00 = np.nanmedian(data[0::2, 0::2][centeryslice, centerxslice])
        median01 = np.nanmedian(data[0::2, 1::2][centeryslice, centerxslice])
        median10 = np.nanmedian(data[1::2, 0::2][centeryslice, centerxslice])
        median11 = np.nanmedian(data[1::2, 1::2][centeryslice, centerxslice])

        logger(
            name,
            "normalizing 00, 01, 10, and 11 pixels by %.1f, %.1f, %.1f, and %.1f."
            % (median00, median01, median10, median11),
        )

        q = (median00 - median11) / (median00 + median11)
        u = (median01 - median10) / (median01 + median10)
        p = np.sqrt(np.square(q) + np.square(u))
        logger(
            name,
            "apparent polarization in flat is q = %+.3f u = %+.3f p = %.3f."
            % (q, u, p),
        )
        if p > pmax:
            logger(name, "rejecting image: apparent polarization is too high.")
            continue

        data[0::2, 0::2] /= median00
        data[0::2, 1::2] /= median01
        data[1::2, 0::2] /= median10
        data[1::2, 1::2] /= median11

        headerlist.append(header)
        datalist.append(data)

    logger(name, "averaging %d flats with rejection." % (len(datalist)))

    flatdata, flatsigma = horno.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    ############################################################################

    logger(name, "making mask.")

    maskdata = np.ones(flatdata.shape, dtype="float32")

    logger(name, "masking nan values.")
    maskdata[np.isnan(flatdata)] = 0

    logger(name, "masking inf values.")
    maskdata[np.isinf(flatdata)] = 0

    logger(name, "masking globally low pixels.")
    maskdata[np.where(flatdata < 0.80)] = 0

    logger(name, "masking locally high or low pixels.")
    low = horno.image.medianfilter(flatdata, 7)
    high = flatdata / low
    maskdata[np.where(high < 0.9)] = 0
    maskdata[np.where(high > 1.1)] = 0

    logger(name, "masking pixels with at least two masked neighbors.")
    # Grow the mask so that any pixel with at least 2 neigboring bad pixels is also bad.
    grow = horno.image.uniformfilter(maskdata, size=3)
    maskdata[np.where(grow <= 7 / 9)] = 0

    logger(name, "fraction of masked pixels is %.5f." % (1 - np.nanmean(maskdata)))
    centeryslice = slice(int(maskdata.shape[0] * 1 / 4), int(maskdata.shape[0] * 3 / 4))
    centerxslice = slice(int(maskdata.shape[1] * 1 / 4), int(maskdata.shape[1] * 3 / 4))
    logger(
        name,
        "fraction of masked pixels in center is %.5f."
        % (1 - np.nanmean(maskdata[centeryslice, centerxslice])),
    )

    horno.image.show(maskdata, zrange=True)

    ############################################################################

    logger(name, "making flat with mask.")

    maskeddatalist = []
    for data in datalist:
        data[np.where(maskdata == 0)] = np.nan
        data /= np.nanmedian(data)
        maskeddatalist.append(data)

    logger(name, "averaging %d flats with rejection." % (len(maskeddatalist)))
    flatdata, flatsigma = horno.image.clippedmeanandsigma(
        maskeddatalist, sigma=3, axis=0
    )

    mean, sigma = horno.image.clippedmeanandsigma(flatdata, sigma=5)
    logger(name, "flat is %.2f ± %.3f." % (mean, sigma))

    sigma = horno.image.clippedmean(flatsigma, sigma=5) / math.sqrt(len(maskeddatalist))
    logger(name, "estimated noise in flat is %.4f." % sigma)

    global _flatdata
    _flatdata = flatdata
    horno.image.show(_flatdata, zrange=True)
    writeflat(flatpath, name="makeflat")

    ############################################################################

    logger(name, "finished.")

    return


def makeobjects(
    fitspaths,
    fitspathsslice=None,
    dooffset=True,
    dotrim=True,
    dodark=True,
    doflat=True,
    name="makeobjects",
    logger=None,
    verbose=True,
):

    logger = horno.log.getlogger(logger)

    ############################################################################

    if verbose:
        logger(name, "making objects %s." % fitspaths)

    ############################################################################

    fitspathlist = horno.path.getrawfitspaths(
        fitspaths, fitspathsslice=fitspathsslice, name=name, logger=logger
    )

    if len(fitspathlist) == 0:
        logger(name, "ERROR: no object files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = bake(
            fitspath,
            name="makeobjects",
            dooffset=dooffset,
            dotrim=dotrim,
            dodark=dodark,
            doflat=doflat,
            verbose=verbose,
        )

        headerlist.append(header)
        datalist.append(data)

    ############################################################################

    if verbose:
        logger(name, "finished.")

    return headerlist, datalist
