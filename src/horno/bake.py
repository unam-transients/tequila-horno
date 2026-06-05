import math
import os.path

import numpy as np

import horno.fits
import horno.image
import horno.instrument
import horno.path
import horno.log

################################################################################

_darkdata = {}
_flatdata = None


def _darkpath(exposuretime, detectortemperature, tag):
    if tag is None:
        tagtext = ""
    else:
        tagtext = "-%s" % tag
    return "dark-%d@%+d%s.fits" % (exposuretime, detectortemperature, tagtext)


def _flatpath(tag):
    if tag is None:
        tagtext = ""
    else:
        tagtext = "-%s" % tag
    return "flat%s.fits" % (tagtext)


def _setdarkdata(exposuretime, detectortemperature, data):
    exposuretime = int(exposuretime)
    detectortemperature = int(detectortemperature)
    key = "%d@%+d" % (exposuretime, detectortemperature)
    global _darkdata
    _darkdata[key] = data


def _getdarkdata(exposuretime, detectortemperature):
    exposuretime = int(exposuretime)
    detectortemperature = int(detectortemperature)
    key = "%d@%+d" % (exposuretime, detectortemperature)
    if key not in _darkdata:
        return np.nan
    else:
        return _darkdata[key]


def _setflatdata(data):
    global _flatdata
    _flatdata = data


def _getflatdata():
    if _flatdata is None:
        raise RuntimeError("no flat is available.")
    return _flatdata


################################################################################


def readdark(
    exposuretime,
    detectortemperature,
    tag=None,
    name="readdark",
    logger=None,
    verbose=True,
):
    logger = horno.log.getlogger(logger)
    path = _darkpath(exposuretime, detectortemperature, tag)
    if not os.path.exists(path):
        raise RuntimeError("%s not found." % path)
    if verbose:
        logger(name, "reading %s." % (path))
    data = horno.fits.readproductdata(path, name=name, logger=logger, verbose=verbose)
    _setdarkdata(exposuretime, detectortemperature, data)
    return data


def readflat(tag=None, name="readflat", logger=None, verbose=True):
    logger = horno.log.getlogger(logger)
    path = _flatpath(tag)
    if not os.path.exists(path):
        raise RuntimeError("%s not found." % path)
    if verbose:
        logger(name, "reading %s." % (path))
    data = horno.fits.readproductdata(path, name=name, logger=logger, verbose=verbose)
    _setflatdata(data)
    return data


def writedark(
    exposuretime,
    detectortemperature,
    tag=None,
    name="writedark",
    logger=None,
    verbose=True,
):
    logger = horno.log.getlogger(logger)
    path = _darkpath(exposuretime, detectortemperature, tag)
    data = _getdarkdata(exposuretime, detectortemperature)
    if verbose:
        logger(name, "writing %s." % (path))
    horno.fits.writeproduct(path, data, name=name, logger=logger, verbose=verbose)
    return


def writeflat(
    tag=None,
    name="writeflat",
    logger=None,
    verbose=True,
):
    logger = horno.log.getlogger(logger)
    path = _flatpath(tag)
    data = _getflatdata()
    if verbose:
        logger(name, "writing %s." % (path))
    horno.fits.writeproduct(path, data, name=name, logger=logger, verbose=verbose)
    return


################################################################################


def usefakedark(
    exposuretime, detectortemperature, name="usefakedark", logger=None, verbose=True
):
    logger = horno.log.getlogger(logger)
    if verbose:
        logger(name, "using fake dark.")
    data = 0
    _setdarkdata(exposuretime, detectortemperature, data)
    return data


def usefakeflat(name="usefakeflat", logger=None, verbose=True):
    logger = horno.log.getlogger(logger)
    if verbose:
        logger(name, "using fake flat.")
    data = 1
    _setflatdata(data)
    return data


################################################################################


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
    header, data = horno.fits.readraw(
        fitspath, name=name, logger=logger, verbose=verbose
    )

    # Set invalid pixels to nan.
    data[np.where(data == horno.instrument.datamax(header))] = np.nan

    if (
        dooffset
        and horno.instrument.offsetyslice(header) is not None
        and horno.instrument.offsetxslice(header) is not None
    ):
        if verbose:
            logger(name, "removing offset.")
        offset = np.nanmedian(
            data[
                horno.instrument.offsetyslice(header),
                horno.instrument.offsetxslice(header),
            ]
        )
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

    if dodark:
        if verbose:
            logger(name, "subtracting dark.")
        data -= _getdarkdata(
            horno.instrument.exposuretime(header),
            horno.instrument.detectortemperature(header),
        )

    if doflat:
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


def makedark(
    fitspaths,
    exposuretime,
    detectortemperature=None,
    fitspathsslice=None,
    tag=None,
    name="makedark",
    logger=None,
    verbose=True,
    show=True,
):

    logger = horno.log.getlogger(logger)

    if verbose:
        logger(
            name,
            "making %.0f second dark at %+.0f C from %s."
            % (exposuretime, detectortemperature, fitspaths),
        )

    fitspathlist = horno.path.getrawfitspaths(
        fitspaths,
        exposuretime=exposuretime,
        detectortemperature=detectortemperature,
        fitspathsslice=fitspathsslice,
        name=name,
        logger=logger,
        verbose=False,
    )

    if len(fitspathlist) == 0:
        raise RuntimeError("no dark files found.")
        return

    if verbose:
        logger(name, "%d dark files found." % len(fitspathlist))

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        fitsbasename = os.path.basename(fitspath)
        header = horno.fits.readrawheader(
            fitspath, name="makedark", logger=logger, verbose=verbose
        )
        # Reject darks taken in daylight, civil twilight, or evening nautical
        # twilight, as these have noticeable light leaks.
        startskystate = header["SSNSK"]
        endskystate = header["ESNSK"]
        startsunha = float(header["SSNHA"])
        if startskystate == "daylight" or endskystate == "daylight":
            if verbose:
                logger(name, "rejecting %s as it was taken in daylight." % fitsbasename)
            continue
        if startskystate == "civiltwilight" or endskystate == "civiltwilight":
            if verbose:
                logger(
                    name,
                    "rejecting %s as it was taken in civil twilight." % fitsbasename,
                )
            continue
        if (
            startskystate == "nauticaltwilight" or endskystate == "nauticaltwilight"
        ) and startsunha > 0:
            if verbose:
                logger(
                    name,
                    "rejecting %s as it was taken in evening nautical twilight."
                    % fitsbasename,
                )
            continue
        if verbose:
            logger(
                name,
                "accepting %s." % fitsbasename,
            )
        header, data = bake(
            fitspath,
            dooffset=True,
            dotrim=True,
            name="makedark",
            logger=logger,
            verbose=verbose,
        )
        headerlist.append(header)
        datalist.append(data)

    if len(datalist) == 0:
        raise RuntimeError("no files accepted.")

    if verbose:
        logger(name, "averaging %d darks with rejection." % len(datalist))
    datamean, datasigma = horno.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    mean, sigma = horno.image.clippedmeanandsigma(datamean, sigma=5)
    if verbose:
        logger(name, "dark is %.2f ± %.2f DN." % (mean, sigma))

    sigma = horno.image.clippedmean(datasigma, sigma=5) / math.sqrt(len(datalist))
    if verbose:
        logger(name, "estimated noise in dark is %.2f DN." % sigma)

    if show:
        horno.image.show(data, zscale=True)

    _setdarkdata(exposuretime, detectortemperature, data)
    writedark(
        exposuretime,
        detectortemperature,
        tag,
        name="makedark",
        logger=logger,
        verbose=verbose,
    )

    if verbose:
        logger(name, "finished.")

    return data


def makeflat(
    fitspaths,
    fitspathsslice=None,
    tag=None,
    pmax=0.2,
    name="makeflat",
    logger=None,
    verbose=True,
):

    logger = horno.log.getlogger(logger)

    ############################################################################

    if verbose:
        logger(name, "making flat %s." % fitspaths)

    ############################################################################

    logger(name, "making flat without mask.")

    fitspathlist = horno.path.getrawfitspaths(
        fitspaths, fitspathsslice=fitspathsslice, name=name, logger=logger
    )

    if len(fitspathlist) == 0:
        raise RuntimeError("no flat files found.")
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
            if verbose:
                logger(
                    name,
                    "rejected %s: no valid data in center."
                    % os.path.basename(fitspath),
                )
                continue
        median = np.nanmedian(data[centeryslice, centerxslice])
        if verbose:
            logger(name, "median in center is %.2f DN." % median)
        if median > horno.instrument.flatmax(header):
            if verbose:
                logger(name, "rejecting image: median in center is too high.")
            continue
        if verbose:
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

        if verbose:
            logger(
                name,
                "normalizing 00, 01, 10, and 11 pixels by %.1f, %.1f, %.1f, and %.1f."
                % (median00, median01, median10, median11),
            )

        q = (median00 - median11) / (median00 + median11)
        u = (median01 - median10) / (median01 + median10)
        p = np.sqrt(np.square(q) + np.square(u))
        if verbose:
            logger(
                name,
                "apparent polarization in flat is q = %+.3f u = %+.3f p = %.3f."
                % (q, u, p),
            )
        if p > pmax:
            if verbose:
                logger(name, "rejecting image: apparent polarization is too high.")
            continue

        data[0::2, 0::2] /= median00
        data[0::2, 1::2] /= median01
        data[1::2, 0::2] /= median10
        data[1::2, 1::2] /= median11

        headerlist.append(header)
        datalist.append(data)

    if verbose:
        logger(name, "averaging %d flats with rejection." % (len(datalist)))

    flatdata, flatsigma = horno.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    ############################################################################

    if verbose:
        logger(name, "making mask.")

    maskdata = np.ones(flatdata.shape, dtype="float32")

    logger(name, "masking nan values.")
    maskdata[np.isnan(flatdata)] = 0

    if verbose:
        logger(name, "masking inf values.")
    maskdata[np.isinf(flatdata)] = 0

    if verbose:
        logger(name, "masking globally low pixels.")
    maskdata[np.where(flatdata < 0.80)] = 0

    if verbose:
        logger(name, "masking locally high or low pixels.")
    low = horno.image.medianfilter(flatdata, 7)
    high = flatdata / low
    maskdata[np.where(high < 0.9)] = 0
    maskdata[np.where(high > 1.1)] = 0

    if verbose:
        logger(name, "masking pixels with at least two masked neighbors.")
    # Grow the mask so that any pixel with at least 2 neigboring bad pixels is also bad.
    grow = horno.image.uniformfilter(maskdata, size=3)
    maskdata[np.where(grow <= 7 / 9)] = 0

    if verbose:
        logger(name, "fraction of masked pixels is %.5f." % (1 - np.nanmean(maskdata)))
    centeryslice = slice(int(maskdata.shape[0] * 1 / 4), int(maskdata.shape[0] * 3 / 4))
    centerxslice = slice(int(maskdata.shape[1] * 1 / 4), int(maskdata.shape[1] * 3 / 4))
    if verbose:
        logger(
            name,
            "fraction of masked pixels in center is %.5f."
            % (1 - np.nanmean(maskdata[centeryslice, centerxslice])),
        )

    horno.image.show(maskdata, zrange=True)

    ############################################################################

    if verbose:
        logger(name, "making flat with mask.")

    maskeddatalist = []
    for data in datalist:
        data[np.where(maskdata == 0)] = np.nan
        data /= np.nanmedian(data)
        maskeddatalist.append(data)

    if verbose:
        logger(name, "averaging %d flats with rejection." % (len(maskeddatalist)))
    flatdata, flatsigma = horno.image.clippedmeanandsigma(
        maskeddatalist, sigma=3, axis=0
    )

    mean, sigma = horno.image.clippedmeanandsigma(flatdata, sigma=5)
    if verbose:
        logger(name, "flat is %.2f ± %.3f." % (mean, sigma))

    sigma = horno.image.clippedmean(flatsigma, sigma=5) / math.sqrt(len(maskeddatalist))
    if verbose:
        logger(name, "estimated noise in flat is %.4f." % sigma)

    global _flatdata
    _flatdata = flatdata
    horno.image.show(_flatdata, zrange=True)
    if verbose:
        writeflat(tag=tag, name="makeflat", logger=logger, verbose=verbose)

    ############################################################################

    if verbose:
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
