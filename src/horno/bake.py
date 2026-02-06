import math
import os.path

import numpy as np

import horno.fits
import horno.image
import horno.instrument
import horno.path

_darkdata = None
_flatdata = None


def readdark(exposuretime, path="dark-{exposuretime:.0f}.fits", name="readdark"):
    global _darkdata
    path = path.format(exposuretime=exposuretime)
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _darkdata = horno.fits.readproductdata(path)
    else:
        raise RuntimeError("no dark found.")
    return _darkdata


def readflat(path="flat.fits", name="readflat"):
    global _flatdata
    path = path.format()
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _flatdata = horno.fits.readproductdata(path)
    else:
        raise RuntimeError("no flat found.")
    return _flatdata


def writedark(path="dark-{exposuretime:.0f}.fits", exposuretime=None, name="writebias"):
    path = path.format(exposuretime=exposuretime)
    print("%s: writing %s." % (name, path))
    horno.fits.writeproduct(path, _darkdata, exposuretime=exposuretime)
    return


def writeflat(path="flat.fits", name="writeflat"):
    path = path.format()
    print("%s: writing %s." % (name, path))
    horno.fits.writeproduct(path, _flatdata)
    return


def bake(
    fitspath,
    name="bake",
    dotrim=False,
    dodark=False,
    doflat=False,
    dosky=False,
    dowindow=False,
    dorotate=False,
    nwindow=None,
    nmargin=0,
):

    print("%s: reading %s." % (name, os.path.basename(fitspath)))
    header, data = horno.fits.readraw(fitspath)

    # Set invalid pixels to nan.
    data[np.where(data == horno.instrument.datamax(header))] = np.nan

    if (
        dotrim
        and horno.instrument.trimyslice(header) is not None
        and horno.instrument.trimxslice(header) is not None
    ):
        print("%s: trimming." % (name))
        data = data[
            horno.instrument.trimyslice(header), horno.instrument.trimxslice(header)
        ]

    if dodark and _darkdata is not None:
        print("%s: subtracting dark." % (name))
        data -= _darkdata

    if doflat and _flatdata is not None:
        print("%s: dividing by flat." % (name))
        data /= _flatdata

    if dosky:
        median = np.nanmedian(data)
        print("%s: subtracting median sky of %.1f DN." % (name, median))
        # data -= np.nanmedian(data, axis=0, keepdims=True)
        # data -= np.nanmedian(data, axis=1, keepdims=True)
        data -= np.nanmedian(data, keepdims=True)

    if dorotate:
        print("%s: rotating to standard orientation." % (name))
        data = horno.instrument.dorotate(header, data)

    if nwindow is not None:

        print("%s: windowing to %d by %d." % (name, nwindow, nwindow))

        assert nwindow <= data.shape[0]
        assert nwindow <= data.shape[1]
        ylo = int((data.shape[0] - nwindow) / 2)
        yhi = ylo + nwindow
        xlo = int((data.shape[1] - nwindow) / 2)
        xhi = xlo + nwindow
        data = data[ylo:yhi, xlo:xhi].copy()

    return header, data


def usefakebias():
    global _biasdata
    _biasdata = None
    return _biasdata


def usefakedark():
    global _darkdata
    _darkdata = None
    return _darkdata


def usefakeflat():
    global _flatdata
    _flatdata = None
    return _flatdata


def makedark(
    fitspaths, exposuretime, darkpath="dark-{exposuretime}.fits", fitspathsslice=None
):

    print("makedark: making %.0f second dark from %s." % (exposuretime, fitspaths))

    fitspathlist = horno.path.getrawfitspaths(
        fitspaths, exposuretime=exposuretime, fitspathsslice=fitspathsslice
    )

    if len(fitspathlist) == 0:
        print("ERROR: no dark files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = bake(fitspath, name="makedark", dotrim=True)
        headerlist.append(header)
        datalist.append(data)

    print("makedark: averaging %d darks with rejection." % len(datalist))
    global _darkdata
    _darkdata, darksigma = horno.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    mean, sigma = horno.image.clippedmeanandsigma(_darkdata, sigma=5)
    print("makedark: dark is %.2f ± %.2f DN." % (mean, sigma))

    sigma = horno.image.clippedmean(darksigma, sigma=5) / math.sqrt(len(datalist))
    print("makedark: estimated noise in dark is %.2f DN." % sigma)

    horno.image.show(_darkdata, zscale=True)

    writedark(darkpath, exposuretime=exposuretime, name="makedark")

    print("makedark: finished.")

    return


def makeflat(fitspaths, flatpath="flat.fits", fitspathsslice=None):

    ############################################################################

    print("makeflat: making flat %s." % fitspaths)

    ############################################################################

    print("makeflat: making flat without mask.")

    fitspathlist = horno.path.getrawfitspaths(fitspaths, fitspathsslice=fitspathsslice)

    if len(fitspathlist) == 0:
        print("ERROR: no flat files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = bake(
            fitspath,
            name="makeflat",
            dotrim=True,
            dodark=True,
        )
        centeryslice = slice(int(data.shape[0] * 1 / 4), int(data.shape[0] * 3 / 4))
        centerxslice = slice(int(data.shape[1] * 1 / 4), int(data.shape[1] * 3 / 4))
        if np.isnan(data[centeryslice, centerxslice]).all():
            print(
                "makedark: rejected %s: no valid data in center."
                % os.path.basename(fitspath)
            )
            continue
        median = np.nanmedian(data[centeryslice, centerxslice])
        print("makeflat: median in center is %.2f DN." % median)
        if median > horno.instrument.flatmax(header):
            print("makeflat: rejecting image: median in center is too high.")
            continue
        print("makeflat: accepted %s." % os.path.basename(fitspath))

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

        print(
            "makeflat: normalizing 00, 01, 10, and 11 pixels by %.1f, %.1f, %.1f, and %.1f."
            % (median00, median01, median10, median11)
        )

        meanmedian = 0.25 * (median00 + median01 + median10 + median11)
        q = (median00 - median11) / meanmedian
        u = (median01 - median10) / meanmedian
        print(
            "makeflat: apparent polarization in flat is q = %+.3f u = %+.3f." % (q, u)
        )

        data[0::2, 0::2] /= median00
        data[0::2, 1::2] /= median01
        data[1::2, 0::2] /= median10
        data[1::2, 1::2] /= median11

        headerlist.append(header)
        datalist.append(data)

    print("makeflat: averaging %d flats with rejection." % (len(datalist)))

    flatdata, flatsigma = horno.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    ############################################################################

    print("makeflat: making mask.")

    maskdata = np.ones(flatdata.shape, dtype="float32")

    print("makeflat: masking nan values.")
    maskdata[np.isnan(flatdata)] = 0

    print("makeflat: masking inf values.")
    maskdata[np.isinf(flatdata)] = 0

    print("makeflat: masking globally low pixels.")
    maskdata[np.where(flatdata < 0.80)] = 0

    print("makeflat: masking locally high or low pixels.")
    low = horno.image.medianfilter(flatdata, 7)
    high = flatdata / low
    maskdata[np.where(high < 0.9)] = 0
    maskdata[np.where(high > 1.1)] = 0

    print("makeflat: masking pixels with at least two masked neighbors.")
    # Grow the mask so that any pixel with at least 2 neigboring bad pixels is also bad.
    grow = horno.image.uniformfilter(maskdata, size=3)
    maskdata[np.where(grow <= 7 / 9)] = 0

    print("makeflat: fraction of masked pixels is %.5f." % (1 - np.nanmean(maskdata)))
    centeryslice = slice(int(maskdata.shape[0] * 1 / 4), int(maskdata.shape[0] * 3 / 4))
    centerxslice = slice(int(maskdata.shape[1] * 1 / 4), int(maskdata.shape[1] * 3 / 4))
    print(
        "makeflat: fraction of masked pixels in center is %.5f."
        % (1 - np.nanmean(maskdata[centeryslice, centerxslice]))
    )

    horno.image.show(maskdata, zrange=True)

    ############################################################################

    print("makeflat: making flat with mask.")

    maskeddatalist = []
    for data in datalist:
        data[np.where(maskdata == 0)] = np.nan
        data / np.nanmedian(data)
        maskeddatalist.append(data)

    print("makeflat: averaging %d flats with rejection." % (len(maskeddatalist)))
    flatdata, flatsigma = horno.image.clippedmeanandsigma(
        maskeddatalist, sigma=3, axis=0
    )

    mean, sigma = horno.image.clippedmeanandsigma(flatdata, sigma=5)
    print("makeflat: flat is %.2f ± %.3f." % (mean, sigma))

    sigma = horno.image.clippedmean(flatsigma, sigma=5) / math.sqrt(len(maskeddatalist))
    print("makeflat: estimated noise in flat is %.4f." % sigma)

    global _flatdata
    _flatdata = flatdata
    horno.image.show(_flatdata, zrange=True)
    writeflat(flatpath, name="makeflat")

    ############################################################################

    print("makeflat: finished.")

    return



def makeobjects(fitspaths, fitspathsslice=None):

    ############################################################################

    print("makeobjects: making objects %s." % fitspaths)

    ############################################################################

    fitspathlist = horno.path.getrawfitspaths(fitspaths, fitspathsslice=fitspathsslice)

    if len(fitspathlist) == 0:
        print("ERROR: no object files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = bake(
            fitspath,
            name="makeobjects",
            dotrim=True,
            dodark=True,
            doflat=True,
        )

        headerlist.append(header)
        datalist.append(data)

    ############################################################################

    print("makeobjects: finished.")

    return headerlist, datalist
