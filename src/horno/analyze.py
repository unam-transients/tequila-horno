import math

import numpy as np

import photutils.detection
import photutils.aperture
import photutils.psf

import horno.image
import horno.instrument
import horno.log


def pointsource(
    header,
    data,
    y=None,
    x=None,
    yslice=None,
    xslice=None,
    finderthreshold=10,
    fwhm=10,
    show=False,
    verbose=False,
    name="pointsource",
    logger=None,
):
    logger = horno.log.getlogger(logger)

    if yslice is None:
        yslice = slice(0, None)
    if xslice is None:
        xslice = slice(0, None)
    data = data[yslice, xslice]

    data00 = data[0::2, 0::2]
    data01 = data[0::2, 1::2]
    data10 = data[1::2, 0::2]
    data11 = data[1::2, 1::2]
    datasum = data00 + data01 + data10 + data11

    mean, median, sigma = horno.image.sigmaclippedstats(datasum)
    finder = photutils.detection.DAOStarFinder(
        fwhm=fwhm,
        threshold=finderthreshold * sigma,
        n_brightest=1,
        exclude_border=True,
    )
    findertable = finder(datasum)
    if findertable is None:
        logger(name, "no sources found")
        return 0, 0, 0, 0, 0, 0, 0, 0, 0
    xcenter = findertable["x_centroid"].value[0]
    ycenter = findertable["y_centroid"].value[0]
    if verbose:
        logger(name, "ycenter = %6.1f xcenter = %6.1f" % (ycenter, xcenter))

    if False:
        # Experimentally determine the centroid shifts.
        for _data in [data00, data01, data10, data11]:
            findertable = finder(_data)
            x = findertable["xcentroid"].value[0]
            y = findertable["ycentroid"].value[0]
            dx = x - xcenter
            dy = y - ycenter
            logger(name, "dy = %+.2f dx = %+.2f" % (dy, dx))

    position = [[xcenter, ycenter]]
    fwhm = photutils.psf.fit_fwhm(
        datasum - median, xypos=position, fit_shape=41, fwhm=15, mask=np.isnan(datasum)
    )
    fwhm = fwhm[0]
    if verbose:
        logger(name, "fwhm = %.1f pixel (%.2f arcsec)" % (fwhm, fwhm * 0.15))

    if show:

        n = 100

        xstart = int(xcenter + 0.5) - n // 2
        ystart = int(ycenter + 0.5) - n // 2

        if xstart < 0:
            xstart = 0
        elif xstart >= datasum.shape[1]:
            xstart = datasum.shape[1] - n

        if ystart < 0:
            ystart = 0
        elif ystart >= datasum.shape[0]:
            ystart = datasum.shape[0] - n

        positioninwindow = [[xcenter - xstart, ycenter - ystart]]
        datawindow = datasum[ystart : ystart + n, xstart : xstart + n]

        if math.isnan(fwhm):
            apertureradius = []
        else:
            apertureradius = [fwhm, 2 * fwhm, 3 * fwhm]

        horno.image.show(
            datawindow,
            small=True,
            aperturexy=positioninwindow,
            apertureradius=apertureradius,
            zrange=True,
        )

    if math.isnan(fwhm):
        return 0, 0, 0, 0, 0, 0, 0, 0, 0

    objectradius = 1.0 * fwhm
    skyinnerradius = 2.0 * fwhm
    skyouterradius = 3.0 * fwhm

    def photometerone(data, y, x):

        position = [[x, y]]

        objectaperture = photutils.aperture.CircularAperture(position, objectradius)
        skyaperture = photutils.aperture.CircularAnnulus(
            position, skyinnerradius, skyouterradius
        )

        objectstatistics = photutils.aperture.ApertureStats(data, objectaperture)
        skystatistics = photutils.aperture.ApertureStats(data, skyaperture)

        objectsum = objectstatistics.sum
        objectnpix = objectstatistics.sum_aper_area.value

        skylevel = skystatistics.mean
        skystd = skystatistics.std

        objectsum = objectsum - objectnpix * skylevel

        objectstd = np.sqrt(
            max(0, objectsum) * horno.instrument.gain(header)
        ) / horno.instrument.gain(header)

        totalstd = np.sqrt(np.square(objectstd) + objectnpix * np.square(skystd))

        return objectsum[0], totalstd[0]

    dcenter = 0.25
    n00, sn00 = photometerone(
        data00,
        ycenter + dcenter,
        xcenter + dcenter,
    )
    n01, sn01 = photometerone(
        data01,
        ycenter + dcenter,
        xcenter - dcenter,
    )
    n10, sn10 = photometerone(
        data10,
        ycenter - dcenter,
        xcenter + dcenter,
    )
    n11, sn11 = photometerone(
        data11,
        ycenter - dcenter,
        xcenter - dcenter,
    )

    if verbose:
        logger(name, "n00 = %.3e n00 / sn00 = %.1f" % (n00, n00 / sn00))
        logger(name, "n01 = %.3e n01 / sn01 = %.1f" % (n01, n01 / sn01))
        logger(name, "n10 = %.3e n10 / sn10 = %.1f" % (n10, n10 / sn10))
        logger(name, "n11 = %.3e n11 / sn11 = %.1f" % (n11, n11 / sn11))

    nmean = 0.25 * (n00 + n01 + n10 + n11)
    snmean = np.sqrt(
        0.25 * (np.square(sn00) + np.square(sn00) + np.square(sn10) + np.square(sn11))
    )

    q = (n00 - n11) / (n00 + n11)
    u = (n01 - n10) / (n01 + n10)

    sq = (
        2
        / np.square(n00 + n11)
        * np.sqrt(np.square(n00 * sn11) + np.square(n11 * sn00))
    )
    su = (
        2
        / np.square(n01 + n10)
        * np.sqrt(np.square(n01 * sn10) + np.square(n10 * sn01))
    )

    if verbose:
        logger(
            name,
            "nmean = %.2e ± %.2e  q = %+.4f ± %.4f  u = %+.4f ± %.4f"
            % (nmean, snmean, q, sq, u, su),
        )

    return ycenter, xcenter, fwhm, nmean, snmean, q, sq, u, su
