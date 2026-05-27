import math

import numpy as np

import photutils.detection
import photutils.aperture
import photutils.psf

import horno.image
import horno.log

def pointsource(
    header,
    data,
    y=None,
    x=None,
    yslice=None,
    xslice=None,
    finderthreshold=10,
    show=False,
    verbose=False,
    name="pointsource",
    logger=None
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
        fwhm=15, threshold=finderthreshold * sigma, brightest=1
    )
    findertable = finder(datasum)
    xcenter = findertable["xcentroid"].value[0]
    ycenter = findertable["ycentroid"].value[0]
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

    if not math.isnan(fwhm):

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
            skynpix = skystatistics.sum_aper_area.value

            objectsum = objectsum - objectnpix * skylevel

            gain = 2
            objectstd = np.sqrt(objectsum * gain) / gain

            return objectsum[0]

        dcenter = 0.25
        s00 = photometerone(
            data00,
            ycenter + dcenter,
            xcenter + dcenter,
        )
        s01 = photometerone(
            data01,
            ycenter + dcenter,
            xcenter - dcenter,
        )
        s10 = photometerone(
            data10,
            ycenter - dcenter,
            xcenter + dcenter,
        )
        s11 = photometerone(
            data11,
            ycenter - dcenter,
            xcenter - dcenter,
        )

        smean = 0.25 * (s00 + s01 + s10 + s11)
        q = (s00 - s11) / smean
        u = (s01 - s10) / smean
        if verbose:
            logger(name, "smean = %.3f q = %+.3f u = %+.3f" % (smean, q, u))
