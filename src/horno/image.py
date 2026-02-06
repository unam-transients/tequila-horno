import math
import warnings

import numpy as np

import astropy.stats
import astropy.visualization

import scipy.ndimage

import matplotlib.pyplot as plt

def sigma_clipped_stats(data, sigma=3.0, axis=None):
    """
    Return sigma-clipped statistics of the given data.

    This behaves exactly as:
    
        astropy.stats.sigma_clipped_stats(
            data, sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )

    except that it converts all ndarrays to float32 before returning them.

    Furthermore, for the common case of clipping a stack of 2D arrays, it does
    so row by row, which is much more efficient in terms of memory use.

    :param data: The data of which to calculate the statistics.
    :param sigma: The number of standard deviations for the upper and lower
        clipping limits. Defaults to 3.0
    :param axis: The axis along with to clip the data. Defaults to None.
    :return: The mean, median, and standard deviation of the data.
    """

    if not isinstance(data, np.ndarray):
        data = np.array(data)

    with warnings.catch_warnings():

        warnings.simplefilter("ignore", Warning)

        if axis == 0 and len(data.shape) == 3:

            ny = data.shape[1]
            nx = data.shape[2]

            meanimage = np.full([ny, nx], np.nan, dtype="float32")
            medianimage = np.full([ny, nx], np.nan, dtype="float32")
            sigmaimage = np.full([ny, nx], np.nan, dtype="float32")

            for iy in range(ny):
                meanrow, medianrow, sigmarow = astropy.stats.sigma_clipped_stats(
                        data[:, iy, :],
                        sigma=sigma,
                        axis=0,
                        cenfunc="median",
                        stdfunc="mad_std",
                    )
                
                meanimage[iy, :] = meanrow
                medianimage[iy, :] = medianrow
                sigmaimage[iy, :] = sigmarow
                                
            mean = meanimage
            median = medianimage
            sigma = sigmaimage

        else:

            mean, median, sigma = astropy.stats.sigma_clipped_stats(
                data, sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
            )

    if isinstance(mean, np.ndarray):
        mean = mean.astype("float32")
    if isinstance(median, np.ndarray):
        median = median.astype("float32")
    if isinstance(sigma, np.ndarray):
        sigma = sigma.astype("float32")

    return mean, median, sigma


def clippedmean(data, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = sigma_clipped_stats(data, sigma=sigma, axis=axis)
    return mean


def clippedsigma(data, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = sigma_clipped_stats(data, sigma=sigma, axis=axis)
    return sigma


def clippedmeanandsigma(data, sigma=3.0, axis=None):
    mean, median, sigma = sigma_clipped_stats(data, sigma=sigma, axis=axis)
    return mean, sigma


def medianfilter(data, size):
    return scipy.ndimage.median_filter(data, size)


def uniformfilter(data, size):
    return scipy.ndimage.filters.uniform_filter(data, size=size, mode="nearest")


def show(
    data, zrange=False, zscale=False, contrast=0.25, zmin=None, zmax=None, small=False
):

    if zmin is not None and zmax is not None:
        interval = astropy.visualization.ManualInterval(zmin, zmax)
    elif zrange:
        interval = astropy.visualization.MinMaxInterval()
    else:
        interval = astropy.visualization.ZScaleInterval(contrast=contrast)
    stretch = astropy.visualization.LinearStretch()
    norm = astropy.visualization.ImageNormalize(
        data, interval=interval, stretch=stretch
    )

    ny = data.shape[0]
    nx = data.shape[1]
    nmax = max(ny, nx)
    
    if np.max(data.shape) > 1000:
        tickinterval = 100
    else:
        tickinterval = int(math.pow(2, int(math.log2(nmax / 16))))
    ticks = list(np.linspace(-tickinterval * (nmax // 2 // tickinterval), +tickinterval * (nmax // 2 // tickinterval), 1 + 2 * (nmax // 2 // tickinterval)))


    if small:
        plt.figure(figsize=(5, 5))
    else:
        plt.figure(figsize=(10, 10))
    plt.imshow(data, origin="lower", norm=norm, extent=[-nx/2, +nx/2, -ny/2, +ny/2])
    plt.xticks(ticks, rotation=90)
    plt.yticks(ticks)
    plt.colorbar(fraction=0.046, pad=0.035)
    plt.show()
