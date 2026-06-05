import glob
import horno.fits
import horno.instrument


def getrawfitspaths(
    fitspaths,
    exposuretime=None,
    detectortemperature=None,
    fitspathsslice=None,
    name=None,
    logger=None,
    verbose=True,
):
    """
    Return an expanded and filtered list of FITS paths.

    The ``fitspaths`` argument will be expanded by :func:`glob.glob`. This
    expansion must give a list of names of FITS files or compressed FITS files.
    If ``exposuretime`` is not ``None``, then files which do not have that
    exposure time are eliminated from the list. Finally, if ``fitspathsslice``
    is not ``None``, then this list is sliced.

    :param fitspaths: A pattern to be expanded by :func:`glob.glob`. The
        expansion must give a list of names of FITS files or compressed FITS
        files.

    :param exposuretime: An exposure time as a number or ``None``. If it is not
        ``None``, then files which do not have that exposure time are eliminated
        from the expanded list.

    :param fitspathsslice: A slice, either of the string ``"firsthalf"`` or
        ``"secondhalf"``, or ``None``. If it is not ``None``, then the expanded
        and filtered list is sliced before being returned. The special values
        ``"firsthalf"`` or ``"secondhalf"`` refer, as might be expected, to
        slices containing the first half and second half of the file names.

    :return: An expanded, filtered, and sliced list of FITS file name.
    """

    def flatten(lists):
        """
        Return a flattened list of lists.
        """
        return [item for sublist in lists for item in sublist]

    if isinstance(fitspaths, str):
        fitspaths = [fitspaths]
    fitspaths = sorted(flatten(glob.glob(fitspath) for fitspath in fitspaths))

    if exposuretime is not None:
        fitspaths = list(
            fitspath
            for fitspath in fitspaths
            if horno.instrument.exposuretime(
                horno.fits.readrawheader(fitspath, name=name, logger=logger, verbose=verbose)
            )
            == exposuretime
        )
    if detectortemperature is not None:
        fitspaths = list(
            fitspath
            for fitspath in fitspaths
            if horno.instrument.detectortemperature(
                horno.fits.readrawheader(fitspath, name=name, logger=logger, verbose=verbose)
            )
            == detectortemperature
        )
    if fitspathsslice == "firsthalf":
        fitspathsslice = slice(None, len(fitspaths) // 2)
    elif fitspathsslice == "secondhalf":
        fitspathsslice = slice(len(fitspaths) // 2, None)
    if fitspathsslice is not None:
        fitspaths = fitspaths[fitspathsslice]
    return fitspaths
