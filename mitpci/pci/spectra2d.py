import MDSplus as mds
import random_data as rd


class TwoDimensionalAutoSpectralDensity(
        rd.spectra2d.TwoDimensionalAutoSpectralDensity):
    '''A class for estimating the PCI's 2-dimensional autospectral density
    given the corresponding complex-valued correlation function.

    Attributes:
    -----------
    Sxx - array_like, (`Nxi`, `Nf`)
        An array of the estimated autospectral density as a
        function of:

            - xi = 1 / wavelength (1st index, `Nxi`), and
            - frequency (2nd index, `Nf`).

        `Sxx` is normalized such that integrating over all of `Sxx`
        yields the total power in the signal PCI signal.

        [Sxx] = rad^2 / (Hz * m)

    xi - array_like, (`Nxi`,)
        The inverse-spatial grid. Note that xi = (1 / wavelength) such that
        the wavenumber k is related to xi via k = (2 * pi * xi).
        [xi] = 1 / m

    f - array_like, (`Nf`,)
        The frequency grid.
        [f] = [self.Fs]

    dxi - float
        The spacing of the inverse-spatial grid.
        [dxi] = 1 / m

    df - float
        The spacing of the frequency grid.
        [df] = Hz

    Fs - float
        The temporal sampling rate.
        [Fs] = samples / second

    Fs_spatial - float
        The "spatial sampling rate", as determined by the spacing between
        adjacent points in the correlation function `corr`, provided
        at input.
        [Fs_spatial] = 1 / [corr.separation]

    spatial_method - string
        The method used to estimate the spatial content of the
        autospectral density, `self.Gxx`.

    p - int
        The order of the Burg AR. Only present if `spatial_method == 'burg'`.

    spatial_window - :py:func, a function of an integer
        The tapering window applied to the spatial dimension of the
        correlation function `corr` prior to taking the spatial Fourier
        transform. Only present (and only applied) if
        `spatial_method == 'fourier'`.

    The additional attributes:

        {`Npts_overlap`, `Npts_per_ens`, `Npts_per_real`,
        `Nreal_per_ens`, `detrend`, `dt`, `t`, `window`}

    are described in the documentation for :py:class:`SpatialCrossCorrelation
    <random_data.array.SpatialCrossCorrelation>`.

    Methods:
    --------
    Type `help(TwoDimensionalAutoSpectralDensity)` in the IPython console
    for a listing.

    '''
    def __init__(
            self, corr, spatial_method='burg',
            burg_params=rd.spectra2d.default_burg_params,
            fourier_params=rd.spectra2d.default_fourier_params):
        '''Create instance of `TwoDimensionalAutoSpectralDensity` class.

        Input parameters:
        -----------------
        corr - :py:class:`ComplexCorrelationFunction
                <mitpci.pci.ComplexCorrelationFunction>` instance
            A PCI `ComplexCorrelationFunction` instance.

        spatial_method - string
            The method to use when estimating the spatial spectral
            density. Valid methods are:

                - 'fourier': use FFT-based estimation, or
                - 'burg': use Burg AR estimation.

            Specification of other values will raise a ValueError.

        burg_params - dict
            A dictionary containing the parameters of relevance
            for Burg spectral estimation. Valid dictionary keys
            are: {'p', 'Nxi'} where

            - p: int, order of Burg AR spectral-density estimate,
            - Nxi: int, number of points in the two-sided *spatial*
                spectral-density estimate at each frequency.

            See documentation for :py:class:`BurgAutoSpectralDensity
            <random_data.spectra.parametric.BurgAutoSpectralDensity>`
            for more information on these parameters.

        fourier_params - dict
            A dictionary containing the parameters of relevance
            for Fourier spectral estimation. Valid dictionary keys
            are: {'window'} where

            - window: tapering function to be applied to spatial
                dimension of correlation function prior to calculating
                the FFT; should be a function of the window length,
                such as `np.hanning`. If `None`, do not apply a window
                to the spatial dimension of the correlation function.
                Although the correlation function typically smoothly
                tapers to zero on its own, the application of a
                window can still suppress leakage, at the cost of
                marginally decreased resolution.

        '''
        self.shot = corr.shot

        rd.spectra2d.TwoDimensionalAutoSpectralDensity.__init__(
            self, corr, spatial_method='burg',
            burg_params=rd.spectra2d.default_burg_params,
            fourier_params=rd.spectra2d.default_fourier_params)


def kmax(shot):
    '''Get PCI's maximum wavenumber for DIII-D shot number `shot`.

    Input parameters:
    -----------------
    shot - int
        The DIII-D shot number.

    Returns:
    --------
    kmax - float
        The PCI's maximum wavenumber for `shot`. Note that
        kmax > 0 indicates that PCI detector element #1 maps
        to the outermost major radius.
        [kmax] = rad / m

    '''
    tree = mds.Tree('pci', shot=shot, mode='ReadOnly')
    node = tree.getNode('.OPTICS.KMAX')
    kmax = node.getData().data()

    return kmax
