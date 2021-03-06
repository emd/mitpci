import numpy as np
import MDSplus as mds
import random_data as rd


# Unit conversions
_cm_per_m = 100.
_ms_per_s = 1000.


class TwoDimensionalAutoSpectralDensity(
        rd.spectra2d.TwoDimensionalAutoSpectralDensity):
    '''A class for estimating the PCI's 2-dimensional autospectral density
    given the corresponding complex-valued correlation function.

    Attributes:
    -----------
    Sxx - array_like, (`Nk`, `Nf`)
        An array of the estimated autospectral density as a
        function of:

            - wavenumber (1st index, `Nk`), and
            - frequency (2nd index, `Nf`).

        `Sxx` is normalized such that integrating over all of `Sxx`
        yields the total power in the signal PCI signal.

        [Sxx] = rad^2 / (Hz * (rad / m))

    k - array_like, (`Nk`,)
        The wavenumber grid.
        [k] = rad / m

    f - array_like, (`Nf`,)
        The frequency grid.
        [f] = Hz

    dk - float
        The spacing of the wavenumber grid.
        [dk] = rad / m

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
            burg_params={'p': 5, 'Nk': 1000},
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
            are: {'p', 'Nk'} where

            - p: int, order of Burg AR spectral-density estimate,
            - Nk: int, the number of points in the wavenumber grid;
                note that the spatial estimate is *two-sided*.

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

        if str.lower(spatial_method) == 'burg':
            burg_params['Nxi'] = burg_params['Nk']

        # Compute spectrum
        rd.spectra2d.TwoDimensionalAutoSpectralDensity.__init__(
            self, corr, spatial_method=spatial_method,
            burg_params=burg_params,
            fourier_params=fourier_params)

        # Get maximum wavenumber
        self._kmax = kmax(self.shot)

        if self._kmax < 0:
            # kmax < 0 indicates that PCI detector element #1
            # maps to the innermost major radius of the beam.
            # To make contact with the usual convention that
            # k_R > 0 indicates propagation in the positive
            # major-radial direction, we need to flip `Sxx`.
            self.Sxx = np.flipud(self.Sxx)

            # Because `self.xi` is *not exactly* symmetric
            # about zero, we also need to flip and negate it
            self.xi = -self.xi[::-1]

        # Create wavenumber grid
        xi2k_conversion = np.abs(self._kmax) / np.max(np.abs(self.xi))
        self.k =  self.xi * xi2k_conversion
        self.dk = self.k[1] - self.k[0]

        del self.xi, self.dxi

    def plotSpectralDensity(self, klim=None, flim=None, vlim=None,
                            units='lab',
                            cmap='viridis', interpolation='none',
                            fontsize=16, title=None,
                            xlabel=None, ylabel=None,
                            cblabel=None, cborientation='vertical',
                            ax=None, fig=None, geometry=111):
        'Plot magnitude of spectral density on log scale.'
        if str.lower(units) == 'lab':
            kmult = 1. / _cm_per_m
            kunits = 'cm^{-1}'

            fmult = 1. / _ms_per_s
            funits = 'kHz'

            Sxxmult = _cm_per_m * _ms_per_s
            Sxxunits = 'rad^2 / (kHz \cdot cm^{-1})'
        else:
            kmult = 1.
            kunits = 'm^{-1}'

            fmult = 1.
            funits = 'Hz'

            Sxxmult = 1
            Sxxunits = 'rad^2 / (Hz \cdot m^{-1})'

        if xlabel is None:
            xsymbol = 'k'
            xlabel = (r'$\mathregular{%s \; [%s]}$'
                      % (xsymbol, kunits))

        if ylabel is None:
            ysymbol = 'f'
            ylabel = (r'$\mathregular{%s \; [%s]}$'
                      % (ysymbol, funits))

        if cblabel is None:
            cbsymbol = '|S_{xx}(k,f)|'
            cblabel = (r'$\mathregular{%s \; [%s]}$'
                       % (cbsymbol, Sxxunits))

        ax = rd.spectra.nonparametric._plot_image(
            kmult * self.k, fmult * self.f, Sxxmult * np.abs(self.Sxx.T),
            xlim=klim, ylim=flim, vlim=vlim,
            norm='log', cmap=cmap, interpolation=interpolation,
            title=title, xlabel=xlabel, ylabel=ylabel,
            cblabel=cblabel, cborientation=cborientation,
            fontsize=fontsize,
            ax=ax, fig=fig, geometry=geometry)

        return ax


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
