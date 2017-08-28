import random_data as rd


# Unit conversions
_ms_per_s = 1000.


class ComplexCorrelationFunction(rd.array.SpatialCrossCorrelation):
    '''A class for computing the complex-valued correlation
    function corresponding to the PCI measurements.

    Attributes:
    -----------
    shot - int
        The DIII-D shot number.

    Gxy - array_like, (`L`, `Nf`)
        The ensemble-averaged complex correlation as a function of:

            - detector-element separation (1st index, `L`), and
            - frequency (2nd index, `Nf`),

        where ensemble averaging has been done

            - in time (i.e. application of ergodic theorem), and
            - at each measurement *separation*.

        The indexing in `L` is such that correlation estimates are
        ordered sequentially from most-negative element separation
        to most-positive element.

        [Gxy] = [Ph_pci.x]^2 / [Fs], where `Ph_pci` and `Fs` are
            provided at initialization

    separation - array_like, (`L`,)
        The element separation.
        [separation] = [Ph_pci.detector_elements], where `Ph_pci`
            is provided at initialization

    equalize - bool
        If True, the powers in the signals provided at initialization
        were equalized prior to performing any spectral calculations.

    The additional attributes:

        {`detrend`, `df`, `dt`, `f`, `Fs`, `Npts_overlap`,
        `Npts_per_ens`, `Npts_per_real`, `Nreal_per_ens`, `t`}

    are described in the documentation for :py:class:`CrossSpectralDensity
    <random_data.spectra.CrossSpectralDensity>`.

    Methods:
    --------
    Type `help(ComplexCorrelationFunction)` in the IPython console
    for a listing.

    '''
    def __init__(self, Ph_pci, tlim=None, equalize=True,
                 print_status=True, **csd_kwargs):
        '''Create an instance of the `ComplexCorrelationFunction` class.

        Input parameters:
        -----------------
        Ph_pci - :py:class:`Phase <mitpci.pci.Phase>` instance

        tlim - array_like, (2,) or None
            If not None, then only use the portion of `Ph_pci`
            that sits between `min(tlim)` and `max(tlim)`.
            [tlim] = s

        equalize - bool
            If True, scale the power from each detector channel
            to be equal to that of of the channel with the most power.
            Note that this may e.g. increase the effective bit noise
            in channels with initially low power.

        print_status - bool
            If True, print status of computations.

        csd_kwargs - any valid keyword arguments for
            :py:class:`CrossSpectralDensity
                <random_data.spectra.CrossSpectralDensity>`.

            The signal sample rate `Fs` and initial timestamp `t0`
            must be specified; spectral-estimation parameters can
            also be specified. For example, use

                corr = ComplexCorrelationFunction(
                    Ph_pci, ..., Nreal_per_ens=100)

            to indicate that `Nreal_per_ens` realizations should be
            averaged over to obtain the ensemble spectral estimate.

            Note that additional parameters of relevance to spectral
            estimation (such as windowing, window overlap, etc.) are
            specified via the keyword packing `**csd_kwargs`. See the
            `CrossSpectralDensity` documentation for further details.

        '''
        self.shot = Ph_pci.shot

        csd_kwargs['Fs'] = Ph_pci.Fs
        csd_kwargs['t0'] = Ph_pci.t0

        rd.array.SpatialCrossCorrelation.__init__(
            self, Ph_pci.x, Ph_pci.detector_elements,
            tlim=tlim, equalize=equalize,
            print_status=print_status, **csd_kwargs)

    def plotNormalizedCorrelationFunction(
            self, xlim=None, flim=None, vlim=[-1, 1],
            units='lab',
            cmap='viridis', interpolation='none', fontsize=16,
            xlabel=None, ylabel=None, no_nan=False):
        'Plot normalized correlation function, Gxy(delta, f) / Gxy(0, f).'
        # Determine y-axis units
        if str.lower(units) == 'lab':
            fmult = 1. / _ms_per_s
            funits = 'kHz'
        else:
            fmult = 1.
            funits = 'Hz'

        # Determine labels, if needed
        if xlabel is None:
            xsymbol = r'$\mathregular{\delta}$'
            xlabel = 'element separation, %s' % xsymbol

        if ylabel is None:
            ysymbol = 'f'
            ylabel = (r'$\mathregular{%s \; [%s]}$'
                      % (ysymbol, funits))

        # Convert frequency to desired units for plotting
        self.f *= fmult

        # Plot
        rd.array.SpatialCrossCorrelation.plotNormalizedCorrelationFunction(
            self, xlim=xlim, flim=flim, vlim=vlim,
            cmap=cmap, interpolation=interpolation, fontsize=fontsize,
            xlabel=xlabel, ylabel=ylabel, no_nan=no_nan)

        # Convert frequency back to original units
        self.f /= fmult

        return
