Python tools for retrieving and analyzing signals
digitized by the mitpci system.


Background:
===========
The mitpci system digitizes data from MIT's (a) phase contrast imaging (PCI)
and (b) heterodyne interferometer systems. The systems share a beam path
through the DIII-D vessel and a digitizer, but they use different
interference schemes and detectors to make their respective measurements.

As these systems share a digitization system, the mechanics underlying
signal retrieval from both systems are identical. This module's first
priority is to provide Python tools for robust and flexible signal retrieval.

The second priority of this module is to provide specialized routines
for analysis of the respective signals.
There are several specialized analysis routines for
the heterodyne interferometer, including:

* compensation of demodulator imperfections (i.e. DC offsets,
  amplitude imbalance, and phase imbalance between the I and Q signals),
  and
* correlation with the toroidally separated V2 interferometer
  for toroidal mode-number identification.

There are also several specialized analysis routines for the PCI, including:

* estimation of the complex-valued, spatial cross-correlation function and
* estimation of the corresponding two-dimensional autospectral density
  with Burg autoregression and Fourier spatial methods.

The use of these routines is discussed below.


Installation:
=============

For use on GA's Iris cluster:
-----------------------------
To use the `mitpci` package, change to the directory
you'd like to download the source files to, e.g.

    $ cd $HOME/python/projects

and retrieve the source files from github by typing

    $ git clone https://github.com/emd/mitpci.git

The created `mitpci` directory defines the
package's top-level directory. The `mitpci` package
depends on several additional packages:

* [random_data](https://github.com/emd/random_data),
* [filters](https://github.com/emd/filters),
* [bci](https://github.com/emd/bci),
* [magnetics](https://github.com/emd/magnetics), and
* [fit_ellipse](https://github.com/ndvanforeest/fit_ellipse).

These packages should similarly be cloned, e.g.

    $ git clone https://github.com/emd/random_data.git
    $ git clone https://github.com/emd/filters.git
    $ git clone https://github.com/emd/bci
    $ git clone https://github.com/emd/magnetics
    $ git clone https://github.com/ndvanforeest/fit_ellipse.git

Now, package management is cleanly handled on Iris via
[modules](https://diii-d.gat.com/diii-d/Iris#Environment_modules).
The above packages all have corresponding modulefiles
[here](https://github.com/emd/modulefiles).
Change to the directory you'd like to download the modulefiles to, e.g.

    $ cd $HOME

and retrieve the modulefiles by typing

    $ git clone https://github.com/emd/modulefiles

If your directory structure *differs* from that suggested above,
you will need to slightly modify the corresponding modulefiles.
Specifically, at the top of the `mitpci`
[modulefile](https://github.com/emd/modulefiles/blob/master/mitpci),
there is a TCL variable named `mitpci_root`;
this must be altered to point at the
top-level directory of the cloned `mitpci` package.
Similarly, the TCL variable `modulefiles_dir`
must be altered to point at the directory containing
the modulefiles for
`random_data`,
`filters`,
`bci`,
`magnetics`, and
`fit_ellipse`.
That's it! You shouldn't need to change anything else in
the `mitpci` modulefile.
(Of course, similar changes will also need to be made
to the `<package>_root` TCL variable in the modulefile
for the
`random_data`,
`filters`,
`bci`,
`magnetics`, and
`fit_ellipse`
modules).
The `mitpci` module can
then be loaded, unloaded, etc., as is discussed in the
above-linked Iris documentation.

The modulefile also defines a series of automated tests
for the `mitpci` package. Run these tests at the command line
by typing

    $ test_mitpci

If the tests return "OK", the installation should be working. Note that
the `mitpci` module and its associated tests require access to the
mitpci MDSplus server, which is located behind GA's firewall; if you
are not within the firewall, the module will err and the tests will fail
due to the inability to read from the mitpci MDSplus server.

For use on other systems:
-------------------------
First, install the `mitpci` dependencies that are *not* on PyPI:

* [random_data](https://github.com/emd/random_data),
* [filters](https://github.com/emd/filters),
* [bci](https://github.com/emd/bci),
* [magnetics](https://github.com/emd/magnetics), and
* [fit_ellipse](https://github.com/ndvanforeest/fit_ellipse).

(Installation instructions provided in each of the links).

Then, define an environmental variable `$pci_path` specifying
the appropriate MDSplus server's tree-path definitions
(`hermit.gat.com::/trees/pci`)
by, for example, adding the following to your `.bashrc`

    $ export pci_path='hermit.gat.com::/trees/pci'

(While data is digitized on `mitpci`, it should (ideally)
always be transferred to `hermit` prior to analysis;
digitization and writing is very resource intensive, and
an untimely request for data retrieval from `mitpci` could cause
sufficient loading to result in data loss or a system crash).

Now, change to the directory you'd like to download the source files to
and retrieve the source files from github by typing

    $ git clone https://github.com/emd/mitpci.git

Change into the `mitpci` top-level directory by typing

    $ cd mitpci

For accounts with root access, install by running

    $ python setup.py install

For accounts without root access (e.g. a standard account on GA's Venus
cluster), install locally by running

    $ python setup.py install --user

To test your installation, run

    $ nosetests tests/

If the tests return "OK", the installation should be working. Note that
the `mitpci` module and its associated tests require access to the
mitpci MDSplus server, which is located behind GA's firewall; if you
are not within the firewall, the module will err and the tests will fail
due to the inability to read from the mitpci MDSplus server.


Use:
====


Raw signals:
------------
To retrieve a signal digitized on channel `channel` on the
"mitpci" digitizer system from DIII-D shot `shot` between
times `tlim`, use the `mitpci.signal.Signal` class

```python
import mitpci

shot = 167342
channel = 8        # channel corresponding to PCI measurement
tlim = [1.0, 2.5]  # [tlim] = s

sig = mitpci.signal.Signal(shot, channel, tlim=tlim)

```

The digitized signal can then be accessed via `sig.x`
(raw signal, in bits), while the timebase can be generated via `sig.t()`.
The signal can be converted to volts via `sig.x * sig.volts_per_bit`.
Note that the `mitpci.signal.Signal` class allows retrieval of both
PCI and heterodyne interferometer measurements.


Interferometer I&Q:
-------------------
A specialized `mitpci.interferometer.Lissajous` class
also exists for retrieval of the heterodyne interferometer's
in-phase (I) and quadrature (Q) signals
(plotting Q vs. I produces a Lissajous figure), e.g.

```python
L = mitpci.interferometer.Lissajous(shot, tlim=tlim)

```

Note that the `mitpci.interferometer.Lissajous` class
automatically compensates for demodulator imperfections
(i.e. DC offsets, amplitude imbalance, and phase imbalance
between I and Q).
The post-processed I and Q objects are accessed via
`L.I` and `L.Q`, respectively.
Both `L.I` and `L.Q` are *instances* of the `mitpci.signal.Signal` class,
so the signal is accessed via e.g. `L.I.x`, and
the timebase can be generated via e.g. `L.I.t()`.
Here, however, the signals have units of volts rather than bits.


Interferometer-measured phase:
------------------------------
The interferometer-measured phase can be computed
from a `Lissajous` instance using
the `mitpci.interferometer.Phase` class e.g.

```python
Ph_int = mitpci.interferometer.Phase(L)

```

By default, a zero-delay, high-pass filter with
a 6-dB (in energy) knee at 10 kHz and
a transition width of 5 kHz is applied to the computed phase;
this removes both
low-frequency vibrational contributions and
bulk-plasma contributions
to the measured phase,
leaving (predominantly) the fluctuating, plasma-induced phase.
If *no* filtering is desired, simply specify
`Ph_int = mitpci.interferometer.Phase(L, filt=None)`.
Other zero-delay FIR filters can be designed with
my [filters](https://github.com/emd/filters) package and
applied to the measured phase using the `filt` keyword argument.


Interferometer-measured mode-number spectra:
--------------------------------------------
Spectral information can be readily computed and visualized using the
[random_data package](https://github.com/emd/random_data).
In particular, toroidal mode numbers can be measured by
correlating the MIT heterodyne interferometer
with the toroidally separated V2 interferometer, as follows:

```python
import random_data as rd

# Spectral-estimation parameters
Tens = 5e-3         # Ensemble time length, [Tens] = s
Nreal_per_ens = 5   # Number of realizations per ensemeble

# Determine mode numbers from interferometers
TorCorr_int = mitpci.interferometer.ToroidalCorrelation(
    Ph_int, Tens=Tens, Nreal_per_ens=Nreal_per_ens)

TorCorr_int.plotModeNumber(
    xlabel='$t \, [\mathrm{s}]$',
    ylabel='$f \, [\mathrm{Hz}]$',
    flim=[0, 300e3],
    all_positive=True)

```

![interferometer_mode_numbers](https://raw.githubusercontent.com/emd/mitpci/master/figs/interferometer_mode_numbers.png)

Note that the 45-degree toroidal separation of the V2 and MIT
interferometers allows identification of 8 distinct
toroidal mode numbers. The exact mode-number range, however,
depends on the mode's rotation:

* for unknown rotation, set `all_positive` to `False`,
  which will plot mode numbers between -3 <= n <= 4;

* for positive rotation (counterclockwise when viewing
  the vacuum vessel from above), set `all_positive` to
  `True`, which will plot mode numbers between 0 <= n <= 7.

The interferometer-measured toroidal mode numbers can be
compared to magnetic measurements.
For example, toroidal mode numbers can be readily extracted
from magnetic measurements using my
[magnetics package](https://github.com/emd/magnetics), as follows:

```python
import magnetics

TorSigs = magnetics.signal.ToroidalSignals(shot, tlim=tlim)

A = rd.array.Array(
    TorSigs.x, TorSigs.locations, Fs=TorSigs.Fs, t0=TorSigs.t0,
    Tens=Tens, Nreal_per_ens=Nreal_per_ens)

fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)

cmap_mag = magnetics.colormap.positive_mode_numbers()[0]

A.plotModeNumber(
    ax=axes[1],
    title='magnetics',
    xlabel='$t \, [\mathrm{s}]$',
    ylabel='',
    cblabel='$n$',
    mode_number_lim=[0, 11],
    cmap=cmap_mag)

TorCorr_int.plotModeNumber(
    ax=axes[0],
    title='interferometers',
    xlabel='$t \, [\mathrm{s}]$',
    ylabel='$f \, [\mathrm{Hz}]$',
    flim=[0, 300e3],
    all_positive=True)

```

![interferometer_and_magnetic_mode_numbers](https://raw.githubusercontent.com/emd/mitpci/master/figs/interferometer_and_magnetic_mode_numbers.png)


PCI-measured phase:
-------------------
The PCI-measured phase from each digitizer channel
can be easily retrieved via the `mitpci.pci.Phase` class e.g.

```python
shot = 171521
tlim = [1.2, 1.5]  # [tlim] = s
digitizer_channels = np.arange(16) + 1

Ph_pci = mitpci.pci.Phase(shot, digitizer_channels, tlim=tlim)

```

Even though all 16 digitizer channels are specified above,
the `mitpci.pci.Phase` class will *only* retrieve data
from the subset of `digitizer_channels` that are actually
mapped to functioning PCI channels.
By default, a zero-delay, high-pass filter with
a 6-dB (in energy) knee at 10 kHz and
a transition width of 5 kHz is applied to the computed phase
(this same filter is applied by default to the interferometer data).
If *no* filtering is desired, simply specify
`Ph_pci = mitpci.pci.Phase(..., filt=None)`.
Other zero-delay FIR filters can be designed with
my [filters](https://github.com/emd/filters) package and
applied to the measured phase using the `filt` keyword argument.
As the PCI only measures the plasma-induced phase
up to a calibration constant,
a calibration constant (`rad_per_bit`, in the class's `__init__()`)
from a recent PCI-interferometer audio cross-calibration
is applied to the PCI signal.


PCI complex-valued, spatial cross-correlation function:
-------------------------------------------------------
Typically, the digitized PCI signal corresponds to measurements
from nonuniformly spaced detector elements
(due to a combination of limited digitizer channels/capacity and
damaged detector elements).
Fortunately, the correlation function
(from which the two-dimensional autospectral density can be computed)
can still be estimated from nonuniformly spaced samples.
For example, the PCI complex-valued, spatial cross-correlation function
can be easily estimated and visualized as follows
(using the `Ph_pci` instance from above):

```python
import matplotlib.pyplot as plt

# Spectral-estimation parameters:
# -------------------------------
tlim = [1.33, 1.38]  # [tlim] = s; this defines ensemble's temporal bounds
Nreal_per_ens = 500  # number of realizations in the ensemble

# Compute and plot correlation function:
# --------------------------------------
corr = mitpci.pci.ComplexCorrelationFunction(
    Ph_pci, tlim=tlim, Nreal_per_ens=Nreal_per_ens)

corr.plotNormalizedCorrelationFunction()
plt.show()

```

![pci_correlation_function](https://raw.githubusercontent.com/emd/mitpci/master/figs/pci_correlation_function.png)

Note that this is a plot of the *normalized* complex-valued, spatial
cross-correlation function. That is, at each frequency `f`,
the correlation function `Gxy(delta, f)` is normalized to `Gxy(0, f)`
(i.e. its value at zero separation (`delta = 0`) and frequency `f`).
This allows for easy visualization of the correlation function's structure
even if `Gxy(delta, f)` varies by several orders of magnitude as
the frequency varies.
Note further the conspicuous gaps at various element separations;
for the currently digitized detector elements
(viewable from `Ph_pci.detector_elements`),
no such correlation pairs correspond to such separations, so
we *cannot* estimate the correlation function at these separations.
The `mitpci.pci.ComplexCorrelationFunction` class provides
a method to linearly interpolate across these gaps.
For example, to linearly interpolate across all gaps
less than or equal to 2 detector-element spacings, use

```python
corr.interpolate(2)

```

Careful inspection of the normalized correlation function, however,
reveals that the correlation function can vary substantially
over a detector-element separation or two, so
such linear interpolation should not be accepted blindly.
(Presumably, which detector elements are digitized
could be optimized to provide better coverage
of the correlation function and prevent such gaps;
such investigations are relegated to future work).
If the interpolation produces undesired results,
the previous correlation function can be restored via
`corr.unInterpolate()`.


PCI two-dimensional autospectral density estimate:
--------------------------------------------------
The two-dimensional autospectral density can be estimated
from the complex correlation function as follows
(proceed with the correlation function from above
after having linearly interpolated across all gaps
that are less than or equal to 2 detector-element spacings):

```python
# Compute 2d autospectral density. Temporal estimates are
# via Welch's method of overlapped periodograms, while
# spatial estimates are via a Burg autoregressive method
# with `p` poles and `Nk` points.
asd2d = mitpci.pci.TwoDimensionalAutoSpectralDensity(
    corr, spatial_method='burg', burg_params={'p': 5, 'Nk': 1000})

flim = [10, 1500]  # [flim] = kHz
asd2d.plotSpectralDensity(flim=flim)

```

![pci_2d_spectrum](https://raw.githubusercontent.com/emd/mitpci/master/figs/pci_2d_spectrum.png)

Spatial resolution typically increases with higher `p`, but
higher `p` can produce "pole splitting" and can introduce
other numerical artifacts. Spatial estimates can also be made
via a Fourier method (`spatial_method = 'fourier'` along
with its own associated `fourier_params`), but
the limited number of spatial samples typically
results in the Fourier estimates performing
much more poorly than the Burg estimates.
