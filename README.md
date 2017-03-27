Python tools for retrieving signals digitized by the mitpci system.


Background:
===========
The mitpci system digitizes data from MIT's (a) phase contrast imaging (PCI)
and (b) heterodyne interferometer systems. The systems share a beam path
through the DIII-D vessel and a digitizer, but they use different
interference schemes and detectors to make their respective measurements.

As these systems share a digitization system, the mechanics underlying
signal retrieval from both systems are identical. This module (`mitpci`)
aims to provide Python tools for flexible signal retrieval and analysis.


Installation:
=============

... on GA's Iris cluster:
-------------------------
Package management is cleanly handled on Iris via
[modules](https://diii-d.gat.com/diii-d/Iris#Environment_modules).
The `mitpci` package has a corresponding modulefile
[here](https://github.com/emd/modulefiles).

To use the `mitpci` package, change to the directory
you'd like to download the source files to and
retrieve the source files from github by typing

    $ git clone https://github.com/emd/mitpci.git

The created `mitpci` directory defines the
package's top-level directory. The `mitpci` package
depends on four additional packages:

* [random_data](https://github.com/emd/random_data),
* [bci](https://github.com/emd/bci),
* [magnetics](https://github.com/emd/magnetics), and
* [distinct_colours](https://personal.sron.nl/~pault/).

These packages and their modulefiles should be similarly cloned.

Now, at the top of the corresponding
[modulefile](https://github.com/emd/modulefiles/blob/master/mitpci),
there is a TCL variable named `mitpci_root`;
this must be altered to point at the
top-level directory of the cloned `mitpci` package.
Similarly, the TCL variable `modulefiles_dir`
must be altered to point at the directory containing
the modulefiles for `random_data`, `bci`, `magnetics`, and
`distinct_colours`.
That's it! You shouldn't need to change anything else in
the modulefile. The `mitpci` module can
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

... elsewhere:
--------------
Define an environmental variable `$pci_path` specifying
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
The MIT heterodyne interferometer's in-phase (I) and quadrature (Q) signals
can be readily loaded as follows:

```python
import mitpci

# Load data from MIT interferometer
shot = 167342
tlim = [1.0, 2.5]  # [tlim] = s
D = mitpci.interferometer.Demodulated(shot, tlim=tlim)

```
Note that the `Demodulated` class automatically compensates for
demodulator imperfections (i.e. DC offsets and amplitude imbalances
between I and Q).

Spectral information can then be readily computed and visualized using the
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
xcorr_int = mitpci.interferometer.ToroidalCorrelation(
    D, Tens=Tens, Nreal_per_ens=Nreal_per_ens)

xcorr_int.plotModeNumber(
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
from magnetic measurements using the
[magnetics package](https://github.com/emd/magnetics), as follows:

```python
import magnetics

torsigs = magnetics.signal.ToroidalSignals(shot, tlim=tlim)

A = rd.array.Array(
    torsigs.x, torsigs.locations, Fs=torsigs.Fs, t0=torsigs.t0,
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

xcorr_int.plotModeNumber(
    ax=axes[0],
    title='interferometers',
    xlabel='$t \, [\mathrm{s}]$',
    ylabel='$f \, [\mathrm{Hz}]$',
    flim=[0, 300e3],
    all_positive=True)

```

![interferometer_and_magnetic_mode_numbers](https://raw.githubusercontent.com/emd/mitpci/master/figs/interferometer_and_magnetic_mode_numbers.png)
