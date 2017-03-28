Python tools for retrieving signals digitized by the mitpci system.


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
Currently, there are no specialized analysis routines for the PCI
implemented in this module.
However, there are several specialized analysis routines for
the heterodyne interferometer, including:

* compensation of demodulator imperfections (i.e. DC offsets and
  amplitude imbalances between the I and Q signals), and
* correlation with the toroidally separated V2 interferometer
  for toroidal mode-number identification.

The use of these routines is discussed below.


Installation:
=============

... on GA's Iris cluster:
-------------------------
To use the `mitpci` package, change to the directory
you'd like to download the source files to, e.g.

    $ cd $HOME/python/projects

and retrieve the source files from github by typing

    $ git clone https://github.com/emd/mitpci.git

The created `mitpci` directory defines the
package's top-level directory. The `mitpci` package
depends on three additional packages:

* [random_data](https://github.com/emd/random_data),
* [bci](https://github.com/emd/bci), and
* [magnetics](https://github.com/emd/magnetics).

These packages should similarly be cloned, e.g.

    $ git clone https://github.com/emd/random_data.git
    $ git clone https://github.com/emd/bci
    $ git clone https://github.com/emd/magnetics

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
the modulefiles for `random_data`, `bci`, and `magnetics`.
That's it! You shouldn't need to change anything else in
the `mitpci` modulefile.
(Of course, similar changes will also need to be made
to the `<package>_root` TCL variable in the modulefile
for the `random_data`, `bci`, and `magnetics` modules).
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

... elsewhere:
--------------
First, install the `mitpci` dependencies that are *not* on PyPI:

* [random_data](https://github.com/emd/random_data),
* [bci](https://github.com/emd/bci), and
* [magnetics](https://github.com/emd/magnetics).

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

However, a specialized `mitpci.interferometer.Demodulated` class
also exists for retrieval of the heterodyne interferometer's
in-phase (I) and quadrature (Q) signals, e.g.

```python
D = mitpci.interferometer.Demodulated(shot, tlim=tlim)

```

Note that the `mitpci.interferometer.Demodulated` class
automatically compensates for demodulator imperfections
(i.e. DC offsets and amplitude imbalances between I and Q).
The post-processed I and Q objects are accessed via
`D.I` and `D.Q`, respectively.
Both `D.I` and `D.Q` are *instances* of the `mitpci.signal.Signal` class, so
the raw digitized signal (in bits) is accessed via e.g. `D.I.x`, and
the timebase can be generated via e.g. `D.I.t()`.
The phase is readily calculated via `D.getPhase()`.

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
from magnetic measurements using my
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
