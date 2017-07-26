# Standard library imports
import numpy as np
import matplotlib.pyplot as plt

# Related 3rd-party imports
import random_data as rd

# Intra-package imports
from .signal import Phase


def check_centering(shot, tlim, channels=np.arange(3, 17)):
    fig, axes = plt.subplots(4, 4, sharex=True, sharey=True)

    Tens = 5e-3
    Nreal_per_ens = 10
    flim = [10e3, 2e6]

    # Reference channel, in middle of detector
    ch_ref = 9
    Ph = Phase(shot, ch_ref, tlim=tlim)
    asd_ref = rd.spectra.AutoSpectralDensity(
        Ph.x, Fs=Ph.Fs, t0=Ph.t0,
        Tens=Tens, Nreal_per_ens=Nreal_per_ens)

    find = np.where(np.logical_and(
        asd_ref.f >= flim[0],
        asd_ref.f <= flim[1]))[0]

    # Loop through each channel and plot relative to reference channel
    for i, ch in enumerate(channels):
        Ph = Phase(shot, ch, tlim=tlim)
        asd = rd.spectra.AutoSpectralDensity(
            Ph.x, Fs=Ph.Fs, t0=Ph.t0,
            Tens=Tens, Nreal_per_ens=Nreal_per_ens)

        row = i % 4
        col = i // 4

        axes[row, col].loglog(
            asd_ref.f[find],
            np.mean(asd_ref.Gxx[find, :], axis=-1),
            c='r')
        axes[row, col].loglog(
            asd.f[find],
            np.mean(asd.Gxx[find, :], axis=-1),
            c='k')

        axes[row, col].set_title('ch. % i' % ch)

    labels00 = [
        'reference (ch. %i)' % ch_ref,
        'ch. %i' % channels[0]
    ]

    axes[0, 0].legend(labels00, loc='upper right')
    axes[0, 0].set_xlim(flim)
    axes[0, 0].set_ylim([1e-16, 1e-8])

    plt.show()

    return
