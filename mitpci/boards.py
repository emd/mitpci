'''This module implements a class for characterizing the trigger offset
between the two boards of the mitpci digitizer.

'''


# Standard library imports
import numpy as np

# Related 3rd-party imports
import random_data as rd

# Intra-package imports
from .pci.signal import Phase


class TriggerOffset(object):
    '''
    '''
    def __init__(self, Ph_pci, Nreal_per_ens=100):
        '''

        For example, as of 9/7/17, the PCI's digitizer to detector
        mapping is:

                digitizer channel  <-->  detector element
                -----------------        ----------------
                   7 (board 7)                  15
                   8 (board 7)                  16
                   9 (board 8)                  17
                  12 (board 8)                  18

        such that the detector

               correlation pairs     <-->  element separation
            -----------------------        ------------------
             7 & 8  (board 7 only)                 1
             8 & 9  (board 7 & 8)                  1
             9 & 12 (board 8 only)                 1


        '''
        # Ensure we are working with correct object
        if type(Ph_pci) is not Phase:
            raise ValueError('`Ph_pci` must be of type %s' % Phase)

        self.shot = Ph_pci.shot
        self.Nreal_per_ens = Nreal_per_ens

        # Check topology and parse results
        res = _check_topology(
            Ph_pci.detector_elements, Ph_pci.digizer_board,
            d1='DT216_7', d2='DT216_8')

        self.topological_direction = res[0]
        self.ind_b7_min = res[1]
        self.ind_b7_max = res[2]
        self.ind_b8_min = res[3]
        self.ind_b8_max = res[4]

    def _getOffset(self, Ph_pci):
        'Get trigger offset of board 8 relative to board 7'
        Tens = Ph_pci.t()[-1] - Ph_pci.t0

        # Compute cross-spectral densities with ordering
        # as determined by `self.topological_direction`, where
        #
        # - csd_b7 = cross-spectral density from board 7 only
        # - csd_b78 = cross-spectral density from board 7 & 8, and
        # - csd_b8 = cross-spectral density from board 8 only.
        #
        if self.topological_direction > 0:
            csd_b7 = rd.spectra.CrossSpectralDensity(
                Ph_pci.x[self.ind_b7_min, :],
                Ph_pci.x[self.ind_b7_max, :],
                Fs=Ph_pci.Fs, t0=Ph_pci.t0,
                Tens=Tens, Nreal_per_ens=self.Nreal_per_ens)
            csd_b78 = rd.spectra.CrossSpectralDensity(
                Ph_pci.x[self.ind_b7_max, :],
                Ph_pci.x[self.ind_b8_min, :],
                Fs=Ph_pci.Fs, t0=Ph_pci.t0,
                Tens=Tens, Nreal_per_ens=self.Nreal_per_ens)
            csd_b8 = rd.spectra.CrossSpectralDensity(
                Ph_pci.x[self.ind_b8_min, :],
                Ph_pci.x[self.ind_b8_max, :],
                Fs=Ph_pci.Fs, t0=Ph_pci.t0,
                Tens=Tens, Nreal_per_ens=self.Nreal_per_ens)
        else:
            csd_b7 = rd.spectra.CrossSpectralDensity(
                Ph_pci.x[self.ind_b7_max, :],
                Ph_pci.x[self.ind_b7_min, :],
                Fs=Ph_pci.Fs, t0=Ph_pci.t0,
                Tens=Tens, Nreal_per_ens=self.Nreal_per_ens)
            csd_b78 = rd.spectra.CrossSpectralDensity(
                Ph_pci.x[self.ind_b7_min, :],
                Ph_pci.x[self.ind_b8_max, :],
                Fs=Ph_pci.Fs, t0=Ph_pci.t0,
                Tens=Tens, Nreal_per_ens=self.Nreal_per_ens)
            csd_b8 = rd.spectra.CrossSpectralDensity(
                Ph_pci.x[self.ind_b8_max, :],
                Ph_pci.x[self.ind_b8_min, :],
                Fs=Ph_pci.Fs, t0=Ph_pci.t0,
                Tens=Tens, Nreal_per_ens=self.Nreal_per_ens)

        # - explanation for offset
        # - weighting by:
        #     - coherence: threshold
        #     - normalized difference: threshold
        # - plotting

        return


def _check_topology(elements, domains, d1='1', d2='2'):
    '''Check the topology of `elements` in relation to `domains`.

    The checks are performed within the context of the needs of
    :py:class:`TriggerOffset <mitpci.boards.TriggerOffset>`

    Parameters:
    -----------
    elements - array_like, (4,)
        Elements that are mapped onto `domains`.
        A ValueError is raised if `elements` does not have length 4.
        [elements] = arbitrary units

    domains - array_like, (4,)
        An array of *strings* that denotes the domains that `elements`
        are mapped onto.

        `domains` may only have 2 values -- those specified by `d1` and
        `d2`. Further, `d1` and `d2` must appear exactly twice in
        `domains`. If any of the above conditions are not met, a
        ValueError is raised.

        Finally, if the elements belonging to each domain are
        *not* simply connected, a ValueError is raised.

    d1 - string
        The first allowable value of `domains`.

    d2 - string
        The second allowable value of `domains`.

    Returns:
    --------
    res - tuple, where

    res[0] = topological_direction - int
        Values may be +/-1, with a positive value indicating that
        `elements` increases when moving from domain `d1` to `d2`
        (and a negative value indicating that `elements` decreases
        when moving from domain `d1` to `d2`).

    res[1] = ind_d1_min - int
        `elements[ind_d1_min]` gives the minimum value of `elements`
        that is mapped onto domain `d1`.

    res[2] = ind_d1_max - int
        `elements[ind_d1_max]` gives the maximum value of `elements`
        that is mapped onto domain `d1`.

    res[3] = ind_d2_min - int
        `elements[ind_d2_min]` gives the minimum value of `elements`
        that is mapped onto domain `d2`.

    res[4] = ind_d2_max - int
        `elements[ind_d2_max]` gives the maximum value of `elements`
        that is mapped onto domain `d2`.

    '''
    if len(elements) != 4:
        raise ValueError('`elements` must have length 4')
    elif len(domains) != 4:
        raise ValueError('`domains` must have length 4')
    elif np.sum(domains == d1) != 2:
        raise ValueError('`d1` must appear exactly 2x in `domains`')
    elif np.sum(domains == d2) != 2:
        raise ValueError('`d2` must appear exactly 2x in `domains`')

    # Get indices corresponding to domains 1 and 2
    ind_d1 = np.where(domains == d1)[0]
    ind_d2 = np.where(domains == d2)[0]

    # Determine the minimum and maximum elements
    # that are mapped onto each domain
    d1_min = np.min(elements[ind_d1])
    d1_max = np.max(elements[ind_d1])
    d2_min = np.min(elements[ind_d2])
    d2_max = np.max(elements[ind_d2])

    # Check topology. We want simply connected domains; for example
    #
    #       elements:   0 | 1 | 2 | 3
    #       --------------------------
    #       domains:      1   |   2
    #
    # or
    #
    #       elements:   3 | 2 | 1 | 0
    #       --------------------------
    #       domains:      1   |   2
    #
    # However, non-simply connected domains, for example, looks like
    #
    #       elements:   0 | 1 | 2 | 3
    #       --------------------------
    #       domains:    1 | 2 | 1 | 2
    #
    # from which the following "if-then" logic readily follows.
    if (d1_max > d2_min) and (d2_max > d1_min):
        raise ValueError('`elements` not simply connected to `domains`')

    # Define a `topological direction` such that a positive value
    # corresponds to `elements` increasing as we move from domain
    # `d1` to domain `d2`
    topological_direction = np.sign(d2_min - d1_max)

    if topological_direction > 0:
        spacing = np.array([
            d1_max - d1_min,
            d2_min - d1_max,
            d2_max - d2_min])
    else:
        spacing = np.array([
            d1_max - d1_min,
            d1_min - d2_max,
            d2_max - d2_min])

    # Ensure that `spacing` is uniform
    if len(np.unique(spacing)) != 1:
        raise ValueError('Correlation pairs do not have uniform spacing')

    # Get indices corresponding to minimum and maximum
    # elements that are mapped onto each domain
    ind_d1_min = np.where(elements == d1_min)[0][0]
    ind_d1_max = np.where(elements == d1_max)[0][0]
    ind_d2_min = np.where(elements == d2_min)[0][0]
    ind_d2_max = np.where(elements == d2_max)[0][0]

    return (topological_direction,
            ind_d1_min, ind_d1_max,
            ind_d2_min, ind_d2_max)
