from nose import tools
import numpy as np
import mitpci


def test__check_topology():
    d1 = '1'
    d2 = '2'

    # Tests for blatantly wrong inputs:
    # ---------------------------------
    elements = np.arange(4)
    domains = np.array([d1, d1, d2, d2])

    # Incorrect length for argument passed to `elements`
    tools.assert_raises(
        ValueError,
        mitpci.boards._check_topology,
        *[elements[:-1], domains],
        **{'d1': d1, 'd2': d2})

    # Incorrect length for argument passed to `domains`
    tools.assert_raises(
        ValueError,
        mitpci.boards._check_topology,
        *[elements, domains[:-1]],
        **{'d1': d1, 'd2': d2})

    # Incorrect amount of `d1` in argument passed to `domains`
    domains_corrupted = domains.copy()
    domains_corrupted[2] = d1
    tools.assert_raises(
        ValueError,
        mitpci.boards._check_topology,
        *[elements, domains_corrupted],
        **{'d1': d1, 'd2': d2})

    # Incorrect amount of `d2` in argument passed to `domains`
    domains_corrupted = domains.copy()
    domains_corrupted[1] = d2
    tools.assert_raises(
        ValueError,
        mitpci.boards._check_topology,
        *[elements, domains_corrupted],
        **{'d1': d1, 'd2': d2})

    # Tests for simply connected topology:
    # ------------------------------------
    elements = np.arange(4)
    domains = np.array([d1, d2, d1, d2])  # *not* simply connected
    tools.assert_raises(
        ValueError,
        mitpci.boards._check_topology,
        *[elements, domains],
        **{'d1': d1, 'd2': d2})

    # Tests for uniform spacing:
    # --------------------------
    elements = np.array([0, 1, 3, 4])     # non-uniform spacing
    domains = np.array([d1, d1, d2, d2])  # simply connected
    tools.assert_raises(
        ValueError,
        mitpci.boards._check_topology,
        *[elements, domains],
        **{'d1': d1, 'd2': d2})

    # Tests for correct return values:
    # --------------------------------
    elements = np.arange(4)
    domains = np.array([d1, d1, d2, d2])

    # Test 1: baseline
    res = mitpci.boards._check_topology(
        elements, domains, d1=d1, d2=d2)

    topological_direction = res[0]
    ind_d1_min = res[1]
    ind_d1_max = res[2]
    ind_d2_min = res[3]
    ind_d2_max = res[4]

    tools.assert_equal(topological_direction, 1)
    tools.assert_equal(ind_d1_min, 0)
    tools.assert_equal(ind_d1_max, 1)
    tools.assert_equal(ind_d2_min, 2)
    tools.assert_equal(ind_d2_max, 3)

    # Test 2: reverse `elements` relative to baseline
    res = mitpci.boards._check_topology(
        elements[::-1], domains, d1=d1, d2=d2)

    topological_direction = res[0]
    ind_d1_min = res[1]
    ind_d1_max = res[2]
    ind_d2_min = res[3]
    ind_d2_max = res[4]

    tools.assert_equal(topological_direction, -1)
    tools.assert_equal(ind_d1_min, 1)
    tools.assert_equal(ind_d1_max, 0)
    tools.assert_equal(ind_d2_min, 3)
    tools.assert_equal(ind_d2_max, 2)

    # Test 3: reverse `domains` relative to baseline
    res = mitpci.boards._check_topology(
        elements, domains[::-1], d1=d1, d2=d2)

    topological_direction = res[0]
    ind_d1_min = res[1]
    ind_d1_max = res[2]
    ind_d2_min = res[3]
    ind_d2_max = res[4]

    tools.assert_equal(topological_direction, -1)
    tools.assert_equal(ind_d1_min, 2)
    tools.assert_equal(ind_d1_max, 3)
    tools.assert_equal(ind_d2_min, 0)
    tools.assert_equal(ind_d2_max, 1)

    # Test 4: reverse `elements` and `domains` relative to baseline
    res = mitpci.boards._check_topology(
        elements[::-1], domains[::-1], d1=d1, d2=d2)

    topological_direction = res[0]
    ind_d1_min = res[1]
    ind_d1_max = res[2]
    ind_d2_min = res[3]
    ind_d2_max = res[4]

    tools.assert_equal(topological_direction, 1)
    tools.assert_equal(ind_d1_min, 3)
    tools.assert_equal(ind_d1_max, 2)
    tools.assert_equal(ind_d2_min, 1)
    tools.assert_equal(ind_d2_max, 0)

    return


# def test__checkChannels():
#     shot = 171521
#     tlim = [1.00, 1.05]  # [tlim] = s
# 
#     # Pass incorrect argument type to `TriggerOffset`
#     digitizer_channels = [7, 8, 9, 12]
#     Ph_pci = mitpci.pci.Phase(shot, digitizer_channels, tlim=tlim)
#     tools.assert_raises(
#         ValueError,
#         mitpci.boards.TriggerOffset,
#         *[Ph_pci.x])
# 
#     # Too many channels from board 7
#     digitizer_channels = [6, 7, 8, 9, 12]
#     Ph_pci = mitpci.pci.Phase(shot, digitizer_channels, tlim=tlim)
#     tools.assert_raises(
#         ValueError,
#         mitpci.boards.TriggerOffset,
#         *[Ph_pci])
# 
#     # Too many channels from board 8
#     digitizer_channels = [7, 8, 9, 12, 13]
#     Ph_pci = mitpci.pci.Phase(shot, digitizer_channels, tlim=tlim)
#     tools.assert_raises(
#         ValueError,
#         mitpci.boards.TriggerOffset,
#         *[Ph_pci])
# 
#     return
