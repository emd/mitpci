Python tools for retrieving signals digitized by the mitpci system.


Background:
-----------
The mitpci system digitizes data from MIT's (a) phase contrast imaging (PCI)
and (b) heterodyne interferometer systems. The systems share a beam path
through the DIII-D vessel and a digitizer, but they use different
interference schemes and detectors to make their respective measurements.

As these systems share a digitization system, the mechanics underlying
signal retrieval from both systems are identical. This module (`mitpci`)
aims to provide Python tools for flexible signal retrieval.


Installation:
-------------
Define an environmental variable `$pci_path` specifying
the appropriate MDSplus servers' tree-path definitions
(`mitpci.gat.com::/trees/pci` and `hermit.gat.com::/trees/pci`)
by, for example, adding the following to your `.bashrc`

    $ export pci_path=mitpci.gat.com::/trees/pci;hermit.gat.com::/trees/pci

(Data is digitized and first available on `mitpci`, but
it is later transferred to `hermit` for long(er)-term storage).

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
the `mitpci` module and its associated tests requires access to the
mitpci MDSplus server, which is located behind GA's firewall; if you
are not within the firewall, the module will err and the tests will fail
due to the inability to read from the mitpci MDSplus server.
