try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'mitpci',
    'version': '0.1',
    'packages': ['mitpci'],
    'install_requires': ['numpy', 'MDSplus', 'nose'],
    'author': 'Evan M. Davis',
    'author_email': 'emd@mit.edu',
    'url': '',
    'description': 'Interact with signals from `mitpci` digitizer system'
}

setup(**config)
