try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'mitpci',
    'version': '0.3',
    'packages': ['mitpci'],
    'install_requires': [
        'numpy', 'matplotlib', 'scipy', 'MDSplus', 'nose'],
        # 'random_data', 'bci', 'magnetics'],  # <- Not on PyPI
    'author': 'Evan M. Davis',
    'author_email': 'emd@mit.edu',
    'url': '',
    'description': 'Interact with signals from `mitpci` digitizer system'
}

setup(**config)
