from setuptools import setup, find_packages
import logging
logger = logging.getLogger(__name__)

version = '0.1.0'

try:
    with open('README.md', 'r') as f:
        long_desc = f.read()
except:
    logger.warning('Could not open README.md.  long_description will be set to None.')
    long_desc = None

setup(
    name = 'neutralatomcompilation',
    packages = find_packages(),
    version = version,
    description = 'compiler for emerging Neutral Atom / Rydberg QC',
    long_description = long_desc,
    long_description_content_type = 'text/markdown',
    #author = '',
    #author_email = '',
    url = '',
    download_url = '',
    keywords = ['quantum computing', 'neutral atoms', 'compiler'],
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    install_requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'qutip',
        'qiskit',
        'networkx',
    ],
    extras_require = {
    },
)

