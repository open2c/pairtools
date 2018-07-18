Quickstart
==========

Installation
------------

Requirements:

- python 3.4 and higher
- unix sort
- bgzip
- pbgzip (optional)
- samtools >= 1.4
- Python packages ``Cython``, ``numpy``,  ``click``, ``nose``

We highly recommend using the conda package manager to install scientific packages and their dependencies. To get ``conda``, you can download either the full `Anaconda <https://www.continuum.io/downloads>`_ Python distribution which comes with lots of data science software or the minimal `Miniconda <http://conda.pydata.org/miniconda.html>`_ distribution which is just the standalone package manager plus Python. In the latter case, you can install the packages as follows:

::

    $ conda install samtools sort bgzip pbgzip Cython numpy pip
    $ pip install click nose


Install the latest version of pairtools using pip.

::

    $ pip install git+https://github.com/mirnylab/pairtools

