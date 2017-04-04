wcsMatrixSolver

Uses matrix math to quickly generate wcs solutions from an input file. Not yet designed for general use, but getting there. If you need help, email the author.

1. Installation
===============

1a. Dependencies
----------------

Correct function depends on a lot of things:
* numpy
* astropy
* trippy
* matplotlib
* numdisplay (from stsci)

You also need SExtractor available on the command line. This can often be installed with through standard linux repository systems (yum, apt-get, etc.) or can be simply downloaded and compiled. [Source is available here].(http://www.astromatic.net/software/sextractor)


1b. Installation
----------------

For now, just clone the repo, and make the destination available in your PYTHONPATH and PATH environment variables. Eventually this will probable go to pypi.

2. Usage
========

Basic usage is
    > wcsMatrixSolver.py fitsFileName

then provide appropriate other options (use -h to see these). If image is in the region of the tsv you load, all is well.

Again, email the author. This is left ambiguous on purpose... for now...
