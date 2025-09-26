# spectroscopy

This repository contains Python scripts and an example Google Colab (Jupyter) notebook for analyzing stellar spectra recorded by the NRES instruments installed on Las Cumbres Observatory 1-meter telescopes.

The code was developed during my 2025 Summer Project, funded by the Maricopa Community Colleges. As a continuation of this project, more notebooks will be developed and a set of webpages will be developed to introduce concepts of stellar spectra to interested students and faculty. The scripts and notebook are commented for myself and for others who wish to understand the methods I used.

The Python scripts currently handle the following:
* loading NRES FITS files into memory
* loading PHOENIX synthetic spectra from the FTP server into memory
* normalizing intensities of PHOENIX spectra
* cross-matching NRES spectra with PHOENIX spectra to derive radial velocity
* correcting radial velocity at telescope to Solar System barycenter velocity and date

Future implementations:
* searching the effective temperature (Teff) and metallicity (Z) grid of PHOENIX spectra
* correct telescope radial velocities using observations of standard stars
* automatic plotting of radial velocity curves, suitable for publication
* simultaneous cross-matching NRES spectra to fit double-lined stellar spectra
