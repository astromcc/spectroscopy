# spectroscopy

This repository contains Python code and other resources for analyzing stellar spectra recorded by the NRES instruments on the Las Cumbres Observatory 1-meter telescopes.

A set of webpages will be developed to introduce concepts of stellar spectra to interested students and faculty.

The Python scripts currently handle the following:
* loading NRES FITS files into memory
* loading PHOENIX synthetic spectra from the FTP server into memory
* normalizing intensities of PHOENIX spectra to match NRES data
* cross-matching NRES spectra with PHOENIX spectra to derive radial velocity
* correcting radial velocity at telescope to Solar System barycenter velocity and date

Future implementations:
* searching the effective temperature (Teff) and metallicity (Z) grid of PHOENIX spectra
* correct telescope radial velocities using observations of standard stars
* automatic plotting of radial velocity curves, suitable for publication
* simultaneous cross-matching NRES spectra to fit double-lined stellar spectra
