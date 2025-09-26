# spectroscopy

This repository contains Python scripts for analyzing stellar spectra recorded by the NRES instruments installed on Las Cumbres Observatory 1-meter telescopes.

The code was developed during my 2025 Summer Project, funded by the Maricopa Community Colleges. As a continuation of this project, example notebooks will be developed and a set of webpages will be written to introduce concepts of stellar spectra to interested students and faculty. The scripts are heavily commented.

The script "nres_data.py" loads NRES FITS files into memory.

The script "phoenix_data.py" loads PHOENIX synthetic spectra from the FTP server into memory.

The script "spec_tools.py"
* normalizes intensities of PHOENIX spectra
* cross-matches NRES spectra with PHOENIX spectra to derive radial velocity
* corrects radial velocity at telescope to Solar System barycenter velocity and date

Future implementations:
* searching the effective temperature (Teff) and metallicity (Z) grid of PHOENIX spectra
* correct telescope radial velocities using observations of standard stars
* automatic plotting of radial velocity curves, suitable for publication
* simultaneous cross-matching NRES spectra to fit double-lined stellar spectra
