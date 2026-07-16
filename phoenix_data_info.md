# class PhoenixTemplate

Represent a single PHOENIX synthetic stellar spectrum

Methods to load and normalize PHOENIX stellar spectrum template from FITS files

Files are stored in a shared Google Drive folder accessible to all MCC faculty and students

PHOENIX spectrum comes in two parts:
* a -wavelength- file common to all PHOENIX stellar spectra
* a -flux- file with unique Teff, logg, [Fe/H]

Spectral fluxes are normalized (range of 0...1) using an approximation of
continuum which is determined from a rolling maximum of flux array

## Attributes:

* filepath (str) -- directory path to PHOENIX files<br>
    (default is '/content/drive/MyDrive/astdata/PHOENIX Spectra/')

* filename (str) -- name of PHOENIX flux file to be loaded

* wavelengthfile (str) -- name of PHOENIX wavelength file to be loaded<br>
    (default is 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

* wave_header (FITS header) -- FITS header of wavelength template file

* wavelengths (float array) -- wavelength basis for other spectral arrays

* data_header (FITS header) -- FITS header of flux template file

* optical (Boolean array) -- mask of optical wavelengths (3500-9500 A)

* fluxes (float array) -- unnormalized synthetic fluxes for model star

* continuum (float array) -- estimated continuum spectrum

* norm_fluxes (float array) -- continuum normalized fluxes 

* Teff (float) -- effective temperature of PHOENIX model star

* logg (float) -- log of surface gravity of PHOENIX model star

* FeH (float) -- metallicity of PHOENIX model star

## Methods:

* read_wavefile() -- Read wavelength file for PHOENIX spectrum<br>
    returns wave_header, optical, wavelengths

* read_fluxfile() -- Read flux file for PHOENIX spectrum<br>
    returns data_header, fluxes, Teff, logg, FeH

* normalize_fluxes() -- Estimates continuum, then normalizes fluxes<br>
    returns continuum, norm_fluxes