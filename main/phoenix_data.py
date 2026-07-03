
'''A class to represent a single PHOENIX synthetic stellar spectrum.

Load and normalize PHOENIX stellar spectrum template from FITS files.
Files are stored in a shared Google Drive folder accessible to all
MCC faculty and students.

PHOENIX spectrum comes in two parts: a -wavelength- file common to all PHOENIX
stellar spectra, and a -flux- file with unique Teff, logg, [Fe/H].

Spectral fluxes are normalized (0..1) using an approximation of
continuum which is determined from a rolling maximum of flux array.

Attributes
----------
filepath (str) -- directory path to PHOENIX files
    (default is '/content/drive/MyDrive/astdata/PHOENIX Spectra/')

filename (str) -- name of PHOENIX flux file to be loaded

wavelengthfile (str) -- name of PHOENIX wavelength file to be loaded
    (default is 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

wave_header (FITS header) -- FITS header of wavelength template file

wavelengths (float array) -- wavelength basis for other spectral arrays

data_header (FITS header) -- FITS header of flux template file

optical (Boolean array) -- mask of optical wavelengths (3500-9500 A)

fluxes (float array) -- unnormalized synthetic fluxes for model star

continuum (float array) -- estimated continuum spectrum

norm_fluxes (float array) -- continuum normalized fluxes 

Teff (float) -- effective temperature of PHOENIX model star

logg (float) -- log of surface gravity of PHOENIX model star

FeH (float) -- metallicity of PHOENIX model star


Methods
-------
read_wavefile() -- Read wavelength file for PHOENIX spectrum.
    returns wave_header, optical, wavelengths

read_fluxfile() -- Read flux file for PHOENIX spectrum.
    returns data_header, fluxes, Teff, logg, FeH

normalize_fluxes() -- Estimates continuum, then normalizes fluxes.
    returns continuum, norm_fluxes


*** Author: Kevin Healy, Mesa CC
'''

import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import splrep, splev

from astropy import units as u
from astropy.io import fits
from astropy.modeling.models import BlackBody

class PhoenixTemplate:
    # Retrieve PHOENIX template spectra from FTP server
    # Spectrum in two parts: wavelength grid in a standard file
    # and fluxes for specific stellar parameters in another file
    # Wavelength file relevant for ALL PHOENIX spectrum templates.

    # From Husser et al. (2013) article, filename format is:
    # Teff is 5 digits, padded with leading zero
    #      = 02300 - 07000 K, step = 100 K
    #      = 07000 - 12000 K, step = 200 K
    # logg is float w/ 2 decimal places
    #      =  0.00 -  6.00, step = 0.5
    # For [Fe/H] = 0.0, not all Teff, logg models are available! 
    #         0.00 available for Teff = 2300 -  5600 K
    #         0.50 available for Teff = 2300 -  5900 K
    #         1.00 available for Teff = 2300 -  8200 K
    #         1.50 available for Teff = 2300 -  9200 K
    #         2.00 available for Teff = 2300 - 12000 K
    #         all higher logg's available for all Teff

    # [Fe/H] is signed float w/ 1 decimal place
    #        = -4.0 - -2.0, step = 1.0
    #        = -2.0 - +1.0, step = 0.5
    # [a/Fe] adds additional string to signify alpha
    #        = -0.2 - +1.2, step = 0.2

    ### Assume, for now, that [Fe/H] and [a/Fe] are both 0.0
    # TO-DO: Permit other input values of [Fe/H] and [a/Fe]

    # Spectra are sampled at:
    #   500 -  3000 A, step = 0.1 A
    #  3000 - 25000 A, R ~ 500,000
    # 25000 - 55000 A, R ~ 100,000

    def __init__(self, Teff, logg,
            # Set default server, path, and filenames for [Fe/H] = 0.0
            FeH = 0.0,
            ftp_server = 'ftp://phoenix.astro.physik.uni-goettingen.de/',
            wave_path  = 'HiResFITS/',
            flux_path  = 'HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/',
            wave_file = 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',
            flux_file_sfx = ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"):
        self.Teff = Teff
        self.str_Teff = None
        self.logg = logg
        self.str_logg = None
        self.FeH = FeH
        self.str_FeH = "0.0"
        self.ftp_server = ftp_server
        self.wave_file = wave_file
        self.wave_uri = ftp_server + wave_path + wave_file
        self.flux_path = flux_path
        self.flux_file_sfx = flux_file_sfx
        self.flux_uri = None
        self.wave_header = None
        self.wavelengths = None
        self.data_header = None
        self.optical = None
        self.fluxes = None
        self.norm_fluxes = None

    def format_Teff(self, debug=False):
        '''Create formatted Teff string for filename for PHOENIX flux file.
        debug=True to print formatted string.
        '''
        min_Teff =  2300
        cut_Teff =  7000
        max_Teff = 12000
        temp = self.Teff

        if temp < min_Teff:
            self.str_Teff = "02300"
            print("Lowest available Teff for PHOENIX spectra is 2300 K.")
        elif min_Teff <= temp <= cut_Teff:
            # Round to nearest 100 K, then pad with zeroes to 5 characters
            rounded = round(self.Teff, -2)
            self.str_Teff = f"{int(rounded):05d}"
        elif cut_Teff < temp <= max_Teff:
            # Round to nearest 200 K, then pad with zeroes to 5 characters
            rounded = round(temp / 200.0) * 200.0
            self.str_Teff = f"{int(rounded):05d}"
        elif temp > max_Teff:
            self.str_Teff = "12000"
            print("Highest available Teff for PHOENIX spectra is 12000 K.")

        if debug:
            print(f"Matching Teff = {self.str_Teff}")

    def format_logg(self, debug=False):
        '''Create formatted logg string for filename for PHOENIX flux file.
        debug=True to print formatted string.
        '''
        min_logg = 0.00
        max_logg = 6.00
        logg = self.logg

        if logg < min_logg:
            self.str_logg = "0.00"
            print("Lowest available log(g) for PHOENIX spectra is 0.00.")
        elif min_logg <= logg <= max_logg:
            # Rounding to nearest 0.5, displaying 2 decimal places
            rounded = round(2.0 * logg) / 2.0
            self.str_logg = f"{rounded:.2f}"
        elif(self.logg > max_logg):
            self.str_logg = "6.00"
            print("Highest available log(g) for PHOENIX spectra is 6.00.")
        
        if debug:
            print(f"Matching log(g) = {self.str_logg}")

    def read_wavefile(self):    
        '''Tries to read the wavelength file for the PHOENIX spectrum.
        Loads only the 'optical' part of the wavelength grid (3500-9500 A)
        to save memory.
        '''
        # Load only optical region of wavelength grid, report possible errors
        print(f"** Trying to read {self.wave_uri}")
        try:
            with fits.open(self.wave_uri, use_fsspec=True, fsspec_kwargs={"anon": True}) as wave_file:
                self.wave_header = wave_file[0].header
                # Select wavelength range between 3000 A and 10000 A
                self.optical = np.logical_and(wave_file[0].data >= 4500, wave_file[0].data <= 7000)
                # Access wavelength data in the primary HDU
                self.wavelengths = wave_file[0].data[self.optical]
                print(f"   Wavelength file read successfully.")
        except FileNotFoundError:
            print(f"Error: File not found at {self.wave_uri}")
        except Exception as e:
            print(f"An error occurred while reading FITS file {self.wave_uri}: {e}")

    def read_fluxfile(self):
        '''Tries to read the chosen PHOENIX spectrum file matching Teff and logg.
        TO DO: add functionality to read other metallicities.
        Uses .format_Teff() and .format_logg() to find best matches for Teff and logg.
        TO DO: add filters for low logg values not available for hotter Teff's
        '''
        # Define flux grid file FTP location for specific Teff, logg, FeH
        self.format_Teff()
        self.format_logg()
        flux_file = "lte" + self.str_Teff + "-" + self.str_logg + "-" + self.str_FeH + self.flux_file_sfx
        self.flux_uri = self.ftp_server + self.flux_path + flux_file
        print(f"** Trying to read {self.flux_uri}")
        try:
            with fits.open(self.flux_uri, use_fsspec=True, fsspec_kwargs={"anon":True}) as spec_file:
                self.data_header = spec_file[0].header
                self.fluxes = spec_file[0].data[self.optical]
                self.Teff = self.data_header['PHXTEFF']
                self.logg = self.data_header['PHXLOGG']
                self.FeH  = self.data_header['PHXM_H']
                print(f"   Flux file read successfully: Teff = {self.Teff}, log(g) = {self.logg:.1f}, [Fe/H] = {self.FeH:.1f}")
        except FileNotFoundError:
            print(f"Error: File not found at {self.flux_uri}")
        except Exception as e:
            print(f"An error occurred while reading FITS file {self.flux_uri}: {e}")

    def normalize_fluxes(self):
        '''Flatten continuum of PHOENIX spectrum to range of 0..1.
        '''
        wls = self.wavelengths
        flx = self.fluxes
        Teff = self.Teff
        # Create Planck spectrum B_lambda for Teff
        bb = BlackBody(temperature = Teff * u.K, scale = 1.0 * u.erg / (u.cm**2 * u.s * u.AA * u.sr))
        bb_model = bb(wls)
        # Normalize template spectrum by Planck spectrum
        flx_bb = flx / bb_model
        # Then divide template spectrum by maximum to scale 0 to 1
        norm_flx = flx_bb / max(flx_bb)
        # Identify peaks of normalized spectrum to fit spline to
        peaklist, _ = find_peaks(norm_flx, distance=30000)
        wlspks = wls[peaklist]
        flxpks = norm_flx[peaklist]
        # Fit spline to peaks, then evaluate on wavelength grid
        tck = splrep(wlspks, flxpks, k=3)
        spline_fit = splev(wls, tck)
        self.norm_fluxes = norm_flx / spline_fit
        print("** Model flux normalized.")

    def prepare_template(self):
        '''Single method to load and normalize PHOENIX template.
        '''
        self.read_wavefile()
        self.read_fluxfile()
        self.normalize_fluxes()
