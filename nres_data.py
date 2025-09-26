
'''A class to represent a single NRES spectral observation.

Loads and processes an LCO NRES spectrum from its FITS file.

NRES FITS file contains fluxes & uncertainties, wavelengths, and masks
for each echelle order.

Fluxes are both unnormalized and normalized for blaze angle.
Normalized orders are preferred for cross-correlation and other tasks.

Attributes
----------
filepath (str) -- directory path to NRES file

filename (str) -- name of NRES file to be loaded

primary_header (FITS header) -- main header for observation

spectrum_header (FITS header) -- data header for echelle orders

orders () -- integer for each echelle order

wavelengths () -- wavelengths for each echelle order 

flux_mask () -- Boolean masks of ??? (TO-DO: ask LCO about this.)

fluxes () -- fluxes normalized for blaze angle for each echelle order

flux_uncertainty () -- formal uncertainty of fluxes for each order

obstime (Astropy Time object) -- midpoint of exposure at telescope

skycoords (SkyCoord object) -- coordinates of observed star

location (EarthLocation object) -- description of telescope site

bjd_tdb (Astropy Time object) -- midpoint of exposure corrected for
    light-travel-time to Solar System barycenter

rvbarycorr (Astropy velocity) -- radial velocity correction between
    telescope site on Earth and Solar System barycenter 


Methods
-------
identify_fiber() -- Read header to locate fiber for sky object.
    returns integer 0 or 2

read_fluxfile() -- Read NRES spectrum file.
    returns primary_header, spectrum_header, orders,
            wavelengths, flux_mask, fluxes, flux_uncertainty

get_earthlocation() -- Identify telescope site, create EarthLocation.
    returns location

get_observingparams() -- Calculate exposure midpoint, assign SkyCoord
    returns obstime, skycoords

get_barycorrections() -- Calculate light-travel-time, bary RV correction
    returns bjd_tdb, rvbarycorr

join_orders(first_order, last_order) -- Join echelle orders into spectrum.
    returns joined_wave, joined_flux, joined_func

*** Author: Kevin Healy, Mesa CC
'''

from astropy.io import fits
from astropy import units
from astropy import constants
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation
import numpy as np

class NRESObservation:
    '''Class for NRES data and derived properties.'''
    def __init__(self, filepath, filename):
        self.bjd_tdb = None
        self.filepath = filepath
        self.filename = filename
        self.joined_wave = None
        self.joined_flux = None
        self.joined_func = None
        self.primary_header = None
        self.spectrum_header = None
        self.orders = None
        self.flux_mask = None
        self.fluxes = None
        self.flux_uncertainty = None
        self.obssite = None
        self.obstime = None
        self.skycoords = None
        self.location = None
        self.rvbarycorr = None
        self.fitted_vel = None
        self.fitted_uvel = None
        self.radialvelocity = None

        self.wavelengths = None

    def identify_fiber(self):
        '''Identify fiber order and return fiber for sky target.
        '''
        fibers = self.primary_header['OBJECTS'].split('&')
        if fibers[0] == 'none':
            sky_fiber = 2 # fiber order is 'none&thar&STAR'
        elif fibers[2] == 'none':
            sky_fiber = 0 # fiber order is 'STAR&thar&none'
        else:
            sky_fiber = 9 # unknown fiber order
        return sky_fiber    

    def read_fluxfile(self):
        '''Try to read the chosen NRES spectrum file.'''
        try:
            with fits.open(self.filepath + self.filename) as spec_file:
                # Retrieve the headers from the first two HDUs
                self.primary_header  = spec_file['PRIMARY'].header
                self.spectrum_header = spec_file['SPECTRUM'].header
                
                # Retrieve the data for the science target
                data = spec_file['SPECTRUM'].data
                sky_fiber = self.identify_fiber()
                if sky_fiber != 9:
                    data = data[data['fiber'] == sky_fiber]
                    print(f"NRES fiber for sky is {sky_fiber}.")
                else:
                    print("Unknown NRES fiber order! Empty spectra will result.")

                # Retrieve array of echelle orders
                self.orders = data['order']

                # Retrieve arrays for wavelengths and blaze-normalized fluxes, 68x4096
                self.wavelengths = data['wavelength']
                self.flux_mask   = data['mask']
                self.fluxes      = data['normflux']
                self.flux_uncertainty = data['normuncertainty']
                print(f"Reading spectrum file {self.filename}.")
        except FileNotFoundError:
            print(f"Error: File not found at {self.filepath}{self.filename}")
        except Exception as e:
            print(f"An error occurred while reading FITS file {self.filename}: {e}")

    def get_earthlocation(self):
        '''Try to identify telescope site and create EarthLocation'''
        # Retrieve observatory site from FITS header
        obssite = self.primary_header['SITE']
        # Try to match the 'SITE' description in the header.
        # As of June 2025, there are 4 NRES 1-meter telescopes
        # Cerro Tololo, McDonald,  SAAO,  Wise
        #    'ctio'    'mcdonald' 'SAAO' 'wise'
        if ('Cerro' in obssite):
            obssite = 'ctio'
            self.location = EarthLocation.of_site(obssite)
        elif ('McDonald' in obssite):
            obssite = 'mcdonald'
            self.location = EarthLocation.of_site(obssite)
        elif ('SAAO' in obssite):
            obssite = 'SAAO'
            self.location = EarthLocation.of_site(obssite)
        elif ('Wise' in obssite):
            obssite = 'wise'
            self.location = EarthLocation.of_site(obssite)
        else:
            print("Value of 'SITE' in FITS header is not recognized.")
        self.obssite = obssite

    def get_observingparams(self):
        '''Calculate midpoint time of exposure,
        assign SkyCoord object for observing target
        '''
        # Retrieve length and date/time of exposure from FITS header,
        # then calculate: 'obstime' = midpoint of exposure
        exptime = self.primary_header['EXPTIME']
        dateobs = self.primary_header['DATE-OBS']
        self.obstime = Time(dateobs) + 0.5 * exptime * units.second

        # Retrieve RA/DEC from FITS header, convert to Astropy SkyCoord object
        ra  = self.primary_header['RA']
        dec = self.primary_header['Dec']
        ra_coord  = Angle(ra, unit='hourangle').deg
        dec_coord = Angle(dec, unit='degree').deg
        self.skycoords = SkyCoord(ra=ra_coord*units.deg, dec=dec_coord*units.deg, frame='icrs')

    def get_barycorrections(self):
        '''Calculate barycentric time and velocity correction.
        'bjd_tdb' = obs time corrected for light-travel to barycenter
        'rvbarycorr' = obs vel corrected for Earth motion wrt barycenter
        '''
        self.bjd_tdb = self.obstime.tdb + self.obstime.light_travel_time(self.skycoords, kind='barycentric', location=self.location)
        print(f"Barycenter time of observation is {self.bjd_tdb}")
        self.rvbarycorr = self.skycoords.radial_velocity_correction(kind='barycentric', location=self.location, obstime=self.obstime)
        self.rvbarycorr = self.rvbarycorr.to('km/s')
        print(f"Barycentric velocity correction is {self.rvbarycorr.value:.3f} km/s.")

    def get_radialvelocity(self):
        '''Calculate barycentric-corrected radial velocity.
        Offset velocity found from spectrum fitting
        by barycentric correction from `get_barycorrections()`
        Form of correction described in documentation:
        https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.radial_velocity_correction
        'radialvelocity' = corrected velocity of star
        '''
        c_kms = constants.c.to('km/s').value
        v_m  = self.fitted_vel
        v_bc = self.rvbarycorr.value
        self.radialvelocity = v_m + v_bc + (v_m * v_bc / c_kms)

    def join_orders(self, first_order, last_order):
        '''Join wavelength, flux, uncertainty arrays for echelle orders.
        Trim echelle spectrum to non-zero wavelengths and fluxes, then
        clip (normalized) fluxes greater than 2, trim overlap in 
        wavelength by selecting smaller uncertainty, joining
        wavelengths into single array, fluxes into single array,
        flux uncertainties into single array.
        '''
        # Identify indices of chosen first and last echelle orders
        low_index  = np.argwhere(self.orders == first_order)[0][0]
        high_index = np.argwhere(self.orders ==  last_order)[0][0]

        # Start join w/ longest-wavelength echelle order (lowest index)...
        idx = low_index
        valid = np.logical_and((self.fluxes[idx] != 0),(self.fluxes[idx] < 1.75))
        order_wave = self.wavelengths[idx]
        order_flux = self.fluxes[idx].astype('float16')
        order_func = self.flux_uncertainty[idx].astype('float16')
        all_wave = order_wave[valid]
        all_flux = order_flux[valid]
        all_func = order_func[valid]
        min_wave = min(order_wave[valid])

        # Then append other wavelengths in index order,
        # trimming off short-wavelength overlap of short-wavelength order
        for idx in range(low_index+1, high_index):
            min_wave = min(all_wave)
            valid = np.logical_and((self.fluxes[idx] != 0),(self.fluxes[idx] < 1.75))
            order_wave = self.wavelengths[idx]
            order_flux = self.fluxes[idx].astype('float16')
            order_func = self.flux_uncertainty[idx].astype('float16')
            valid_obs_wave = order_wave[valid]
            valid_obs_flux = order_flux[valid]
            valid_obs_func = order_func[valid]
            trim = (valid_obs_wave <= min_wave)
            all_wave = np.concatenate((valid_obs_wave[trim], all_wave))
            all_flux = np.concatenate((valid_obs_flux[trim], all_flux))
            all_func = np.concatenate((valid_obs_func[trim], all_func))
        self.joined_wave = all_wave
        self.joined_flux = all_flux
        self.joined_func = all_func
    
    def clean_spectrum(self, sigma=3.0):
        '''Remove elements of flux array significantly different from local mean.
        Need to identify and delete flux values that are far above or below the
        mean of the neighboring values. Something like sigma-clipping, with an
        awareness of how 
        '''