# class NRESObservation

Represents a single NRES spectral observation

Methods to load and process Las Cumbres Oobservatory NRES spectrum

NRES FITS file contains fluxes & uncertainties, wavelengths, and masks
for each echelle order.

Fluxes are both unnormalized and normalized for blaze angle.
Normalized orders are preferred for cross-correlation and other tasks.

## Attributes:

* filepath (str) -- directory path to NRES file

* filename (str) -- name of NRES file to be loaded

* primary_header (FITS header) -- main header for observation

* spectrum_header (FITS header) -- data header for echelle orders

* orders () -- integer for each echelle order

* wavelengths () -- wavelengths for each echelle order 

* flux_mask () -- Boolean masks of ??? (TO-DO: ask LCO about this.)

* fluxes () -- fluxes normalized for blaze angle for each echelle order

* flux_uncertainty () -- formal uncertainty of fluxes for each order

* obstime (Astropy Time object) -- midpoint of exposure at telescope

* skycoords (SkyCoord object) -- coordinates of observed star

* location (EarthLocation object) -- description of telescope site

* bjd_tdb (Astropy Time object) -- midpoint of exposure corrected for
    light-travel-time to Solar System barycenter

* rvbarycorr (Astropy velocity) -- radial velocity correction between
    telescope site on Earth and Solar System barycenter 

## Methods:

* identify_fiber() -- Read header to locate fiber for sky object.<br>
    returns integer 0 or 2

* read_fluxfile() -- Read NRES spectrum file.<br>
    returns primary_header, spectrum_header, orders,
            wavelengths, flux_mask, fluxes, flux_uncertainty

* get_earthlocation() -- Identify telescope site, create EarthLocation.<br>
    returns location

* get_observingparams() -- Calculate exposure midpoint, assign SkyCoord<br>
    returns obstime, skycoords

* get_barycorrections() -- Calculate light-travel-time, bary RV correction<br>
    returns bjd_tdb, rvbarycorr

* join_orders(first_order, last_order) -- Join echelle orders into spectrum<br>
    returns joined_wave, joined_flux, joined_func