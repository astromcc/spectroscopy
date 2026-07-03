'''Performs several operations on stellar spectra and template spectra.

Assumes NRES spectrum as class NRESObservation in nres_load.py
and PHOENIX spectrum as class PhoenixTemplate in phoenix_load.py

Processing mimics steps in rv.py module in banzai-nres
https://github.com/LCOGT/banzai-nres/tree/main/banzai_nres

*** Author: Kevin Healy, Mesa CC
'''
import numpy as np

from astropy import units
from astropy import constants

from scipy.optimize import curve_fit

def interp_template(obs_wave, template_obj):
    '''Trim to observed wavelength range, resample template fluxes.
    '''
    # Identify shortest and longest wavelength in observed spectrum
    wave_min, wave_max = min(obs_wave), max(obs_wave)
    # Mask template spectrum to observed wavelength range
    trim = np.logical_and((template_obj.wavelengths >= wave_min - 1.0),
                          (template_obj.wavelengths <= wave_max + 1.0))
    trimmed_wave = template_obj.wavelengths[trim]
    trimmed_flux = template_obj.norm_fluxes[trim]
    # Then interpolate template fluxes to Doppler-shifted wavelengths
    flux_interp = np.interp(obs_wave, trimmed_wave, trimmed_flux)
    return flux_interp

def doppler(wavelengths, v_kms):
    '''Shift wavelength array by small-v approximation Doppler shift.
    Positive velocity means shifted wavelengths are stretched.
    Velocity must be expressed in km/sec for unit consistency.
    '''
    c_kms = constants.c.to('km/s').value
    shifted_waves = wavelengths * (1.0 + v_kms / c_kms)
    return shifted_waves

def gaussian(x, A, mu, sigma, b):
    '''Gaussian function to fit peak of inner-product array.
    '''
    f = A*np.exp(-(x-mu)**2/(2*sigma**2)) + b
    return f

def velocity_scan(velocity_range, nres_obj, template_obj):
    '''Calculate inner product of two spectra over velocity range.
    This is a good measure of the fit between two spectra.
    Also calculate formal uncertainties of inner product from
    propogation of uncertainty of observed fluxes assuming zero
    uncertainty in PHOENIX spectrum. Uncertainty of inner product
    is sq. root of sum of squares of normalized flux uncertainties
    = sqrt(sum(squares(obs_func/obs_flux)))
    '''
    innerprods = np.empty(len(velocity_range))
    unc_innerprods = np.empty(len(velocity_range))
    for idx in range(len(velocity_range)):
        # Doppler shift *observed* wavelengths *to* template wavelengths
        # so Doppler velocity is opposite sign of search velocity.
        velocity = velocity_range[idx]
        new_wave = doppler(nres_obj.joined_wave, -velocity)
        flux_interp = interp_template(new_wave, template_obj)
        innerprods[idx] = np.inner(nres_obj.joined_flux, flux_interp)
        unc_innerprods[idx] = np.sqrt(np.sum((nres_obj.joined_func/nres_obj.joined_flux)**2))
        norm_innprod = innerprods / max(innerprods)
        norm_innunc  = unc_innerprods / max(innerprods)
    return norm_innprod, norm_innunc

def find_velocity(nres_obj, template_obj):
    '''Search for velocity that best matches two spectra.
    Search is performed twice:
    * Coarse velocity range = -1000 to +1000 km/s, dv = 20 km/s
    * Maximum of inner products from coarse search is identified, then
    * Fine search in window around max of coarse search, dv = fine_dv
    '''
    # Width and spacing of coarse search
    coarse_min = -1000. # km/s
    coarse_max = +1000. # km/s
    coarse_dv = 20. # km/s
    medium_dv = 1.0 # km/s
    fine_dv   = 0.1 # km/s
    
    # Coarse search over wide velocity range, ignore scan uncertainty
    coarse_vels = np.arange(coarse_min, coarse_max, coarse_dv)
    coarse_scan, _ = velocity_scan(coarse_vels, nres_obj, template_obj)

    # Identify largest inner product => best estimate of velocity
    coarse_peak_vel = coarse_vels[np.argmax(coarse_scan)]
    print(f"Coarse search done, best-fit velocity is {coarse_peak_vel:+.1f} km/s.")

    # Medium search around peak of coarse search, ignore scan uncertainty
    medium_vels = np.arange(coarse_peak_vel - 10., coarse_peak_vel + 10., medium_dv)
    medium_scan, _ = velocity_scan(medium_vels, nres_obj, template_obj)

    # Identify largest inner product => best estimate of velocity
    medium_peak_vel = medium_vels[np.argmax(medium_scan)]
    print(f"Medium search done, best-fit velocity is {medium_peak_vel:+.1f} km/s.")

    # Fine search around peak of medium search, narrow velocity spacing
    fine_vels = np.arange(medium_peak_vel - 10., medium_peak_vel + 10., fine_dv)
    fine_scan, fine_unc = velocity_scan(fine_vels, nres_obj, template_obj)

    # Fit fine_scan to determine best velocity and error on velocity
    # 'gaussian' = f(A, mu, sigma, b),
    #              where A = amplitude, b = baseline,
    #                    mu = mean value, sigma = width of curve
    # Best estimate of velocity is 2nd parameter (index = 1)
    # Best estimate of error in velocity is variance of parameter
    params, pcov = curve_fit(gaussian, fine_vels, fine_scan,
                             bounds=([0, fine_vels[0], 0, 0.9], [0.1, fine_vels[-1], 10, 1]),
                             sigma=fine_unc, absolute_sigma=True)
    velocity_est = params[1]
    velocity_unc = np.sqrt(pcov[1][1])

    print(f"Fine search done, best velocity estimate is {velocity_est:+.3f} km/s.")
    print(f"Estimate of velocity uncertainty is {velocity_unc:0.4f} km/s.")
    nres_obj.fitted_vel  = velocity_est
    nres_obj.fitted_uvel = velocity_unc

