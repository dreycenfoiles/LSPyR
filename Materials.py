import numpy as np
from scipy.interpolate import interp1d

def Gold(wavelength):

	energy = np.array([0.64, 0.77, 0.89, 1.02, 1.14, 1.26, 1.39, 1.51, 
    1.64, 1.76, 1.88, 2.01, 2.13, 2.26, 2.38, 2.50, 2.63, 2.75, 2.88, 
    3.00, 3.12, 3.25, 3.37, 3.50, 3.62, 3.74, 3.87, 3.99, 4.12, 4.24, 
    4.36, 4.49, 4.61, 4.74, 4.86, 4.98, 5.11, 5.23, 5.36, 5.48, 5.60, 
    5.73, 5.85, 5.98, 6.10, 6.22, 6.35, 6.47, 6.60])

	wavelengths = 12.3984/energy 

	refractive_indices = np.array([0.92 + 13.78j, 0.56 + 11.21j, 0.43 + 9.519j, 
	0.35 + 8.145j, 0.27 + 7.150j, 0.22 + 6.350j, 0.17 + 5.663j, 
	0.16 + 5.083j, 0.14 + 4.542j, 0.13 + 4.103j, 0.14 + 3.697j, 
	0.21 + 3.272j, 0.29 + 2.863j, 0.43 + 2.455j, 0.62 + 2.081j, 
	1.04 + 1.833j, 1.31 + 1.849j, 1.38 + 1.914j, 1.45 + 1.948j, 
	1.46 + 1.958j, 1.47 + 1.952j, 1.46 + 1.933j, 1.48 + 1.895j, 
	1.50 + 1.866j, 1.48 + 1.871j, 1.48 + 1.883j, 1.54 + 1.898j, 
	1.53 + 1.893j, 1.53 + 1.889j, 1.49 + 1.878j, 1.47 + 1.869j, 
	1.43 + 1.847j, 1.38 + 1.803j, 1.35 + 1.749j, 1.33 + 1.688j, 
	1.33 + 1.631j, 1.32 + 1.577j, 1.32 + 1.536j, 1.30 + 1.497j, 
	1.31 + 1.460j, 1.30 + 1.427j, 1.30 + 1.387j, 1.30 + 1.350j, 
	1.30 + 1.304j, 1.33 + 1.277j, 1.33 + 1.251j, 1.34 + 1.226j, 
	1.32 + 1.203j, 1.28 + 1.188j])

	permittivity = refractive_indices**2
	real_permittivity = permittivity.real
	imag_permittivity = permittivity.imag 
	
	real_fit = interp1d(wavelengths,real_permittivity)
	imag_fit = interp1d(wavelengths,imag_permittivity)

	return np.complex(real_fit(wavelength),imag_fit(wavelength))



