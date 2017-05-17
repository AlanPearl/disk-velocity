import numpy as np
import os
from numpy.lib import recfunctions as rfs
from scipy.special import sph_harm


rfile = os.path.join(os.path.dirname(__file__), 'pmr.csv')
dfile = os.path.join(os.path.dirname(__file__), 'pmd.csv')

#fits
ra_fit = np.genfromtxt(rfile, delimiter=',')
de_fit = np.genfromtxt(dfile, delimiter=',')

#These files are 7 lines long, the first column is the magnitude slice it was
#fit to. After that the columns are leading coefficients for:
#sph(order, degree) = (0,0), (0,1), (1,1)_real, (1,1)_imag, (0,2), (1,2)_real,
# (1,2)_imag, (2,2)_real, (2,2)_imag, (0,3)...
#we do not use negative orders as they add no new information
#the zeroth orders never have imaginary information

harm_degree = 8






def surface_harmonics(point, popt):
	"""returns raw corrections to ra, dec of point=(ra, dec) and the
	sph. harmonics coefficients popt.

	popt is in the order discussed above.
	"""
	#convert data to sph_harmable, work in radians with colatitude
	deg2rad = np.pi / 180.

	longitude, latitude = point

	colatitude = (90-latitude)*deg2rad
	longitude = longitude*deg2rad

	coeff_count = 0

	# the result is the sum of all the sph harmonics of the soln in this spot
	res = 0
	for degree in range(harm_degree + 1):
		#we only use the positive orders because the negatives are just inverses
		#of the positives
		for order in range(degree + 1):

			#Find intrinsic value of spherical harmonic
			plot_harm = sph_harm(order, degree, longitude, colatitude)

			#And apply relevant fit coefficient
			if order == 0:
				res += popt[coeff_count] * plot_harm.real
				coeff_count += 1
			else:
				res += popt[coeff_count] * plot_harm.real
				coeff_count += 1
				res += popt[coeff_count] * plot_harm.imag
				coeff_count += 1

	#Return estimated offset
	return( res )




def recenter(ra, de, jmag):
	"""returns the corrections to PPMXL pmd and pma based on a position
	in degrees and a J magnitude (in mag).

	The corrections are in mas/yr.
	"""
	#the data used in the fits are in shells 0.1 mags wide and half a mag
	#apart [14.0-14.1, 14.5-14.6, ...]
	#that's the basis for these cuts and stuff


	#fringe cases
	if jmag < 14.05: #assume first fit
		#find sum of all orders and degrees of sph harmonics at given point
		#multiplied by the fit coefficients for each order and degree

		#Note fit[shell][zeroth column] is magnitude so we cut it out
		pmr_corr = surface_harmonics([ra, de], ra_fit[0][1:])
		pmd_corr = surface_harmonics([ra, de], de_fit[0][1:])
		return [pmr_corr, pmd_corr]


	if jmag > 17.05: #assume last fit
		pmr_corr = surface_harmonics([ra, de], ra_fit[-1][1:])
		pmd_corr = surface_harmonics([ra, de], de_fit[-1][1:])
		return [pmr_corr, pmd_corr]





	for i in range(len(ra_fit)):
		#find the shells which enclose this particular j
		if ra_fit[i][0]+0.05 >= jmag:

			#linearly interpolate between the prior shell fit and the next
			#based on the j magnitude
			travel_frac = (jmag - (ra_fit[i-1][0] + 0.05)) / 0.5

			pmr_corr_up = surface_harmonics([ra, de], ra_fit[i][1:])
			pmr_corr_down = surface_harmonics([ra, de], ra_fit[i-1][1:])

			pmd_corr_up = surface_harmonics([ra, de], de_fit[i][1:])
			pmd_corr_down = surface_harmonics([ra, de], de_fit[i-1][1:])

			pmr_corr = pmr_corr_down + travel_frac*( pmr_corr_up - pmr_corr_down )
			pmd_corr = pmd_corr_down + travel_frac*( pmd_corr_up - pmd_corr_down )

			return( pmr_corr, pmd_corr )




def _test():

	#fits
	ra_fit = np.genfromtxt('pmr.csv', delimiter=',')
	de_fit = np.genfromtxt('pmd.csv', delimiter=',')

	data_pre = np.genfromtxt('vc_qso.csv', delimiter=',', names=True)

	pmr_new, pmd_new = [], []
	for star in data_pre:
		#correction for each
		pmr_corr, pmd_corr = recenter(
			star['ra'], star['de'], star['jmag']
		)

		pmr_new.append(star['pmr_mas'] - pmr_corr)
		pmd_new.append(star['pmd_mas'] - pmd_corr)

	data_post = rfs.append_fields(
		data_pre, names=['pmr_corr_mas', 'pmd_corr_mas'],
		data=[pmr_new, pmd_new], dtypes=['>f4','>f4'], usemask=False
	)

	np.savetxt(
		'vc_qso_recentered.csv', data_post, delimiter=',', comments='',
		header='pmra,pmde,jmag,ra,de,pmr_mas,pmd_mas,pmr_corr_mas,pmd_corr_mas'
	)


if __name__ == '__main__':

	_test()
