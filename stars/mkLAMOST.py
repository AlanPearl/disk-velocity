#! /usr/bin/python2

import os,sys
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import pyfits
import numpy as np
from alan import writefits, read_argv
import recenter

def distCarlin(distk50, distr50, distkerr, distrerr):
	"""Returns 2-tuple containing numpy arrays of distance and distance error,
	taking the error weighted mean of distk50 and distr50 from
	Jeff Carlin's Bayesian distance estimates.
	"""
	r_helio = (distk50/distkerr**2 + distr50/distrerr**2)
	r_helio_err = 1. / (distkerr**-2 + distrerr**-2)
	r_helio *= r_helio_err
	r_helio_err = np.sqrt(r_helio_err)

	ob = np.logical_or( 
		np.logical_not( r_helio_err < np.inf ) , np.isnan( r_helio ) )

	use_k = np.logical_or( distkerr[ob] <= distrerr[ob] , 
			np.logical_not( distrerr[ob] < np.inf ) )

	r_helio[ob] = distr50[ob]
	r_helio_err[ob] = distrerr[ob]

	r_helio[ob][use_k] = distk50[ob][use_k]
	r_helio_err[ob][use_k] = distkerr[ob][use_k]
	
	return r_helio, r_helio_err


if __name__ == '__main__':
	infile = 'dr3_dr4q1q2_stellar_gr0_ebv_2MASS-dist_shorter_dupecombine_PPMXL_UCAC4_1p5arcsec_vickers_pmcorr.fits'
	outfile = 'LAMOST.fits'
	vars().update(read_argv(sys.argv))
	
	L = pyfits.open(infile)[1].data
	
	outdict = {}
	# Use the error weighted mean of the distance calculations
	outdict['dist'], outdict['dist_err'] = distCarlin(
						L['distk50'],L['distr50'],L['distkerr'],L['distrerr'])
	outdict['pmra_vick'] = L['pmra_vickers_corr']
	outdict['pmde_vick'] = L['pmde_vickers_corr']
	# The signal to noise ratio in each filter must be >5, so we only need
	# the minimum of the three values to determine if we can use the star
	outdict['snr'] = np.min([L['snrr'],L['snrg'],L['snri']], axis=0)
	

	nameorder = [
	'ra', 
	'dec', 
	'dist', 
	'dist_err', 
	'rv', 
	'rv_err', 
	'pmra', 
	'pmde', 
	'epma', 
	'epmd', 
	'pmra_vick', 
	'pmde_vick', 
	'snr', 
	'M_K50', 
	'k0', 
	'j0', 
	'subclass'
	]
	
	for name in nameorder:
		if not name in outdict:
			outdict[name] = L[name]

	print 'Writing to %s' %outfile
	writefits(outfile, outdict, nameorder)