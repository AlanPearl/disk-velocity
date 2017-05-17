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
import recenter # (vickers correction)

def vick_corr(pmra,pmde,ra,dec,J):
	"""
	Returns 2-tuple containing numpy arrays of the Vickers-corrected proper
	motion values (ra and dec, respectively) given the uncorrected values from
	the PPMXL database, their position (ra,dec), and their J magnitude
	"""
	J = list(J)
	
	vickra = [np.float64(0.)] * R.size
	vickde = [np.float64(0.)] * R.size
	for i in xrange(R.size):
		if i%50000==0:
			print "%d\t/ %d Vickers' corrections calculated"%(i,R.size)
		answer = recenter.recenter(ra[i], dec[i], J[i])
		if answer == None:
			vickra[i], vickde[i] = np.nan,np.nan
		else:
			vickra[i], vickde[i] = answer
	vickra = np.array(vickra)
	vickde = np.array(vickde)
	pmra_vick = pmra + vickra
	pmde_vick = pmde + vickde
	
	return pmra_vick, pmde_vick

if __name__ == '__main__':
	infile1 = '../data/LAMOST.fits'
	infile2 = '../data/RAVE_DR5.fits'
	outfile = '../data/LAMOST-RAVE.fits'
	vars().update(read_argv(sys.argv))
	
	L = pyfits.open(infile1)[1].data
	R = pyfits.open(infile2)[1].data
	
	ra = np.concatenate([ L['ra'], R['RAdeg'] ])
	dec = np.concatenate([ L['dec'], R['DEdeg']])
	
	#NOTE: converting pc to kpc for RAVE distances (these are not the Carlin distance calculation, but simply the calculations that came with the RAVE data release)
	dist = np.concatenate([ L['dist'], R['distance']/1000.])
	dist_err = np.concatenate([ L['dist_err'], R['edistance']/1000.])
	
	rv = np.concatenate([ L['rv'], R['HRV'] ])
	rv_err = np.concatenate([ L['rv_err'], R['eHRV'] ])
	
	
	pmra = np.concatenate([ L['pmra'], R['pmRA_PPMXL'] ])
	pmde = np.concatenate([ L['pmde'], R['pmDE_PPMXL'] ])
	epma = np.concatenate([ L['epma'], R['epmRA_PPMXL'] ])
	epmd = np.concatenate([ L['epmd'], R['epmDE_PPMXL'] ])
	
	# No Pearl proper motion correction for RAVE stars :(
	# Assign them NaN values and errors of 100
	pmra_pearl = np.concatenate([ L['pmra_pearl'], [np.nan]*R.size ])
	pmde_pearl = np.concatenate([ L['pmde_pearl'], [np.nan]*R.size ])
	pmra_se = np.concatenate([ L['pmra_se'], [100.0]*R.size ])
	pmde_se = np.concatenate([ L['pmde_se'], [100.0]*R.size ])
	
	# Calculate Vickers-corrected proper motions for RAVE stars
	pmra_vickR, pmde_vickR = vick_corr(R['pmRA_PPMXL'], R['pmDE_PPMXL'], 
					R['RAdeg'], R['DEdeg'], R['Jmag_2MASS'])
	pmra_vick = np.concatenate([ L['pmra_vick'], pmra_vickR ])
	pmde_vick = np.concatenate([ L['pmde_vick'], pmde_vickR ])
	
	snr = np.concatenate([ L['snr'], R['SNR_K'] ])
	
	# No Carlin Magnitude calculations for RAVE stars
	M_K50 = np.concatenate([ L['M_K50'], [np.nan]*R.size ])
	# These values don't really matter if we don't have M_K (no RAVE data for HR diagram)
	k0 = np.concatenate([ L['k0'], [np.nan]*R.size ])
	j0 = np.concatenate([ L['j0'], [np.nan]*R.size ])
	
	# RAVE is not a spectroscopic survey, so they have no subclass (star type)
	# But I'll use this to check which stars came from RAVE
	subclass = np.concatenate([ L['subclass'], ['RAVE']*R.size ])
	
	
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
	'pmra_pearl', 
	'pmde_pearl', 
	'pmra_se', 
	'pmde_se', 
	'snr', 
	'M_K50', 
	'k0', 
	'j0', 
	'subclass'
	]
	
	outdict = {}
	for name in nameorder[:-1]:
		outdict[name] = vars()[name].astype(np.float32)
	outdict[nameorder[-1]] = vars()[nameorder[-1]]
	
	print 'Writing to %s' %outfile
	writefits(outfile, outdict, nameorder)