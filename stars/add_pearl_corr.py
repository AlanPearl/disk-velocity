#! /usr/bin/python2
"""

### add_pearl_corr.py
### Alan Pearl

Uses the table pearl_corr.fits to implement the pearl correction.

Flags: (must begin with a dash, followed by a non-numeral character)
======
	-overwrite	-> Automatically overwrite during all file conflicts

Keywords: (in default order)
=========
	-infile		-> default: carlin_star_catalog.fits
	-outfile	-> default: carlin_star_catalog_pearl.fits
	-corrfile	-> default: pearl_corr.fits
	-shape		-> default: 1440,360 (number of bins in ra,dec)
	-ralim		-> default: 0,360
	-declim		-> default: -15,75

"""

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import pyfits as fits
from alan import writefits, read_argv
import numpy as np

def pearl_corr(stars,pearlcorr,shape,ralim,declim):
	rabinnum, decbinnum = shape
	ra = stars['ra']
	dec = stars['dec']
	pmra = stars['pmra']
	pmde = stars['pmde']
	
	pmra_shift = pearlcorr['pmra_shift']
	pmde_shift = pearlcorr['pmde_shift']
	pmra_se = pearlcorr['pmra_se']
	pmde_se = pearlcorr['pmde_se']
	
	badshift = np.array([0.]).astype(pmra_shift.dtype)
	badse = np.array([100.]).astype(pmra_se.dtype)

	pmra_shift = np.concatenate((pmra_shift,badshift))
	pmde_shift = np.concatenate((pmde_shift,badshift))
	pmra_se = np.concatenate((pmra_se,badse))
	pmde_se = np.concatenate((pmde_se,badse))


	ra_int = float(ralim[1] - ralim[0]) / rabinnum
	dec_int = float(declim[1] - declim[0]) / decbinnum

	ra_index = ((ra - ralim[0]) // ra_int).astype(int)
	dec_index = ((dec - declim[0]) // dec_int).astype(int)

	index = dec_index*rabinnum + ra_index
	out_ra = np.where( np.logical_or( 0 > ra_index , ra_index > rabinnum ) )
	out_dec = np.where( np.logical_or( 0 > dec_index , dec_index > rabinnum ) )
	index[out_ra] = -1
	index[out_dec] = -1

	pmra_shift = pmra_shift[index]
	pmde_shift = pmde_shift[index]
	pmra_se = pmra_se[index]
	pmde_se = pmde_se[index]

	pmra_pearl = pmra + pmra_shift
	pmde_pearl = pmde + pmde_shift
	
	return pmra_pearl, pmde_pearl, pmra_se, pmde_se


if __name__ == '__main__':
	
	infile = 'LAMOST.fits'
	corrfile = 'pearl_corr.fits'
	outfile = 'LAMOST.fits'
	overwrite = False
	# Default:
	# 360 dec bins over range (-15,75); 1440 ra bins over range (0,360)
	shape = '1440,360'
	ralim = '0,360'
	declim = '-15,75'
	
	argdict = read_argv(sys.argv)
	vars().update(argdict)
	shape = shape.split(',')
	shape[0] = int(shape[0])
	shape[1] = int(shape[1])
	ralim = ralim.split(',')
	ralim[0] = float(ralim[0])
	ralim[1] = float(ralim[1])
	declim = declim.split(',')
	declim[0] = float(declim[0])
	declim[1] = float(declim[1])
	
	print 'Reading files...'
	stars = fits.open(infile)[1].data
	pearlcorr = fits.open(corrfile)[1].data
	
	print 'Calculating corrected proper motions...'
	pmra_pearl,pmde_pearl,pmra_se,pmde_se = pearl_corr(stars, pearlcorr, shape, ralim, declim)
	
	outdic = {}
	for name in stars.names:
		outdic[name] = stars[name]
	outdic['pmra_pearl'] = pmra_pearl
	outdic['pmde_pearl'] = pmde_pearl
	outdic['pmra_se'] = pmra_se
	outdic['pmde_se'] = pmde_se
	
	print 'Writing to %s...' %outfile
	writefits(outfile, outdic, overwrite=overwrite)