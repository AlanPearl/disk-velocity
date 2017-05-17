#! /usr/bin/python2
"""

### subsampling.py
### Alan Pearl

Subsamples quasars at various rates, and calculates the average 
magnitude of each bin

Flags: (must come before keywords)
======
	-save[=outfile]
	-keep		-> Doesn't plot anything. Useful for using the 
					data in another script
	-show		-> Shows each histogram at different subsampling levels

Keywords: (in default order)
=========
	-qsofile	-> default: 'qsocoords.fits'
	-outfile	-> default: 'qsocoordsfix.fits'

begin code:
"""

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import pyfits
import numpy as np
from math import pi, sin, cos, sqrt, atan2
from alan import writefits, read_argv
from numpy import delete as np_delete, array as np_array
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

# Define the alan functions here

# Defining variables
infile = 'qagcoords_ltpm30.fits'
outfile = 'subsample_qag_ltpm30.fits'
rates = np.linspace(0,45,181)[1:]
overwrite = False
save = False
show = False
keep = False
corrname = 'pearl'
xtit = 'Sampling Interval (deg)'
ytit = r'Mean Magnitude of Proper Motion (mas yr$^{-1}$)'
xlim = '0,45'
ylim = '0,8'


keys = ['infile', 'outfile', 'corrname', 'xlim', 'ylim']
flagvars = {'save':'outfile'}
argdict = read_argv(sys.argv, keys, flagvars)
vars().update(argdict)

xlim = map(float, xlim.split(','))
ylim = map(float, ylim.split(','))

############################################################


cvals = []

for i in xrange(257): # Blue -> Red
	if i < 96:
		red = .0 + .75*float(i)/96.
		green = 0. + .25*float(i)/96.
		blue = 1. - .25*float(i)/96.
	else:
		red = .75 - .25*float(i-96)/160.
		green = .25 - .25*float(i-96)/160.
		blue = .75 - .75*float(i-96)/160.
	
	cvals.append([red, green, blue])
cmapa = mpl.colors.ListedColormap(cvals)
	

norma = mpl.colors.Normalize(vmin=0, vmax=10, clip=True)



# MAIN FUNCTION
############################################################

t1 = pyfits.open(infile)[1].data

ra = t1['ra']
dec = t1['dec']
dec0 = dec.min()
dec1 = dec.max()

if corrname == '' or corrname.lower() == 'none':
	pmra = t1['pmra']
	pmde = t1['pmde']
else:
	pmra = t1['pmra_' + corrname]
	pmde = t1['pmde_' + corrname]


pma = np.sqrt( pmra**2 + pmde**2 )


means = [pma.mean()]
means_err = [pma.std() / np.sqrt(pma.size-1)]
for rate in rates:
	
	# The last bin may extend beyond the data
	#bins_r = np.arange(0.,360.+rate,rate)
	#bins_d = np.arange(dec0,dec1+rate,rate)
	
	# Upper parts of the data may not be included
	bins_r = np.arange(0.,360.,rate)
	bins_d = np.arange(dec0,dec1,rate)
	
	bins = [bins_r,bins_d]
	ra_cent = (bins_r + .5*rate)[:-1]
	dec_cent = (bins_d + .5*rate)[:-1]
	dec_cent, ra_cent = np.meshgrid(dec_cent,ra_cent)
	
	sample = np.array([ra,dec]).T
	
	sums_r,_ = np.histogramdd(sample,weights=pmra,bins=bins)
	sums_d,_ = np.histogramdd(sample,weights=pmde,bins=bins)
	nums  ,_ = np.histogramdd(sample,bins=bins)
	
	ra_cent = ra_cent.flatten()
	dec_cent = dec_cent.flatten()
	sums_r = sums_r.flatten()
	sums_d = sums_d.flatten()
	nums = nums.flatten()
	good = nums != 0
	
	avgs_r = sums_r * 1
	avgs_d = sums_d * 1
	avgs_r[good] /= nums[good]
	avgs_d[good] /= nums[good]
	avg_mag = np.sqrt( avgs_r[good]**2 + avgs_d[good]**2 )
	
	if show:
		fig,ax = plt.subplots()
		ax.hist2d(ra_cent[good], dec_cent[good], weights=avg_mag, 
					bins=bins, cmin=1e-10, cmap=cmapa, norm=norma)
		plt.show()
	
	avg_avg = avg_mag.mean()
	avg_err = avg_mag.std() / np.sqrt(avg_mag.size-1)
	means.append(avg_avg)
	means_err.append(avg_err)
	
	#print '\tSubsampling rate:', rate, 'degrees'
	#print 'Average pm magnitude:', avg_avg

rates = np.concatenate([[0.], rates])

############################################################

# CURVE FITTING
############################################################

means = np.array(means)
means_err = np.array(means_err)
avg_val = means.mean()
#weights = means_err**-2
#avg_val = np.sum(weights*means) / np.sum(weights)

# Create a best fit line
#p0 = (1,-1,0) # Initial guesses for A, B, and C
#def fitfunc(x, A, B, C): # A*e^(Bx) + C
#	return A * np.exp(B*x) + C
#popt, pcov = curve_fit(fitfunc, rates, means, p0, means_err, absolute_sigma=True)
#A, B, C = popt
#bfline_x = np.linspace(0,45,100)
#bfline_y = fitfunc( bfline_x, A, B, C)

if not keep:
	fig, ax = plt.subplots()
	#ax.plot(bfline_x, bfline_y, linestyle='--', linewidth=1, color='black')
	ax.errorbar(rates, means, yerr=means_err, linestyle='none')
	ax.scatter(rates, means)
	ax.text(30, 7, 'Avg: %.3f' %avg_val)
	ax.set_xlabel(xtit)
	ax.set_ylabel(ytit)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)

	if save:
		plt.savefig(outfile,bbox_inches='tight')
	else:
		plt.show()




