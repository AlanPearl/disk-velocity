#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import matplotlib.pyplot as plt
import numpy as np
import pyfits as fits
from math import sqrt
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt

args = sys.argv
given = set([])

for arg in args[1:]:
	exec(arg)
	given.add(arg.split('=')[0])

if not 'plotnum' in given:
	plotnum = 1
if not 'save' in given:
	save = False
if not 'gridline' in given:
	gridline = False

if plotnum == 1:
	
	if not 'infile' in given:
		infile = '../data/qagcoords_ltpm30.fits'
		#infile = 'qsocoords_dedup.fits'
	if not 'xtit' in given:
		xtit = r'Proper Motion (mas yr$^{-1}$)'
	if not 'ytit' in given:
		ytit = 'Number of Objects'
	if not 'xmajtick' in given:
		xmajtick = 5.
	if not 'ymajtick' in given:
		ymajtick = 1.
	if not 'xlim' in given:
		xlim = [-15,15]
	if not 'ylim' in given:
		ylim = [0,3000]
	if not 'corrname' in given:
		corrname = 'none'
	if not 'binperx' in given:
		binperx = 5
	if not 'plotcomps' in given:
		plotcomps = True
	if not 'plotabs' in given:
		plotabs = True
	
	
	xmajtickloc = MultLoc(xmajtick)
	xmintickloc = MultLoc(ymajtick)
	
	
	t1 = fits.open(infile)[1].data
	if corrname == 'none':
		pmra = t1['pmra']
		pmde = t1['pmde']
	else:
		pmra = t1['pmra_' + corrname]
		pmde = t1['pmde_' + corrname]
	pm = np.sqrt(pmra**2 + pmde**2)
	#pm = t1['jmag']
	
	del t1

	fig, ax = plt.subplots()
	binvals = np.linspace(xlim[0], xlim[1], int(xlim[1] - xlim[0])*binperx + 1)

	
	if plotcomps:
		if plotabs:
			ax.hist(pm, bins=binvals, color='blue', alpha=.5, 
					label=r'$\mathregular{|\mu|}$')
		ax.hist(pmra, bins=binvals, color='red', alpha=.5, 
				label=r'$\mathregular{\mu_{RA}}$')
		ax.hist(pmde, bins=binvals, color='green', alpha=.5, 
				label=r'$\mathregular{\mu_{DEC}}$')
		plt.legend()
	else:
		ax.hist(pm, bins=binvals, color='blue', 
				label=r'$\mathregular{|\mu|}$')

	ax.set_xlabel(xtit)
	ax.set_ylabel(ytit)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	ax.xaxis.set_major_locator(xmajtickloc)
	ax.xaxis.set_minor_locator(xmintickloc)
	
	if gridline:
		ax.axvline(30, linestyle='--', dashes=(3,2), linewidth=1, color='black')
		ax.axvline(0, linestyle='--', dashes=(3,2), linewidth=1, color='black')

	
	if save == False:
		plt.show()
	else:
		plt.savefig(save,bbox_inches='tight')





else:
	
	if not 'infile' in given:
		infile = 'dist_vs_match.fits'
	if not 'xtit' in given:
		xtit = 'Angular Separation (degrees)'
	if not 'ytit' in given:
		ytit = r'Difference in Proper Motion (mas yr$^{-1}$)'
	if not 'xmajtick' in given:
		xmajtick = .5
	if not 'xmintick' in given:
		xmintick = .1
	if not 'numbins' in given:
		numbins = 100
	if not 'xlim' in given:
		xlim = [0,10]
	
	minval, maxval = xlim
	
	xmajtickloc = MultLoc(xmajtick)
	xmintickloc = MultLoc(xmintick)

	t1 = fits.open(infile)[1].data
	dist_vs = t1['dist']
	match_vs = t1['match']
	del t1
	
	halfint = (maxval - minval)/2./numbins
	dist_center = np.linspace(minval + halfint, maxval - halfint, numbins)
	bins = np.linspace(minval, maxval, numbins+1)
	bin_place = np.digitize(dist_vs, bins) - 1
	
	avg_match = [0.] * numbins
	err_list = [0.] * numbins
	for i in xrange(numbins):
		sample = match_vs[ bin_place == i ]
		n = sample.size
		if n < 2:
			continue
			
		avg = sample.mean()
		std = sample.std()
		err = std / np.sqrt(n-1)
	
		avg_match[i] = avg
		err_list[i] = err

	fig, ax = plt.subplots()
	ax.scatter(dist_center, avg_match, c='blue', s=40)
	ax.errorbar(dist_center, avg_match, yerr=err_list, ecolor='blue', linestyle='none')
	
	ax.set_xlim((minval,maxval))
	ax.set_xlabel(xtit)
	ax.set_ylabel(ytit)
	ax.xaxis.set_major_locator(xmajtickloc)
	ax.xaxis.set_minor_locator(xmintickloc)
	if gridline:
		ax.axvline(2.5, linestyle='--', dashes=(3,2), linewidth=1, color='black')
	
	if save == False:
		plt.show()
	else:
		plt.savefig(save,bbox_inches='tight')
