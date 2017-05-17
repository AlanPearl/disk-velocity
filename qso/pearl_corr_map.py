#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from alan import read_argv
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mplcolors
import matplotlib.colorbar as mplcbar
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt

def formatAxis(ax, xlim, ylim, xlabel, ylabel, xmajtick, ymajtick, xmintick, ymintick):
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.xaxis.set_major_locator(MultLoc(xmajtick))
	ax.yaxis.set_major_locator(MultLoc(ymajtick))
	ax.xaxis.set_minor_locator(MultLoc(xmintick))
	ax.yaxis.set_minor_locator(MultLoc(ymintick))
	ax.tick_params(axis='both', which='minor', length=4)
	ax.tick_params(axis='both', which='major', length=8)

def makecmap(n):
	irange = xrange(256)
	#irange = xrange(255,-1,-1)
	
	cmap = []
	red0 = 1.0
	green0 = 0.4
	blue0 = 0.0
	redf = 0.0#1.0
	greenf = 0.7#0.0
	bluef = 0.7#1.0
	if n==1:
		for i in irange:
			if i < 128:
				red = red0 - red0 * (i/128.)**2
				green = green0 - green0 * (i/128.)**2
				blue = blue0 - blue0 * (i/128.)**2
			else:
				red = redf - redf * ((256-i)/128.)**2
				green = greenf - greenf * ((256-i)/128.)**2
				blue = bluef - bluef * ((256-i)/128.)**2
			cmap.append([red,green,blue])
	if n==2:
		for i in irange:
			red = redf - redf * ((256.-i)/256.)
			green = greenf - greenf * ((256.-i)/256.)
			blue = bluef - bluef * ((256.-i)/256.)
			#print red,green,blue
			cmap.append([red,green,blue])
	
	cmap = mplcolors.ListedColormap(cmap)
	if n==1:
		norm = mplcolors.Normalize(vmin=-10, vmax=10, clip=True)
	if n==2:
		cmap = plt.get_cmap('gnuplot_r')
		norm = mplcolors.Normalize(vmin=0, vmax=5, clip=True)
	return cmap,norm

def makecbar(cax, cmap, norm, n):
	cbar = mplcbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
	cbar.set_label(r'mas yr$^{-1}$', weight='medium')
	cbar.set_ticks(MultLoc(5./n/n))
	if n==1:
		cax.yaxis.get_ticklabels()[0].set_visible(False)
	else:
		cax.yaxis.get_ticklabels()[-2].set_label('100')

def makeplot(ax, ra, dec, val, cmap, norm, good):
	bina = np.linspace(0,360,1441)
	bind = np.linspace(-15,75,361)
	bins = [bina,bind]
	ax.hist2d(ra,dec,bins=bins,weights=val,cmap=cmap,norm=norm)#, cmin=1)
	#hist, xbins, ybins = np.histogram2d(ra,dec,bins=bins,weights=val)
	#extent = [xbins.min(),xbins.max(),ybins.min(),ybins.max()]
	#im = ax.imshow(hist[good.reshape((1440,360))].T, interpolation='none', origin='lower', extent=extent, cmap=cmap, norm=norm)


def createFigure(ra, dec, vals, textcolor, save=False):
	good = vals[-1] != 100.
	fig = plt.figure()
	cmap, norm = makecmap(1)
	cmap2, norm2 = makecmap(2)
	n = len(vals)
	gs = gridspec.GridSpec(2,3)
	gs.set_width_ratios([30,30,1])
	cax = plt.subplot(gs[2])
	makecbar(cax, cmap, norm, 1)
	cax2 = plt.subplot(gs[5])
	makecbar(cax2, cmap2, norm2, 2)
	axes = []
	xlim = [360,0]; ylim = [-15,75]
	xmajtick = 60; ymajtick = 30
	xmintick = ymintick = 5
	text = [r'pmra_shift', 
			r'pmde_shift', 
			r'pmra_se', 
			r'pmde_se']
	for i in xrange(n):
		gs_i = i
		if gs_i>1: gs_i += 1
		ax = plt.subplot(gs[gs_i])
		if i < 2:
			makeplot(ax, ra, dec, vals[i], cmap, norm, good)
		else:
			makeplot(ax, ra, dec, vals[i], cmap2, norm2, good)
		if i < 2:
			ax.tick_params(labelbottom='off', color=textcolor, which='both')
			ax.text(350,65,text[i],color=textcolor)
			xlabel = ''
		else:
			ax.tick_params(color=textcolor, which='both')
			ax.text(350,65,text[i],color=textcolor)
			xlabel = r'$\alpha$ (deg)'
		if i%2 == 1:
			ax.tick_params(labelleft='off')
			ylabel = ''
		else:
			ylabel = r'$\delta$ (deg)'
		if i==2:
			ax.xaxis.get_ticklabels()[1].set_visible(False)
		formatAxis(ax, xlim, ylim, xlabel, ylabel, xmajtick, ymajtick, xmintick, ymintick)
		ax.set_xlim(xlim)
	gs.update(wspace=0,hspace=0)
	
	if save==False:
		plt.show()
	else:
		plt.savefig(save,bbox_inches='tight')


if __name__ == '__main__':
	infile = '../data/pearl_corr.fits'
	outfile = False
	vars().update(read_argv(sys.argv))
	t1 = pyfits.open(infile)[1].data
	vals = []
	for i in xrange(4):
		vals.append(t1[t1.names[i+2]])
	
	ra = t1['ra']
	dec = t1['dec']
	createFigure(ra, dec, vals, 'white', save=outfile)