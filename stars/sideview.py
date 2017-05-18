#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from alan import read_argv, prse
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as mplcolors
import matplotlib.colorbar as mplcbar
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt
from matplotlib.gridspec import GridSpec

def get_aspect(xvar,yvar):
	RGC_AVG = 9.
	if xvar == 'T':
		r = RGC_AVG # The average R value
		return 180./np.pi / r
	else:
		return 1

def vvectorPlot(ax,cax):
	RGC_AVG = 9.
	corrname = 'pearl'
	selim = '8,10,10'
	set_pixel = False
	plot_errors = 'False'
	fontsize = 12
	marksize = 15
	errlim = 10000.
	vel_scale = 100. # The scale to divide velocities by before plotting
	vtheta_avg = -210
	Rtit = r'$\mathregular{R}$ (kpc)'
	Ttit = r'$\mathregular{\theta}$ (deg) - minor tick marks every 0.1 kpc'
	Ztit = r'$\mathregular{Z}$ (kpc)'
	veltit = r'$\mathregular{\langle V \rangle}$ (km s$^{-1}$)'
	errtit = r'$\mathregular{\langle V \rangle uncertainty}$ (km s$^{-1}$)'
	Rlim = '8,10'
	Tlim = '11.,-15'
	Zlim = '-2,2'
	Rstep = .2
	Tstep = -1.3
	Zstep = .2
	colors = ['red', 'black', 'blue']
	Rmajtick = 1.
	Rmintick = 0.1
	Zmajtick = 1.
	Zmintick = 0.1
	Tmajtick = 5.
	Tmintick = 1.
	Rfmt = '%d'
	Zfmt = '%d'
	Tfmt = '%d'

	xvar = 'T'
	yvar = 'Z'
	vardic = {'R':'r', 'Z':'z', 'T':'theta'}


	keys = ['corrname', 'selim', 'infile', 'marksize', 'pixel', 'fontsize', 'xlim', 'ylim', 'xmajtick', 'xmintick', 'ymajtick', 'ymintick', 'xfmt', 'yfmt', 'xvar', 'yvar', 'zvar', 'errlim', 'plot_errors']
	flagkeys = ['save']
	keycalled = {}
	for key in keys: keycalled[key] = False
	for key in flagkeys: keycalled[key] = False

	save = False
	hcbar = False
	reversecolor = False
	allstar = False
	tightz = False
	flags = []
	
	argdict = read_argv(sys.argv, keys=keys)
	for key in argdict:
		keycalled[key] = True
		if type(argdict[key]) is str:
			exec("%s = '%s'" %(key,argdict[key]))
		else:
			exec("%s = %s" %(key, str(argdict[key])))
	vtheta_avg = float(vtheta_avg)


	errlim = float(errlim)
	xvar = xvar[0].upper()
	yvar = yvar[0].upper()

	vardic1 = vardic.copy()
	if not keycalled['zvar']:
		try:
			del vardic1[xvar]
			del vardic1[yvar]
			zvar = vardic1.keys()[0]
		except:
			raise myError('There is something wrong with xvar,yvar = %s,%s' %(xvar,yvar))

	vardic_r = {}
	vardic_r[xvar] = 'x'
	vardic_r[yvar] = 'y'
	vardic_r[zvar] = 'z'
	
	if tightz:
		tightz = '_tightz'
	else:
		tightz = ''
	
	if not keycalled['infile']:
		if allstar:
			if xvar=='Z' and yvar=='R':
				infile = '%s-%s%s_%s%s.csv' %(vardic[yvar], vardic[xvar], tightz, corrname, selim)
			else:
				infile = '%s-%s%s_%s%s.csv' %(vardic[xvar], vardic[yvar], tightz, corrname, selim)
		else:
			if xvar=='Z' and yvar=='R':
				infile = '%s-%s%s_%sF%s.csv' %(vardic[yvar], vardic[xvar], tightz, corrname, selim)
			else:
				infile = '%s-%s%s_%sF%s.csv' %(vardic[xvar], vardic[yvar], tightz, corrname, selim)

	if not keycalled['save']:
		if corrname == 'pearl':
			zcorrname = 'zpearl'
		else:
			zcorrname = corrname
		
		if errlim == 10000.:
			zerrlim = 'ALL'
		else:
			zerrlim = str(errlim)
		outfile = 'plots3/%s-%s_%sF%s.png' %(vardic[xvar], vardic[yvar], corrname, selim)

	if not keycalled['xmajtick']:
		xmajtick = vars()[xvar + 'majtick']
	if not keycalled['ymajtick']:
		ymajtick = vars()[yvar + 'majtick']
	if not keycalled['xmintick']:
		xmintick = vars()[xvar + 'mintick']
	if not keycalled['ymintick']:
		ymintick = vars()[yvar + 'mintick']
	if not keycalled['xfmt']:
		xfmt = vars()[xvar + 'fmt']
	if not keycalled['yfmt']:
		yfmt = vars()[yvar + 'fmt']
	if not keycalled['xlim']:
		xlim = vars()[xvar + 'lim']
	if not keycalled['ylim']:
		ylim = vars()[yvar + 'lim']
	xtit = vars()[xvar + 'tit']
	ytit = vars()[yvar + 'tit']

	xmajtick = float(xmajtick)
	xmintick = float(xmintick)
	ymajtick = float(ymajtick)
	ymintick = float(ymintick)
	marksize = float(marksize)

	if keycalled['pixel']:
		pixel = pixel.split(',')
		set_pixel = True
		w = float(pixel[0])
		h = float(pixel[1])

	xlim = xlim.split(',')
	ylim = ylim.split(',')
	for i in [0,1]:
		xlim[i] = float(xlim[i])
		ylim[i] = float(ylim[i])

	# COLORBAR
	###########################################
	if reversecolor:
		irange = xrange(255,-1,-1)
	else:
		irange = xrange(256)

	cmap = []
	red0 = 1.0
	green0 = 0.4
	blue0 = 0.0
	redf = 0.0#1.0
	greenf = 0.7#0.0
	bluef = 0.7#1.0
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

	cmap = mplcolors.ListedColormap(cmap)
	###########################################
	norm = mplcolors.Normalize(vmin=-20, vmax=20, clip=True)
	norm2 = mplcolors.Normalize(vmin=0, vmax=10, clip=True)
	cmap2 = mplcm.gray_r
	ellipsecolor = 'grey'#[32./255.,178./255.,170./255.]
	errcolor = 'black'#'blue'
	#ellipsecolor = 'lime'
	#errcolor = 'green'
	
	coords = prse(infile, delim=',', dictin=True)


	X = np.array(coords[vardic[xvar]])
	Y = np.array(coords[vardic[yvar]])
	xvel1 = np.array(coords['v%s' %vardic[xvar]]) / vel_scale
	xerr1 = np.array(coords['v%s_err' %vardic[xvar]]) / vel_scale
	xse1 = np.array(coords['v%s_se' %vardic[xvar]]) / vel_scale
	yvel = np.array(coords['v%s' %vardic[yvar]]) / vel_scale
	yerr = np.array(coords['v%s_err' %vardic[yvar]]) / vel_scale
	yse = np.array(coords['v%s_se' %vardic[yvar]]) / vel_scale
	zvel = np.array(coords['v%s' %vardic[zvar]])
	zerr = np.array(coords['v%s_err' %vardic[zvar]])
	zse = np.array(coords['v%s_se' %vardic[zvar]])
	xunc1 = np.sqrt( xerr1**2 + xse1**2 )
	yunc = np.sqrt( yerr**2 + yse**2 )
	zunc = np.sqrt( zerr**2 + zse**2 )
	#xunc1 = xerr1 + xse1
	#yunc = yerr + yse
	#zunc = zerr + zse

	Tname = vardic_r['T'] + 'vel'
	if Tname == 'zvel':
		zvel -= vtheta_avg
	elif Tname == 'yvel':
		yvel -= vtheta_avg / vel_scale
	else:
		xvel1 -= vtheta_avg / vel_scale


	if xvar == 'T':
		r = RGC_AVG # The average R value
		aspect = 180./np.pi / r
		xvel = (xvel1) * aspect
		#xerr = xerr1 * aspect
		#xse  = xse1  * aspect
		xunc = xunc1 * aspect
		if not keycalled['xmintick']:
			xmintick = ymintick * aspect
	else:
		aspect = 1
		xvel = xvel1
		#xerr = xerr1
		#xse  = xse1
		xunc = xunc1

	#fig = plt.figure()
	#ax = fig.add_subplot(111, aspect=aspect)

	if set_pixel:
		DPI = fig.get_dpi()
		fig.set_size_inches(w/DPI, h/DPI)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	ax.set_xlabel(xtit, fontsize=fontsize, weight='medium')
	ax.set_ylabel(ytit, fontsize=fontsize, weight='medium')
	#ax.set_axis_bgcolor('black')

	#ax.set_axis_bgcolor('black')
	#ax.scatter (X[mask],Y[mask], c=err[mask], edgecolor='none', s=marksize*1.5, cmap=cmap2, norm=norm2)
	#ax.scatter (X,Y, c=vel, cmap=cmap, s=marksize, norm=norm)
	#ax2.scatter(X,Y, c=err, s=marksize, norm=norm2)


	#ax.quiver(X,Y,xerr,yerr, color='grey', edgecolor='red',
	#	units='xy', scale=100)
	#ax.quiver(X,Y,xse,yse, color='black', edgecolor='blue',
	#	units='xy', scale=100)

	if plot_errors[:1].lower() == 't': # Makes Error Ellipses
		ells = [Ellipse(xy=(X[i]+xvel[i],Y[i]+yvel[i]), zorder=1, linewidth=0, 
			width=xunc[i]*2, height=yunc[i]*2) for i in xrange(X.size)]
		for e in ells:
			ax.add_artist(e)
			e.set_clip_box(ax.bbox)
			e.set_alpha(1)
			e.set_facecolor(ellipsecolor)
		
		ax.errorbar(X+xvel, Y+yvel, yerr=zunc/vel_scale, 
			ecolor=errcolor, linestyle='none', linewidth=2, capsize=0)


	color = cmap(norm(zvel))
	edgecolor = color
	#edgecolor = mplcm.jet(norm2(zunc))
	#color[:,3] = 1. - norm2(zunc)
	ax.scatter(X,Y, c=color, edgecolor=edgecolor, s=marksize, zorder=2)
	if xvar == 'T':
		Q = ax.quiver(X,Y,xvel,yvel, color=color, zorder=3, 
			units='xy', scale=1, angles='xy', scale_units='xy', width=.1)
	else:
		Q = ax.quiver(X,Y,xvel,yvel, color=color, zorder=3, 
			units='xy', scale=1, angles='xy', scale_units='xy', width=.02)
	
	if yvar == 'R' and xvar == 'Z':
		x,y = -1.3*aspect,ylim[1]-0.4
		ax.scatter(x,y, c='black', s=marksize, zorder=4)
		plt.quiverkey(Q,x+.1*aspect,y, .2*aspect, r'20 km s$^{-1}$', coordinates='data', 
			color='black', zorder=5)
	elif yvar == 'R' and xvar == 'T':
		x,y = -1.6*aspect,ylim[1]-0.4
		ax.scatter(x,y, c='black', s=marksize, zorder=4)
		plt.quiverkey(Q,x+.1*aspect,y, .2*aspect, r'20 km s$^{-1}$', coordinates='data', 
			color='black', zorder=5)
	else:
		assert(yvar=='Z' and xvar=='T')
		x,y = -1.6*aspect,1.3
		ax.scatter(x,y, c='black', s=marksize, zorder=4)
		plt.quiverkey(Q,x+.1*aspect,y, .2*aspect, r'20 km s$^{-1}$', coordinates='data', 
			color='black', zorder=5)



	xmajtickloc = MultLoc(xmajtick)
	xmintickloc = MultLoc(xmintick)
	ymajtickloc = MultLoc(ymajtick)
	ymintickloc = MultLoc(ymintick)
	xmajtickfmt = StrFmt(xfmt)
	ymajtickfmt = StrFmt(yfmt)
	ax.xaxis.set_major_locator(xmajtickloc)
	ax.xaxis.set_minor_locator(xmintickloc)
	ax.xaxis.set_major_formatter(xmajtickfmt)
	ax.yaxis.set_major_locator(ymajtickloc)
	ax.yaxis.set_minor_locator(ymintickloc)
	ax.yaxis.set_major_formatter(ymajtickfmt)
	ax.tick_params(axis='both', which='minor', length=4)
	ax.tick_params(axis='both', which='major', length=8)

	ax.axhline(0, linestyle='--', dashes=(3,2), linewidth=1, color='black')
	ax.axhline(8, linestyle='--', dashes=(3,2), linewidth=1, color='black')
	ax.axvline(0, linestyle='--', dashes=(3,2), linewidth=1, color='black')

	texdic = {'R':r'R', 'Z':r'Z', 'T':r'\theta'}
	oops = ''
	if zvar == 'T':
		oops = r'} + %d \mathregular{' %(-vtheta_avg)
		oops = r' + %d ' %(-vtheta_avg)
	veltit = r'$\mathregular{\langle V_' + texdic[zvar] + r'\rangle' + oops + r'}$ (km s$^{-1}$)'
	if hcbar:
		cbar = mplcbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
	else:
		cbar = mplcbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
	cbar.set_label(veltit, fontsize=fontsize, weight='medium')
	cbar.set_ticks(MultLoc(5.))

if __name__ == '__main__':
	save = False
	r79 = False
	v200 = False
	outfile = "sideviewF.png"
	rz_file = "../data/r-z-metabins.csv"
	rtheta_file = "../data/r-theta-metabins.csv"
	thetaz_file = "../data/theta-z-metabins.csv"
	rlim1 = rlim2 = "8,10"
	thetalim2 = thetalim3 = "9,-12"
	zlim1 = "-1.6,1.6"; zlim3 = "-1.3,1.7"
	vtheta_avg = "-210"
	

	flagvars = {'save':'outfile'}
	vars().update(read_argv(sys.argv, flagvars=flagvars))

	fig = plt.figure(figsize=(6.7,12))
	gs = GridSpec(3,2)
	cax1 = plt.subplot(gs[1]); cax1.set_aspect(30)
	ax1 = plt.subplot(gs[0]); ax1.set_aspect(get_aspect('Z','R'))
	cax2 = plt.subplot(gs[3]); cax2.set_aspect(30)
	ax2 = plt.subplot(gs[2]); ax2.set_aspect(get_aspect('T','R'))
	cax3 = plt.subplot(gs[5]); cax3.set_aspect(45)
	ax3 = plt.subplot(gs[4]); ax3.set_aspect(get_aspect('T','Z'))
	
	
	sys.argv = "./sideview.py xvar=z yvar=r plot_errors=True".split()
	sys.argv.append("infile="+rz_file)
	sys.argv.append("xlim="+zlim1)
	sys.argv.append("ylim="+rlim1)
	sys.argv.append("vtheta_avg="+vtheta_avg)
	vvectorPlot(ax1,cax1)

	sys.argv = "./sideview.py xvar=theta yvar=r plot_errors=True".split()
	sys.argv.append("infile="+rtheta_file)
	sys.argv.append("xlim="+thetalim2)
	sys.argv.append("ylim="+rlim2)
	sys.argv.append("vtheta_avg="+vtheta_avg)
	vvectorPlot(ax2,cax2)

	sys.argv = "./sideview.py -reversecolor xvar=theta yvar=z plot_errors=True".split()
	sys.argv.append("infile="+thetaz_file)
	sys.argv.append("xlim="+thetalim3)
	sys.argv.append("ylim="+zlim3)
	sys.argv.append("vtheta_avg="+vtheta_avg)
	vvectorPlot(ax3,cax3)

	gs.set_height_ratios([1.03,1,1.5])
	gs.set_width_ratios([35,1])
	gs.update(wspace=-0.3, hspace=0.2)
	#plt.tight_layout()
	#plt.subplots_adjust(hspace=.3)

	if save:
		plt.savefig(outfile,bbox_inches='tight')
	else:
		plt.show()