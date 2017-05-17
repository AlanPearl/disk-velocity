#! /usr/bin/python2
"""
### Alan Pearl

Functions to recreate Figure 2 (Carlin et al. 2013) with vector plots

Flags: (must come before keywords)
======
	-reversecolor-> reverse the color scale (also flips the z axis so we view from the right)
	-allstar	-> Use all star types instead of only F stars

Keywords: (in default order)
=========
	marksize=	-> default: 15
	pixel=		-> default: {Automatic} format:900,450
	fontsize=	-> defailt: 12
	xlim=		-> default: 11,-14.5    or  7.7,9.9
	ylim=		-> default: -2,2

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
from alan import read_argv, prse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as mplcolors
import matplotlib.colorbar as mplcbar
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt
from matplotlib.gridspec import GridSpec


def cond(z,zlim):
	if zlim[1] > zlim[0]:
		if z < zlim[1]-.01:
			return True
		else:
			return False
	else:
		if z > zlim[1]+.01:
			return True
		else:
			return False

def vvectorPlot2(fig,argdict):
	"""args: R
			 theta
			 z
			 vR
			 vtheta
			 vz
			 vR_err
			 vtheta_err
			 vz_err
			 vR_se
			 vtheta_se
			 vz_se
	"""
	# Defining variables
	set_pixel = False
	plot_errors = 'True'
	fontsize = 12
	marksize = 15
	vel_scale = 100. # The scale to divide velocities by before plotting
	vtheta_avg = -210
	RGC_avg = 9.
	Rtit = r'$\mathregular{R}$ (kpc)'
	Ttit = r'$\mathregular{\theta}$ (deg)'
	Ztit = r'$\mathregular{Z}$ (kpc)'
	veltit = r'$\mathregular{<V>}$ (km s$^{-1}$)'
	errtit = r'$\mathregular{<V> uncertainty}$ (km s$^{-1}$)'
	Rlim = '8,10'
	Tlim = '4.5,-5.9'
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
	zvar = 'R'
	vardic = {'R':'r', 'Z':'z', 'T':'theta'}


	keys = ['xvar','yvar','zvar','marksize', 'pixel', 'fontsize', 'xlim', 'ylim', 'zlim', 'xmajtick', 'xmintick', 'ymajtick', 'ymintick', 'xfmt', 'yfmt', 'plot_errors']
	flagvars = {}
	keycalled = {}
	for key in keys: keycalled[key] = False
	for key in flagvars: keycalled[key] = False

	reversecolor = False
	allstar = False
	flags = []
	args = sys.argv
	args = read_argv(args, keys=keys, flagvars=flagvars)
	for key in args: keycalled[key] = True
	
	#locals().update(args) <--- Not working >:(
	for key in args:
		exec("%s = '%s'" %(key,args[key]))
	vtheta_avg = float(vtheta_avg)

	xvar = xvar[0].upper()
	yvar = yvar[0].upper()
	zvar = zvar[0].upper()

	vardic1 = vardic.copy()
	if not keycalled['zvar']:
		try:
			del vardic1[xvar]
			del vardic1[yvar]
			zvar = vardic1.keys()[0]
		except:
			raise IOError('There is something wrong with xvar,yvar ==> %s,%s' %(xvar,yvar))

	vardic_r = {}
	vardic_r[xvar] = 'x'
	vardic_r[yvar] = 'y'
	vardic_r[zvar] = 'z'


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
	if not keycalled['zlim']:
		zlim = vars()[zvar + 'lim']
	xtit = vars()[xvar + 'tit']
	ytit = vars()[yvar + 'tit']
	ztit = vars()[zvar + 'tit']

	xmajtick = float(xmajtick)
	xmintick = float(xmintick)
	ymajtick = float(ymajtick)
	ymintick = float(ymintick)
	marksize = float(marksize)

	fontsize = float(fontsize)
	if keycalled['pixel']:
		pixel = pixel.split(',')
		set_pixel = True
		w = float(pixel[0])
		h = float(pixel[1])

	zstep = vars()[zvar + 'step']
	xlim = xlim.split(',')
	ylim = ylim.split(',')
	zlim = zlim.split(',')
	for i in [0,1]:
		xlim[i] = float(xlim[i])
		ylim[i] = float(ylim[i])
		zlim[i] = float(zlim[i])

	# COLORBAR
	###########################################
	if reversecolor:
		irange = xrange(255,-1,-1)
		xlim.reverse()
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
	
	X = np.array(argdict[vardic[xvar]])
	Y = np.array(argdict[vardic[yvar]])
	Z = np.array(argdict[vardic[zvar]])
	xvel1 = np.array(argdict['v%s' %vardic[xvar]]) / vel_scale
	xerr1 = np.array(argdict['v%s_err' %vardic[xvar]]) / vel_scale
	xse1 = np.array(argdict['v%s_se' %vardic[xvar]]) / vel_scale
	yvel = np.array(argdict['v%s' %vardic[yvar]]) / vel_scale
	yerr = np.array(argdict['v%s_err' %vardic[yvar]]) / vel_scale
	yse = np.array(argdict['v%s_se' %vardic[yvar]]) / vel_scale
	zvel = np.array(argdict['v%s' %vardic[zvar]])
	zerr = np.array(argdict['v%s_err' %vardic[zvar]])
	zse = np.array(argdict['v%s_se' %vardic[zvar]])
	xunc1 = np.sqrt( xerr1**2 + xse1**2 )
	yunc = np.sqrt( yerr**2 + yse**2 )
	zunc = np.sqrt( zerr**2 + zse**2 )

	Tname = vardic_r['T'] + 'vel'
	if Tname == 'zvel':
		zvel -= vtheta_avg
	elif Tname == 'yvel':
		yvel -= vtheta_avg / vel_scale
	else:
		xvel1 -= vtheta_avg / vel_scale


	z = zlim[0]
	i_main = 0
	subplotNum = 1
	axes = []
	assert(abs(abs(zlim[0]-zlim[1])-1.3*8)<.0001)
	while cond(z,zlim):
		i_main += 1
		
		if xvar == 'T':
			if zvar == 'R':
				r = RGC_avg # z + .5 * zstep
			else:
				r = RGC_avg # The average R value
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

		if zlim[1] > zlim[0]:
			mask = np.where(np.logical_and( z < Z , Z <= z+zstep ))
		else:
			mask = np.where(np.logical_and( z > Z , Z >= z+zstep ))

		ax = fig.add_subplot(4,2,subplotNum, aspect=aspect)
		axes.append(ax)
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

		#print 'Subplot number', subplotNum
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		if subplotNum > 6:
			#print '\tPlotting x labels'
			ax.set_xlabel(xtit, fontsize=fontsize, weight='medium')
			if subplotNum != 8:
				# For some reason, get_major_ticks returns ticks: [-3, -2, -1, 0, 1, 2, 3]
				# Where -3 and 3 seem to already have invisible labels
				ax.xaxis.get_major_ticks()[5].label.set_visible(False)
				#ax.xaxis.get_major_ticks()[3].label1.set_visible(False)
				#ax.xaxis.get_major_ticks()[3].label2.set_visible(False)
		else:
			ax.tick_params(labelbottom='off')
		if subplotNum % 2 == 0:
			ax.tick_params(labelleft='off') 
		else:
			#print '\tPlotting y labels'
			ax.set_ylabel(ytit, fontsize=fontsize, weight='medium')
			if subplotNum != 1:
				ax.yaxis.get_major_ticks()[-2].label.set_visible(False)
		
		if subplotNum > 6:
			subplotNum -= 5
		else:
			subplotNum += 2
		
		X2 = X[mask]
		Y2 = Y[mask]
		xvel2 = xvel[mask]
		yvel2 = yvel[mask]
		zvel2 = zvel[mask]
		xunc2 = xunc[mask]
		yunc2 = yunc[mask]
		zunc2 = zunc[mask]
		

		if plot_errors[:1].lower() == 't': # Makes Error Ellipses
			ells = [Ellipse(xy=(X2[i]+xvel2[i],Y2[i]+yvel2[i]), zorder=1, linewidth=0, 
				width=xunc2[i]*2, height=yunc2[i]*2) for i in xrange(X2.size)]
			for e in ells:
				ax.add_artist(e)
				e.set_clip_box(ax.bbox)
				e.set_alpha(1)
				e.set_facecolor(ellipsecolor)
		
			ax.errorbar(X2+xvel2, Y2+yvel2, yerr=zunc2/vel_scale, 
				ecolor=errcolor, linestyle='none', linewidth=2, capsize=0)


		color = cmap(norm(zvel2))
		edgecolor = color
		#edgecolor = mplcm.jet(norm2(zunc))
		#color[:,3] = 1. - norm2(zunc)
		ax.scatter(X2,Y2, c=color, edgecolor=edgecolor, s=marksize, zorder=2)
		if xvar == 'T':
			Q = ax.quiver(X2,Y2,xvel2,yvel2, color=color, zorder=3, 
				units='xy', scale=1, angles='xy', scale_units='xy', width=.1)
		else:
			Q = ax.quiver(X2,Y2,xvel2,yvel2, color=color, zorder=3, 
				units='xy', scale=1, angles='xy', scale_units='xy', width=.02)

		if xvar == 'R':
			assert(False)
			x,y = 9.6,-1.8
			ax.scatter(x,y, c='black', s=marksize, zorder=4)
			plt.quiverkey(Q,x+.1,y, .2, r' ', coordinates='data', 
				color='black', zorder=5)
		elif yvar == 'R':
			x,y = -1.9*aspect,ylim[0]+0.1
			ax.scatter(x,y, c='black', s=marksize, zorder=4)
			plt.quiverkey(Q,x+.1*aspect,y, .2*aspect, r' ', coordinates='data', 
				color='black', zorder=5)
		else:
			assert(False)
			x,y = -1.9*aspect,-1.8
			ax.scatter(x,y, c='black', s=marksize, zorder=4)
			plt.quiverkey(Q,x+.1*aspect,y, .2*aspect, r' ', coordinates='data', 
				color='black', zorder=5)
		
		if reversecolor:
			ax.text(x+0.65,y+.05, r'20 km s$^{-1}$', fontsize=fontsize*.8)
		else:
			ax.text(x-0.09,y+.05, r'20 km s$^{-1}$', fontsize=fontsize*.8)
		
		if zstep < 0:
			z1 = z + zstep
			z2 = z
		else:
			z1 = z
			z2 = z + zstep
		
		if reversecolor:
			xtext = xlim[0]-2.4#(xlim[1] + xlim[0]) / 2.
		else:
			xtext = xlim[0]+.25#(xlim[1] + xlim[0]) / 2.
		ytext = ylim[1] - (ylim[1] - ylim[0]) / 10.
		
		#ax.text ( xtext, ytext, '%.1f < %s < %.1f' %(z1, ztit, z2), weight='bold', fontsize=fontsize )
		ax.text ( xtext, ytext, r'$'+str(z2)+r'^{\circ} > \theta  > '+str(z1)+r'^{\circ}$', weight='bold', fontsize=fontsize )

		ax.axhline(0, linestyle='--', dashes=(3,2), linewidth=1, color='black')
		ax.axhline(8, linestyle='--', dashes=(3,2), linewidth=1, color='black')
		ax.axvline(0, linestyle='--', dashes=(3,2), linewidth=1, color='black')

		texdic = {'R':r'R', 'Z':r'Z', 'T':r'\theta'}
		oops = r''
		if zvar == 'T':
			oops = r'} + %d \mathregular{' %(-vtheta_avg)
			oops = r' + %d ' %(-vtheta_avg)
		veltit = r'$\mathregular{\langle V_' + texdic[zvar] + r'\rangle' + oops + r'}$ (km s$^{-1}$)'
		
		z += zstep

	plt.subplots_adjust(wspace=0, hspace=0)

	cax, kw = mplcbar.make_axes(tuple(axes), location='top', aspect=40, pad=.05)
	cbar = mplcbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
	cbar.set_label(veltit, fontsize=fontsize, weight='medium')
	cbar.set_ticks(MultLoc(5.))

	if set_pixel:
		DPI = fig.get_dpi()
		fig.set_size_inches(w/DPI, h/DPI)

def get_aspect(xvar,yvar):
	RGC_AVG = 9.
	if xvar == 'T':
		r = RGC_AVG # The average R value
		return 180./np.pi / r
	else:
		return 1

def vvectorPlot(ax,cax):
	# Defining variables
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