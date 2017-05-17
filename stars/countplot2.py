#! /usr/bin/python2
"""

### countplot.py
### Alan Pearl

Plots the binned ra vs. dec of quasars with the color representing their
improvement after the Vickers correction.

Flags: (must come before keywords)
======
	-save[=FILE]-> save plot file

Keywords: (in default order)
=========
	rootname=	-> default: 'x-y-count', options: {x/theta}-{y/z}-{count}
	infile=		-> default: 'coords_pearl2.4.fits'
	marksize=	-> default: 100
	pixel=		-> default: 470,470
	fontsize=	-> defailt: 12
	xlim=		-> default: -9.8,-7.8
	ylim=		-> default: -2,2
	xbinnum=	-> default: 40
	ybinnum=	-> default: 80

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
from alan import prse, read_argv, rd2lb, bindata
import pyfits as fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar as mplcolorbar
import matplotlib.colors as mplcolors
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt

# Defining variables
infile = 'coords_pearl2.4.fits'
rootname = 'x-y-count' # Snapshot file root name
pixel = '470,470'
gridsize = 360
fontsize = 12
marksize = 100
xbinnum = 40
ybinnum = 80
xbinlen = .5
ybinlen = .5
xmintick = 0.1
ymintick = 0.1
xmajtick = 0.5
ymajtick = 0.5
cmax = 250
xfmt = '%.1f'
yfmt = '%.1f'

xtit = r'$\mathregular{X_{GC}}$ (kpc)'
ytit = r'$\mathregular{Y_{GC}}$ (kpc)'
xlim = '-10,-8'
ylim = '-2.,2.'


keys = ['rootname', 'infile', 'binfile', 'marksize', 'pixel', 'fontsize', 'xlim', 'ylim', 'xbinnum', 'ybinnum', 'xmintick', 'ymintick', 'xmajtick', 'ymajtick', 'xfmt', 'yfmt', 'xtit', 'ytit', 'cmax']

save = False
ASPECT = False
flagvars = {'save':'outfile'}
args = sys.argv
argdict = read_argv(args, keys=keys, flagvars=flagvars)
vars().update(argdict)


if not 'outfile' in argdict:
	outfile = '%s.png' %(rootname)

cmax = float(cmax)
xmintick = float(xmintick)
ymintick = float(ymintick)
xmajtick = float(xmajtick)
ymajtick = float(ymajtick)
xbinnum = int(float(xbinnum))
ybinnum = int(float(ybinnum))

if not 'binfile' in argdict:
	binfile = '%s.asc' %(rootname)

marksize = float(marksize)

pixel = pixel.split(',')
w = float(pixel[0])
h = float(pixel[1])

xlim = xlim.split(',')
ylim = ylim.split(',')
for i in [0,1]:
	xlim[i] = float(xlim[i])
	ylim[i] = float(ylim[i])


params = rootname.split('-')

if len(params) == 3:
	if params[0] == 'theta':
		if not 'xtit' in argdict:
			xtit = r'$\mathregular{\theta}$ (deg)'
		if not 'xmintick' in argdict:
			xmintick = 1.
		if not 'xmajtick' in argdict:
			xmajtick = 5.
		if not 'xfmt' in argdict:
			xfmt = '%d'
		if not 'xlim' in argdict:
			xlim = [11.,-14.5]
		if not 'xbinnum' in argdict:
			xbinnum = 102
	
	if params[0] == 'y':
		if not 'xtit' in argdict:
			xtit = r'$\mathregular{Y_{GC}}$ (kpc)'
		if not 'xmajtick' in argdict:
			xmajtick = 1.
		if not 'xfmt' in argdict:
			xfmt = '%d'
		if not 'xlim' in argdict:
			xlim = [-2.,2.]
		if not 'xbinnum' in argdict:
			xbinnum = 80
	
	if params[0] == 'l':
		if not xtit in argdict:
			xtit = r'$l$ (deg)'
		if not ytit in argdict:
			ytit = r'$\mathregular{b}$ (deg)'
		if not 'xfmt' in argdict:
			xfmt = '%d'
		if not 'yfmt' in argdict:
			yfmt = '%d'
	
	if params[1] == 'z':
		if not 'ytit' in argdict:
			ytit = r'$\mathregular{Z_{GC}}$ (kpc)'
else:
	sys.exit('Rootname must be in form {x/y/theta}-{y/z}-{count} (3 Parameters)')


##############

table = fits.open(infile)[1].data
t1 = {}
for name in table.names:
	t1[name] = table[name]
del table
t1['l'], t1['b'] = rd2lb(t1['ra'],t1['dec'], deg=True)

x1, y1 = t1[params[0]], t1[params[1]]

del t1

###############

#x,y,z = prse(binfile)
(x,y),_,_,z = bindata([x1,y1], binlens=[1,1], lims=[[0,360],[-90,90]], mincount=0)

if params[0] == 'l':
	xarray = np.linspace(0,360,xbinnum+1)
	binlen = xarray[1] - xarray[0]
	xarray = np.linspace(0+.5*binlen, 360-.5*binlen, xbinnum)
	yarray = np.linspace(-90,90,ybinnum+1)
	binlen = xarray[1] - xarray[0]
	yarray = np.linspace(-90+.5*binlen, 90-.5*binlen, ybinnum)
else:
	xarray = np.linspace(xlim[0], xlim[1], xbinnum+1)
	binlen = xarray[1] - xarray[0]
	xarray = np.linspace(xlim[0]+.5*binlen, xlim[1]-.5*binlen, xbinnum)
	yarray = np.linspace(ylim[0], ylim[1], ybinnum)
	binlen = yarray[1] - yarray[0]
	yarray = np.linspace(ylim[0]+.5*binlen, ylim[1]-.5*binlen, ybinnum)

zarray = []
for i in xrange(ybinnum):
	zarray.append([])
	for j in xrange(xbinnum):
		zarray[i].append( z[i + j*ybinnum] )

####################

if params[0] == 'theta':
	aspect = 6.437728035177789
else:
	aspect = 1

fig = plt.figure()
if ASPECT:
	ax = fig.add_subplot(111, aspect=aspect)
else:
	ax = fig.add_subplot(111)

DPI = fig.get_dpi()
fig.set_size_inches(w/DPI, h/DPI)
#ax.set_xlim(xlim) # set at the end
#ax.set_ylim(ylim)
ax.set_xlabel(xtit, fontsize=fontsize, weight='medium')
ax.set_ylabel(ytit, fontsize=fontsize, weight='medium')




#v = [1,1000,2500,4000,5000,6000,7000,8000,9000,10000,999999999]
v = [1,25,50,100,250,500,1000,2500,5000,8000,999999999]
#v = [1,10,30,50,80,120,160,200,240,280,320,360,400,999999999]

ax.contour(xarray,yarray,zarray, v, colors='black')

norma = mplcolors.PowerNorm(gamma=.55, vmin=0, vmax=cmax, clip=True)

xbins = np.arange(min(xlim), max(xlim)+.5*abs(xbinlen), abs(xbinlen))
ybins = np.arange(min(ylim), max(ylim)+.5*abs(ybinlen), abs(ybinlen))
ax.hist2d(x1, y1, bins=[xbins,ybins], norm=norma, cmin=1)

xmintickloc = MultLoc(xmintick)
ymintickloc = MultLoc(ymintick)
xmajtickloc = MultLoc(xmajtick)
ymajtickloc = MultLoc(ymajtick)
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
ax.axvline(0, linestyle='--', dashes=(3,2), linewidth=1, color='black')
ax.axvline(180, linestyle='--', dashes=(3,2), linewidth=1, color='black')

area = float(xbinlen*ybinlen)
if params[0] == 'theta':
	area *= 6.437728035177789

normCbar = mplcolors.PowerNorm(gamma=.55, vmin=0, vmax=cmax/area, clip=True)
cax, kw = mplcolorbar.make_axes(ax, location='right', aspect=40, pad=0.01)
cbar = mplcolorbar.ColorbarBase(cax, norm=normCbar,	orientation='vertical')

if params[0] == 'l':
	cbar.set_label(r'N (deg$^{-2}$)', fontsize=fontsize, weight='medium')
else:
	cbar.set_label(r'N (kpc$^{-2}$)', fontsize=fontsize, weight='medium')

ax.set_xlim(xlim)
ax.set_ylim(ylim)


if save:
	plt.savefig(outfile,bbox_inches='tight')
else:
	plt.show()
