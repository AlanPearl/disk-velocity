#! /usr/bin/python2
"""

Display proper motion of quasars at 
each position in the sky.

Flags: (must come before keywords)
======
	-save[=FILE]-> save plot file

Keywords: (in default order)
=========
	corrname=	-> default: fix
	infile=
	marksize=	-> default: 5
	pixel=		-> default: 470,470
	fontsize=	-> defailt: 12
	xlim=		-> default: 80,270
	ylim=		-> default: -70,70
	ra_dec=		-> default: False

"""

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from alan import prse, lb2rd
from math import sqrt
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as mplcolors
import matplotlib.colorbar as mplcolorbar
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt


# Defining variables
delete_bad = False
EBPV = 2.4 # Eliminate Bad Pearl Value (fit_err)
ra_dec = 'False'
corrname = 'fix'
fixcode = 'fix5'
selectcode = 'lt1'
pixel = '1300,700'
ax1tit = r'$\mathregular{\mu_{\alpha} cos \delta}$'
ax2tit = r'$\mathregular{\mu_{\delta}}$'
ztit = r'$\mathregular{\mu}$ (mas yr$^{-1}$)'
ztita = r'$\mathregular{|\mu{}|}$ (mas yr$^{-1}$)'
fontsize = 12
marksize = 10
xlim = '80,270'
ylim = '-70,70'
xlim = '0,360'
ylim = '-90,90'
xmajtickloc = MultLoc(30.)
xmintickloc = MultLoc(1.)
ymajtickloc = MultLoc(30.)
ymintickloc = MultLoc(1.)
gridsize1 = '360'
gridsize2 = '90'

xmajtickfmt = StrFmt('%d')
ymajtickfmt = StrFmt('%d')

keys = ['corrname=', 'infile=', 'marksize=', 'pixel=', 'fontsize=', 'xlim=', 'ylim=', 'ra_dec=']
flagkeys = ['save=']
keycalled = {}
for key in keys: keycalled[key] = False
for key in flagkeys: keycalled[key] = False

save = False
flags = []
args = sys.argv

# Reading arguments
#===========================================================
if len(args)>0:
	args = args[1:]
else:
	args = []

if len(args) > 0:
	while len(args)>0 and args[0][0] == '-':
		flags.append(args[0][1:])
		args = args[1:]
	
	for flag in flags:
		flen = len(flag)
		if flen>3 and 'save' == flag[:4].lower():
			save = True
			if len(flag)>5 and flag[:5] == 'save=':
				outfile = flag[5:]
				keycalled['save='] = True
		else:
			print 'Ignored flag: -' + flag
			
	
	while len(args)>0 \
	and not '=' in args[0] \
	and len(keys)>0:
		exe = "%s'%s'" %(keys[0],args[0])
		exec(exe)
		keys = keys[1:]
		args = args[1:]
	
	for arg in args:
		nokey = True
		for key in keys:
			if key in arg:
				a = arg[len(key):]
				exe = "%s'%s'" %(key, a)
				exec(exe)
				nokey = False
				keycalled[key] = True
		if nokey:
			print 'Ignored argument:', arg
#===========================================================
gridsize1 = int(gridsize1)
gridsize2 = int(gridsize2)
if ra_dec[:1].lower() == 't':
	ra_dec = True
else:
	ra_dec = False

if not keycalled['infile=']:
	infile = 'qsocoords_' + selectcode + '_' + fixcode + '_dedup.fits'

if not keycalled['save=']:
	ofp = (infile[:3], infile[9:-11])
	if corrname == '':
		corrname1 = 'none'
	else:
		corrname1 = corrname
	
	outfile = ofp[0] + ofp[1] + '_' + corrname1


marksize = float(marksize)

pixel = pixel.split(',')
w = float(pixel[0])
h = float(pixel[1])

xlim = xlim.split(',')
ylim = ylim.split(',')
for i in [0,1]:
	xlim[i] = float(xlim[i])
	ylim[i] = float(ylim[i])

if ra_dec:
	xtit = r'$\mathregular{\alpha}$ (deg)'
	ytit = r'$\mathregular{\delta}$ (deg)'
else:
	xtit = r'$\mathregular{l}$ (deg)'
	ytit = r'$\mathregular{b}$ (deg)'



#Color Mapping & Normalization
############################################################
cmap1 = []
red0 = 1.0
green0 = 0.4
blue0 = 0.0
redf = 0.0
greenf = 0.7
bluef = 0.7
for i in xrange(256):
	if i < 128:
		red = red0 - red0 * (i/128.)**2
		green = green0 - green0 * (i/128.)**2
		blue = blue0 - blue0 * (i/128.)**2
	else:
		red = redf - redf * ((256-i)/128.)**2
		green = greenf - greenf * ((256-i)/128.)**2
		blue = bluef - bluef * ((256-i)/128.)**2
	
	cmap1.append([red,green,blue])

cmap1 = mplcolors.ListedColormap(cmap1)
###########################################

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
cmapa = mplcolors.ListedColormap(cvals)



norm1 = mplcolors.Normalize(vmin=-6, vmax=6, clip=True)
norma = mplcolors.Normalize(vmin=0, vmax=10, clip=True)
############################################################


#Reading input files
#===========================================================
t1 = pyfits.open(infile)[1].data

if delete_bad:
	mask1 = t1['fit_err'] < EBPV
else:
	mask1 = (np.array(range(len(t1))),)

if ra_dec:
	l1 = t1['ra'][mask1]
	b1 = t1['dec'][mask1]
else:
	l1 = t1['l'][mask1]
	b1 = t1['b'][mask1]
#pmra1 = [0.] * len(mask1[0])
pmra1 = t1['pmra'][mask1]
#pmde1 = [0.] * len(mask1[0])
pmde1 = t1['pmde'][mask1]
if corrname == '' or corrname.lower() == 'none':
	pmra_c1 = t1['pmra'][mask1]
	pmde_c1 = t1['pmde'][mask1]
else:
	pmra_c1 = t1['pmra_'+corrname][mask1]
	pmde_c1 = t1['pmde_'+corrname][mask1]

pmra_change = [] ; pmde_change = []
for i in xrange(len(pmra_c1)):
	pmra_change_i = pmra_c1[i] - pmra1[i]
	pmde_change_i = pmde_c1[i] - pmde1[i]
	
	pmra_change.append( pmra_change_i )
	pmde_change.append( pmde_change_i )



del t1
#===========================================================

#===========================================================

#lt1count = 0 ; lt2count = 0 ; lt3count = 0 ; lt4count = 0 ; lt5count = 0 ; gt5count = 0
pma = [None] * len(l1)
for i in xrange(len(l1)):
	pmra_c1_i = pmra_c1[i]
	pmde_c1_i = pmde_c1[i]
	pma[i] = sqrt( pmra_c1_i*pmra_c1_i + pmde_c1_i*pmde_c1_i )
	# if pma[i] < 2.: lt1count += 1
	# elif pma[i] < 4.: lt2count += 1
	# elif pma[i] < 6.: lt3count += 1
	# elif pma[i] < 8.: lt4count += 1
	# elif pma[i] < 10.: lt5count += 1
	# else: gt5count += 1
	#if pma[i] < 50.: lt1count += 1
	#elif pma[i] < 70.: lt2count += 1
	#elif pma[i] < 90.: lt3count += 1
	#elif pma[i] < 110.: lt4count += 1
	#elif pma[i] < 130.: lt5count += 1
	#else: gt5count += 1

avgpma = sum(pma)/len(l1)
print 'Average proper motion:', avgpma

#===========================================================

############################################################
###############  PLOTTING DATA  ############################
############################################################

fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True, figsize=(10,4.1))#, subplot_kw=dict(projection='aitoff'))
#DPI = fig.get_dpi()
#fig.set_size_inches(w/DPI, h/DPI)
#figa, axa = plt.subplots()
#DPIa = figa.get_dpi()
#figa.set_size_inches(w/DPIa, h/DPIa)
ax1.set_aspect(1)
ax2.set_aspect(1)

ax1.hexbin(l1, b1, pmra_c1, norm=norm1, cmap=cmap1, gridsize=(gridsize1,gridsize2), linewidth=0.1)
ax2.hexbin(l1, b1, pmde_c1, norm=norm1, cmap=cmap1, gridsize=(gridsize1,gridsize2), linewidth=0.1)

############################################################
############################################################
############################################################

#Formatting axes
#===========================================================
ax1.set_xlim(xlim)
ax1.set_ylim(ylim)
ax1.set_ylabel(ytit, fontsize=fontsize, weight='medium')
ax2.set_xlim(xlim)
ax2.set_ylim(ylim)
ax2.set_xlabel(xtit, fontsize=fontsize, weight='medium')
ax2.set_ylabel(ytit, fontsize=fontsize, weight='medium')
ax1.text(180,63,ax1tit, fontsize=fontsize,ha='center')
ax2.text(180,63,ax2tit, fontsize=fontsize,ha='center')

ax1.xaxis.set_major_locator(xmajtickloc)
ax1.xaxis.set_minor_locator(xmintickloc)
ax1.xaxis.set_major_formatter(xmajtickfmt)
ax1.yaxis.set_major_locator(ymajtickloc)
ax1.yaxis.set_minor_locator(ymintickloc)
ax1.yaxis.set_major_formatter(ymajtickfmt)
#ax1.yaxis.get_major_ticks()[1].label.set_visible(False)
ax1.tick_params(axis='both', which='minor', length=4)
ax1.tick_params(axis='both', which='major', length=8)
ax2.xaxis.set_major_locator(xmajtickloc)
ax2.xaxis.set_minor_locator(xmintickloc)
ax2.xaxis.set_major_formatter(xmajtickfmt)
ax2.yaxis.set_major_locator(ymajtickloc)
ax2.yaxis.set_minor_locator(ymintickloc)
ax2.yaxis.set_major_formatter(ymajtickfmt)
ax2.tick_params(axis='both', which='minor', length=4)
ax2.tick_params(axis='both', which='major', length=8)

if ra_dec:
	centerline = lb2rd( range(360), [0]*360, deg=True)
	upperline = lb2rd( range(360), [15]*360, deg=True)
	lowerline = lb2rd( range(360), [-15]*360, deg=True)
	decline1 = list(lb2rd( [180]*181, range(-90,91), deg=True))
	decline2 = list(lb2rd( [0]*181, range(-90,91), deg=True))
	decline = np.concatenate([decline1, decline2], axis=1)
	
	centerline0 = [x for (x,y) in sorted(zip(centerline[0],centerline[1]))]
	centerline1 = [y for (x,y) in sorted(zip(centerline[0],centerline[1]))]
	upperline0 = [x for (x,y) in sorted(zip(upperline[0],upperline[1]))]
	upperline1 = [y for (x,y) in sorted(zip(upperline[0],upperline[1]))]
	lowerline0 = [x for (x,y) in sorted(zip(lowerline[0],lowerline[1]))]
	lowerline1 = [y for (x,y) in sorted(zip(lowerline[0],lowerline[1]))]
	decline0 = [x for (x,y) in sorted(zip(decline[0],decline[1]))]
	decline1 = [y for (x,y) in sorted(zip(decline[0],decline[1]))]
	
	ax1.plot(centerline0, centerline1, linestyle='--', linewidth=1, color='grey')
	ax1.plot(upperline0, upperline1, linestyle='--', linewidth=1, color='grey')
	ax1.plot(lowerline0, lowerline1, linestyle='--', linewidth=1, color='grey')
	ax1.plot(centerline0, centerline1, linestyle='--', linewidth=1, color='grey')
	ax1.plot(upperline0, upperline1, linestyle='--', linewidth=1, color='grey')
	ax1.plot(lowerline0, lowerline1, linestyle='--', linewidth=1, color='grey')
	ax1.plot(decline0, decline1, linestyle='--', linewidth=1, color='grey')
	ax2.plot(centerline0, centerline1, linestyle='--', linewidth=1, color='grey')
	ax2.plot(upperline0, upperline1, linestyle='--', linewidth=1, color='grey')
	ax2.plot(lowerline0, lowerline1, linestyle='--', linewidth=1, color='grey')
	ax2.plot(centerline0, centerline1, linestyle='--', linewidth=1, color='grey')
	ax2.plot(upperline0, upperline1, linestyle='--', linewidth=1, color='grey')
	ax2.plot(lowerline0, lowerline1, linestyle='--', linewidth=1, color='grey')
	ax2.plot(decline0, decline1, linestyle='--', linewidth=1, color='grey')

else:
	ax1.axhline(0, linestyle='--', dashes=(3,2), linewidth=1, color='grey')
	ax1.axvline(180, linestyle='--', dashes=(3,2), linewidth=1, color='grey')
	ax2.axhline(0, linestyle='--', dashes=(3,2), linewidth=1, color='grey')
	ax2.axvline(180, linestyle='--', dashes=(3,2), linewidth=1, color='grey')

plt.subplots_adjust(hspace=0)

cax, kw = mplcolorbar.make_axes((ax1,ax2), location='right', aspect=40, pad=0.01)
cbar = mplcolorbar.ColorbarBase(cax, cmap=cmap1, norm=norm1,
	orientation='vertical', label=ztit)
cbar.set_label(ztit, fontsize=fontsize, weight='medium')

#===========================================================

#plt.gca().set_rasterized(True)

if save:
	plt.savefig(outfile,bbox_inches='tight')
else:
	plt.show()



