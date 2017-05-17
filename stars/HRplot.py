#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar as mplcolorbar
import matplotlib.colors as mplcolors
import pyfits
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt
from matplotlib import gridspec
from alan import read_argv

ScatterPlot = False
JKcolor = True
infile = 'HR_paramsF.fits'
outfile = 'HRplot.png'
save = False
Rmin = -10
Rmax = 10
GRmin = -100
GRmax = 100
xlim = '-10,10'
ylim = '-5,5'
xbinlen = .03
ybinlen = .1
fontsize = 12
xmintick = .1
ymintick = 1
xtit = r'Color $\mathregular{(J0-K0)}$'
ytit = r'$\mathregular{M_K}$'
ztit = r'N'

keys = ['xlim','ylim','xbinlen','ybinlen','xmintick','ymintick','xtit','ytit']
flagvars = {'save':'outfile'} # Flags: -ScatterPlot, -JKcolor
args = sys.argv
argdict = read_argv(args, keys=keys, flagvars=flagvars)
vars().update(argdict)

xlim = map(float,xlim.split(','))
ylim = map(float,ylim.split(','))
xbinlen = float(xbinlen)
ybinlen = float(ybinlen)
xmintick = float(xmintick)
ymintick = float(ymintick)

t1 = pyfits.open(infile)[1].data
Kmag = t1['M_K50']
#Rmag = t1['rmag']
#g0 = t1['g0']	# ALL REAL: (max,min) = (99.97,-171.16); median = 13.70
#r0 = t1['r0']	# ALL REAL: (max,min) = (99.97,-123.51); median = 13.41
#GR = g0 - r0	# ALL REAL: (max,min) = (98.96,-47.65); median = 0.25
j0 = t1['j0']	# Take out 338 NaN: (max,min) = (16.54,-29.92); median = 12.37
k0 = t1['k0']	# Take out 338 NaN: (max,min) = (16.01,-3.34); median = 12.15
GR = j0 - k0	# Take out 338 NaN: (max,min) = (2.07,-26.58); median = 0.26
#b = t1['b']
del t1

good1 = GR==GR
good2 = good1# np.logical_and(
#	np.logical_and(Rmin <= Rmag, Rmag <= Rmax),
#	np.logical_and(GRmin <= GR, GR <= GRmax)
#)
good = np.array([good1,good2])
good = np.all(good,axis=0)

if not ScatterPlot:
	xbins = np.arange(GRmin,GRmax+xbinlen*.5,xbinlen)
	ybins = np.arange(Rmin,Rmax+ybinlen*.5,ybinlen)
fig = plt.figure()
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[50,1])
ax = fig.add_subplot(gs[1,0])
axHist = fig.add_subplot(gs[0,0], sharex=ax)
cax = fig.add_subplot(gs[1,1])
histY,histX = np.histogram(GR[good], bins=np.arange(xlim[0],xlim[1],xbinlen), range=[xlim[0],xlim[1]])#, normed=True)
histX = histX[:-1]
#histY = (5-histY)
axHist.barh([0]*len(histY),[xbinlen]*len(histY),histY,histX,color='white')
axHist.tick_params(labelbottom='off')
axHist.tick_params(labelleft='off')
axHist.tick_params(labelright='on')
#axHist.yaxis.set_major_locator(MultLoc(15000))
#axHist.yaxis.get_major_ticks()[0].label1.set_visible(False)
axHist.yaxis.set_ticks([15000,30000,45000], minor=False)

#axHist.set_ylabel("N")

if ScatterPlot:
	ax.scatter(GR[good],Kmag[good],alpha=.05, edgecolor='blue',s=.05)
else:
	data_ = ax.hist2d(GR[good],Kmag[good],bins=[xbins,ybins],cmin=1)
	data_ = data_[0][np.logical_not(np.isnan(data_[0]))]
	norm = mplcolors.Normalize(vmin=data_.max(),vmax=data_.min())
#ax[1].hexbin(GR[good],Rmag[good],mincnt=1,gridsize=1000)
print 'Number selected:', GR[good].size
print 'Out of total:', GR.size

ax.set_xlabel(xtit)
#ax[0].set_ylabel(r'$M_K$')
ax.set_ylabel(ytit)
ax.xaxis.set_minor_locator(MultLoc(xmintick))
ax.yaxis.set_minor_locator(MultLoc(ymintick))
ax.set_xlim(xlim)
ax.set_ylim(ylim)

if not ScatterPlot:
	#cax, kw = mplcolorbar.make_axes((axHist,ax), location='right', aspect=40, pad=0.02)
	cbar = mplcolorbar.ColorbarBase(cax, norm=norm, orientation='vertical', label=ztit)
	cbar.set_label(ztit, fontsize=fontsize, weight='medium')



if not 'ylim' in argdict:
	ylim0 = list(ax.get_ylim()); ylim0.reverse()
	#ylim1 = list(ax[1].get_ylim()); ylim1.reverse()
	ax.set_ylim(ylim0)
	#ax[1].set_ylim(ylim1)

plt.subplots_adjust(hspace=0, wspace=.02)
if save:
	plt.savefig(outfile)
else:
	plt.show()