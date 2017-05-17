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
from alan import prse
import numpy as np
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

def plotRotCurve(ax, R, vtheta, vtheta_err, color, label, marker):
	vtheta = np.abs(vtheta)
	ax.scatter(R,vtheta,color=color,label=label,marker=marker)
	ax.errorbar(R,vtheta,yerr=vtheta_err,color=color,ecolor=color)



if __name__ == '__main__':
	fig,ax = plt.subplots()
	infilethin = '../data/r-metabins_thindisk.csv'
	infilethick = '../data/r-metabins_thickdisk.csv'
	
	xlabel = r'$\mathregular{R}$ (kpc)'
	ylabel = r'$-\mathregular{\langle V_\theta \rangle}$ (km s$^{-1}$)'
	xlim = [8,10]
	ylim = [195,230]
	xmajtick = 0.5
	xmintick = 0.1
	ymajtick = 10
	ymintick = 1
	
	coords = prse(infilethin, delim=',', header=True, dictin=True)
	x = coords['r']; y = coords['vtheta']; err = coords['vtheta_err']
	thinlabel = r'$\mathregular{|Z| \leq 0.2}$ kpc $(N = '+str(int(np.sum(coords['count'])))+r')$'
	plotRotCurve(ax, x, y, err, 'black', thinlabel, 'x')
	
	
	coords = prse(infilethick, delim=',', header=True, dictin=True)
	x = coords['r']; y = coords['vtheta']; err = coords['vtheta_err']
	thicklabel = r'$\mathregular{|Z| > 0.2}$ kpc $(N = '+str(int(np.sum(coords['count'])))+r')$'
	plotRotCurve(ax, x, y, err, 'red', thicklabel, 's')
	
	plt.legend(loc=2)#(0.5,0.1))
	formatAxis(ax, xlim, ylim, xlabel, ylabel, xmajtick, ymajtick, xmintick, ymintick)
	
	
	if len(sys.argv)<2:
		plt.show()
	else:
		plt.savefig(sys.argv[1],bbox_inches='tight')