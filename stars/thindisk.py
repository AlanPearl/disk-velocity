#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from alan import prse, read_argv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as mplcolors
import matplotlib.colorbar as mplcbar
from matplotlib.patches import Ellipse
from matplotlib.ticker import MultipleLocator as MultLoc, FormatStrFormatter as StrFmt
from matplotlib.gridspec import GridSpec
from sideview import get_aspect, vvectorPlot

if __name__ == '__main__':
	infile = '../data/r-theta-metabins_thindisk.csv'
	outfile = '../Paper/thindisk.pdf'
	save = False
	flagvars = {'save':'outfile'}
	vars().update(read_argv(sys.argv, flagvars=flagvars))

	fig = plt.figure(figsize=(6.7,4))
	gs = GridSpec(1,2)
	cax = plt.subplot(gs[1]); cax.set_aspect(24.5)
	ax = plt.subplot(gs[0]); ax.set_aspect(get_aspect('T','R'))

	sys.argv = ("./vvector.py infile=%s xvar=theta yvar=r plot_errors=True xlim=9,-12"%infile).split()
	vvectorPlot(ax,cax)

	#gs.set_height_ratios([1.03,1,1.5])
	gs.set_width_ratios([40,1])
	gs.update(wspace=0.05)#, hspace=0.2)
	#plt.tight_layout()
	#plt.subplots_adjust(hspace=.3)

	if save:
		plt.savefig(outfile,bbox_inches='tight')
	else:
		plt.show()