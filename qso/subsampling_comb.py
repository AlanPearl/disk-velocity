#! /usr/bin/python2
"""

### subsampling_comb.py
### Alan Pearl

Combines the subsampling plots. Hard coded

Flags: (must come before keywords)
======
	-save[=outfile]

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
#import matplotlib as mpl
from scipy.optimize import curve_fit

# Define the alan functions here

# Defining variables
outfile = 'subsampling_comb.png'

# Trying to subtract out random error in subsampling2.py
args_pearl = 'subsampling.py -keep infile=../data/qagcoords_ltpm30.fits corrname=pearl'.split(' ')
args_vick = 'subsampling.py -keep infile=../data/qagcoords_ltpm30.fits corrname=vick'.split(' ')
args_none = 'subsampling.py -keep infile=../data/qagcoords_ltpm30.fits corrname=none'.split(' ')

xtit = 'Sampling Interval (degrees)'
ytit = r'Mean Magnitude of Proper Motion (mas yr$^{-1}$)'
xlim = '0,45'
ylim = '0,8'
save = False

keys = ['xlim', 'ylim']
flags = ['save']
flagvars = {'save':'outfile'}
argdict = read_argv(sys.argv, keys, flagvars)
vars().update(argdict)

xlim = map(float, xlim.split(','))
ylim = map(float, ylim.split(','))

############################################################

pearl = {}
vick = {}
none = {}
print 'Starting...'
sys.argv = args_pearl
execfile(sys.argv[0], pearl)
print '1/3...'
sys.argv = args_vick
execfile(sys.argv[0], vick)
print '2/3...'
sys.argv = args_none
execfile(sys.argv[0], none)
print '3/3...'
spaces = [pearl, vick, none]

############################################################

fig, ax = plt.subplots()
i = 0
colors = ['blue', 'red', 'black']
markers = ['o', 'd', 'x']
for space in spaces:
	corrname = space['corrname']
	if corrname == 'pearl':
		labelname = 'Pearl'
	elif corrname == 'vick':
		labelname = 'Vickers'
	else:
		labelname = 'None'

	#ax.plot(bfline_x, bfline_y, linestyle='--', linewidth=1, color='black')
	ax.errorbar(space['rates'], space['means'], ecolor=colors[i],
				yerr=space['means_err'], linestyle='none')
	ax.scatter(space['rates'], space['means'], label=labelname, 
				marker=markers[i], color=colors[i])
	
	# Print averages for numerical comparison
	#ax.text(20, 7-i/1.5, '%s avg: %.3f' %(labelname,
	#									space['avg_val']))
	
	i += 1

ax.set_xlabel(xtit)
ax.set_ylabel(ytit)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.legend()

if save:
	plt.savefig(outfile,bbox_inches='tight')
else:
	plt.show()




