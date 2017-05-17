#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import pyfits
from alan import prse, writefits, read_argv
import numpy as np


infile = 'spatial-bins.csv'
outfile = 'spatial-bins.fits'

args = read_argv(sys.argv)
vars().update(args)

t1 = prse(infile, delim=',', header=True, dictin=True)

t1['count'] = map(int,t1['count'])
for key in t1:
	if key!='count':
		t1[key] = np.array(t1[key], dtype=np.float32)

names=['r', 'z', 'theta', 'count', 'vr', 'vz', 'vtheta', 'vr_err', 'vz_err', 'vtheta_err', 'vr_se', 'vz_se', 'vtheta_se']
writefits(outfile, outdict=t1, nameorder=names)