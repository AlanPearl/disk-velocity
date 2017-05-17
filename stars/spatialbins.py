#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from alan import bindata, read_argv, writeASCII
from plotmodule import *
import numpy as np
import pyfits
import matplotlib.pyplot as plt

def _bindict(coords, binnames, datnames, errnames, mincount=50, verbose=True):
	datarrays = []
	for name in datnames:
		datarrays.append(coords[name])
	errarrays = []
	for name in errnames:
		errarrays.append(coords[name])
	binarrays = []
	binlens = []
	lims = []
	for name in binnames:
		binarrays.append(coords[name])
		if name.lower()=='r':
			binlens.append(Rbinlen)
			lims.append(Rlim)
		elif name.lower()=='theta':
			binlens.append(Tbinlen)
			lims.append(Tlim)
		elif name.lower()=='z':
			binlens.append(Zbinlen)
			lims.append(Zlim)
		else:
			raise ValueError("%s not understood" %name)
	
	centers,vals,errs,N = bindata(binarrays, binlens, lims, datarrays, errarrays, mincount=mincount, verbose=verbose)
	outdict = {}
	for i,name in enumerate(binnames):
		outdict[name] = centers[i]
	for i,name in enumerate(datnames):
		outdict[name] = vals[i]
	for i,name in enumerate(errnames):
		outdict[name] = errs[i]
	
	return outdict

if __name__ == '__main__':
	
	infile = 'coords_LAMOST-RAVE.fits'
	outfile = 'r-theta-z_pearlF8,10,10.csv'
	delim = ','
	binnames = 'r,theta,z'
	datnames = 'vr,vtheta,vz,vr_se,vtheta_se,vz_se'
	errnames = 'vr_err,vtheta_err,vz_err'
	lims = '8.0:10.0,9.7:-12.4,-2.0:2.0'
	binlens = '0.2,1.3,0.2'
	mincount = 10
	verbose = False
	
	# Script arguments:
	keys = ['infile','outfile','binnames','mincount']
	###
	argdict = read_argv(sys.argv, keys=keys)
	vars().update(argdict)
	binnames = binnames.split(',')
	datnames = datnames.split(',')
	errnames = errnames.split(',')
	lims = lims.split(',')
	for i in xrange(len(lims)):
		lims[i] = map(float,lims[i].split(':'))
	binlens = map(float,binlens.split(','))
	mincount = int(mincount)
	###
	
	coords = pyfits.open(infile)[1].data
	binarrays = []; datarrays = []; errarrays = [];
	for name in binnames:
		binarrays.append(coords[name])
	for name in datnames:
		datarrays.append(coords[name])
	for name in errnames:
		errarrays.append(coords[name])
	
	centers,vals,errs,count = bindata(binarrays, binlens, lims, datarrays, errarrays, mincount=mincount, verbose=verbose)
	
	names = binnames+datnames+errnames+['count']
	columns = list(centers)+list(vals)+list(errs)+[list(count)]
	writeASCII(outfile, columns, names, delim=delim)