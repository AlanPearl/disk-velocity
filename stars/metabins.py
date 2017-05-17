#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from alan import prse, read_argv, writeASCII
import math
import numpy as np

def ReduceDimensions(dims,vals,errs,count,mincount=1):
	"""Use 2D lists instead of numpy arrays (count is a 1D list).
	-Arguments
	 	dims: dimensions containing repeated values to combine
		vals: values [to average within new bins]
		errs: uncertainties [to propagate within new bins]
		count: number of data in each bin [to sum within new bins]
	-Returns
		newdims
		newvals		(should be self-explanatory)
		newerrs
		newcount
	"""
	for d in dims:
		assert(len(d)==len(count))
	for v in vals:
		assert(len(v)==len(count))
	for e in errs:
		assert(len(e)==len(count))
		 
	newbindict = {}
	for i in xrange(len(dims[0])):
		loc = []
		for dimnum in xrange(len(dims)):
			loc.append(dims[dimnum][i])
		loc = tuple(loc)
		if loc in newbindict:
			newbindict[loc].append(i)
		else:
			newbindict[loc] = [i]
	newdims = sorted(newbindict.keys())
	newvals = [[] for _ in xrange(len(vals))]
	newerrs = [[] for _ in xrange(len(errs))]
	newcount = []
	for dim in newdims:
		n = len(newbindict[dim])
		newval = [0. for _ in xrange(len(vals))]
		newerr = [0. for _ in xrange(len(errs))]
		newnum = 0
		for i in newbindict[dim]:
			newnum += count[i]
			for j in xrange(len(vals)):
				newval[j] += vals[j][i]
			for j in xrange(len(errs)):
				newerr[j] += errs[j][i]*errs[j][i]
		for i in xrange(len(newval)):
			newval[i] /= n
			newvals[i].append(newval[i])
		for i in xrange(len(newerr)):
			newerr[i] = math.sqrt(newerr[i])/n
			newerrs[i].append(newerr[i])
		newcount.append(newnum)
	
	assert(len(newdims)==len(newbindict))
	assert(len(newdims[0])==len(dims))
	
	# TRANSPOSE THE DIMENSION LIST!!!
	newdimsT = [[] for _ in xrange(len(dims))]
	for dim in newdims:
		assert(len(dim)==len(dims))
		for i in xrange(len(dim)):
			newdimsT[i].append(dim[i])
	newdims = newdimsT
	
	for d in newdims:
		assert(len(d)==len(newbindict))
	for e in newerrs:
		assert(len(e)==len(newbindict))
	for v in newvals:
		assert(len(v)==len(newbindict))
	assert(len(newcount)==len(newbindict))
	
	assert(len(dims)==len(newdims))
	assert(len(vals)==len(newvals))
	assert(len(errs)==len(newerrs))
	
	delind = []
	for i,count in enumerate(newcount):
		if count < mincount:
			delind.append(i)
	delind.reverse()
	for j in delind:
		for i in xrange(len(newdims)):
			del newdims[i][j]
		for i in xrange(len(newvals)):
			del newvals[i][j]
		for i in xrange(len(newerrs)):
			del newerrs[i][j]
		del newcount[j]
	return newdims, newvals, newerrs, newcount

def CombineError(errs1,errs2):
	"""Use 2D lists instead of numpy arrays.
	-Arguments
		errs1: uncertainty lists [to add in quadrature with errs2]
		errs2: uncertainty lists [to add in quadrature with errs1]
	-Returns
		newerrs: sum in quadrature of each errs1 and errs2 list
	"""
	assert(len(errs1)==len(errs2) and len(errs1[0])==len(errs2[0]))
	newerrs = [[] for _ in xrange(len(errs1))]
	for i in xrange(len(errs1)):
		for j in xrange(len(errs1[0])):
			newerrs[i].append(math.sqrt(errs1[i][j]**2 + errs2[i][j]**2))	
	for e in newerrs:
		assert(len(e)==len(errs1[0]))
	assert(len(newerrs)==len(errs1))
	return newerrs

def makezcut(zcut,zgt,zlt,coords):
	if zcut.lower()=='no':
		return coords
	else:
		z = np.array(coords['z'])
		if zcut.lower()=='and':
			where = np.logical_and(z < zlt, z > zgt)
		elif zcut.lower()=='or':
			where = np.logical_or(z < zlt, z > zgt)
		else:
			raise IOError('Invalid value for `zcut` argument')
		for name in coords:
			coords[name] = np.array(coords[name])[where].tolist()
		return coords
		

if __name__ == '__main__':
	infile = 'spatialbins.csv'
	outfile = 'r-z-metabins.csv'
	binnames = 'r,z'
	mincount = 50
	zcut = 'no'
	zgt = zlt = '0'
	
	argdict = read_argv(sys.argv)
	vars().update(argdict)
	mincount = int(float(mincount))
	binnames = binnames.split(',')
	zgt = float(zgt)
	zlt = float(zlt)
	
	coords,names = prse(infile, delim=',', header=True, dictin=True, givenames=True)
	bins3d = ['r','theta','z']
	coords = makezcut(zcut,zgt,zlt,coords)
	
	dims = []
	for name in binnames:
		dims.append(coords[name])
	vals = []; errs1 = []; errs2 = []; senum = 0
	for name in names:
		if '_err'==name[-4:]:
			errs1.append(coords[name])
		elif '_se'==name[-3:]:
			senum += 1
			errs2.append(coords[name])
		elif 'count'==name:
			count = coords[name]
		elif not name in bins3d:
			vals.append(coords[name])
	
	errs = CombineError(errs1,errs2)
	newdims, newvals, newerrs, newcount = ReduceDimensions(dims,vals,errs,count,mincount=mincount)
	
	columns = newdims+newvals+newerrs+([[0 for i in xrange(len(newdims[0]))] for j in xrange(senum)])+[newcount]
	assert(names[-1]=='count')
	
	# Remove the contracted dimension(s) from `names`
	for name in binnames:
		bins3d.remove(name)
	for name in bins3d:
		names.remove(name)
	
	writeASCII(outfile, columns, names)
	