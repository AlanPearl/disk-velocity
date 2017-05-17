#! /usr/bin/python2
# Alan Pearl


from os import remove as file_remove
from os.path import isfile
from warnings import warn
from math import pi, atan, acos, asin, sin, cos
import numpy as np
from pyfits import Column, ColDefs, BinTableHDU

############################################################

def read_argv(argv, keys=[], flagvars={}):
	"""
	Returns a dictionary defining the values of keys given by the arguments in {argv}
		- A flag will create a boolean variable set to True
			- Example flag (before arguments): -save ---> save = True
		- Anything else will create a string variable set by the given argument
			- Example argument: "infile=data.fits" ---> infile = 'data.fits'
		- Optional `keys` gives a list of variables so the equals sign is not needed
		- Optional `flag_vars` maps a flag to another variable to take the input string
	To implement in code, next line should be >>> vars().update(argdict)
	"""
	argdict = {}
	order = 0
	for arg in argv[1:]:
		is_flag = False
		splitarg = arg.split('=')
		
		if arg[0] == '-' and not arg[1:2].isdigit():
			is_flag = True
			flag = splitarg[0].strip('-')
			argdict[flag] = True
			if len(splitarg) != 1:
				if flag in flagvars:
					key = flagvars[flag]
					argdict[key] = splitarg[-1]
				else:
					raise IOError('You must specify the variable name in' + 
						'flagvars[flag] to give the value in given argument: %s' %arg)
		
		elif len(splitarg) == 1:
			if len(keys) > order:
				key = keys[order]
				argdict[key] = arg
				order += 1
			else:
				raise IOError('Cannot do anything with given argument: %s' %arg)
		
		else:
			for key in splitarg[:-1]:
				key = key.strip('-')
				argdict[key] = splitarg[-1]
				if len(keys) > order and key == keys[order]:
					order += 1
		
		
	
	
	return argdict

############################################################

def interactive_filename(outfile='', overwrite=False):
	"""
	Returns the given outfile if that filename is not already used.
	If filename is already used, it prompts the user for a new filename to return.
	"""
	if overwrite:
		if isfile(outfile):
			file_remove(outfile)
	
	while isfile(outfile) or outfile == '':
		if not outfile == '':
			print 'The file %s already exists.' %outfile
		print 'Please enter a new filename (type "o" to overwrite or "q" to quit).'
		outfile1 = raw_input('=> ')
		if outfile1 == 'o':
			file_remove(outfile)
		elif outfile1 == 'q':
			raise KeyboardInterrupt
		elif outfile1 != '':
			outfile = outfile1
	
	return outfile

############################################################

def Tlist(inlist):
	
	nrow = len(inlist)
	if nrow == 0:
		return []
	ncol = len(inlist[0])
	
	outlist = [[None]*nrow for x in xrange(ncol)]
	
	for i in xrange(ncol):
		for j in xrange(nrow):
			outlist[i][j] = inlist[j][i]	
	
	
	return outlist

############################################################

def prse(infile, delim=' ', header=False, dictin=False, givenames=False):
	"""
	Parse an ASCII file separated by spaces and lines, 
	returning a list of lists [columns].
	(Default row format from NEMO's snapprint: x y z vx vy vz)
	Output separates components:
	outlist[0] = list of all x
	outlist[1] = list of all y
	etc...
	"""
	if not isfile(infile):
		raise ValueError('The file %s does not exist.' %(infile))
	
	bodies = open(infile).read().split('\n')
	
	if len(bodies) == 0:
		warn('The file %d is empty.' %infile)
		return
	
	while bodies[-1] == '': # delete any blank lines at the end
		bodies = bodies[:-1]
	while bodies[0] == '': # or at the beginning
		bodies = bodies[1:]
	
	if header or dictin:
		names = bodies[0].split(delim)
		bodies = bodies[1:]
	
	irange = xrange(len(bodies))
	for i in irange:
		bodies[i] = bodies[i].strip()
		bodies[i] = bodies[i].split(delim)
	
	jmax = len(bodies[0])
	for i in irange:
		for j in xrange(jmax):
			try:
				bodies[i][j] = float(bodies[i][j])
			except:
			
				bodies[i][j] = str(bodies[i][j])
	
	params = Tlist(bodies)
	
	if dictin:
		paramdict = {}
		for i,parameter in enumerate(names):
			paramdict[parameter] = params[i]
		if givenames:
			return paramdict, names
		else:
			return paramdict
	
	else:
		if givenames:
			return params, names
		else:
			return params
	
############################################################


def _writefits_helper(i,columns,outlist,name):
	try:
		dtype = outlist.dtype
	except:
		outlist = np.array(outlist)
		dtype = outlist.dtype
	
	if dtype.kind == 'S':
		dtype = '%dA' %(dtype.itemsize)
	
	columns[i] = Column(name=name, format=dtype, array=outlist)

def writefits(outfile, outdict=None, nameorder=None, outarray=None, outnames=None, overwrite=True):
	"""
	Writes a fits file from one of the following:
		1) Python dictionary `outdict` - uses keys as header names
			- Optional `nameorder` arg: list of header names in desired order
		2) Python List of columns `outarray`
			- Requires `outnames` arg: list of header names in desired order
	"""
	
	if outdict == None:
		numcol = len(outarray)
		if not numcol == len(outnames):
			raise ValueError('writefits: {outarray} and {outnames} must be same length')
		
		columns = [None] * numcol
		for i,outlist in enumerate(outarray):
			_writefits_helper(i, columns, outlist, outnames[i])
	
	else:
		columns = [None] * len(outdict)
		if nameorder == None:
			for i,name in enumerate(outdict):
				_writefits_helper(i, columns, outdict[name], name)
		else:
			for i,name in enumerate(nameorder):
				_writefits_helper(i, columns, outdict[name], name)
		
	
	cols = ColDefs(columns)
	new_table = BinTableHDU.from_columns(cols)
	
	outfile = interactive_filename(outfile, overwrite=overwrite)
	new_table.writeto(outfile)

############################################################

def bin_dens_r(infile, binsize, nbin, bodymass=1e-05):
	
	x,y,z,vx,vy,vz = prse(infile)
	r = []
	for i in xrange(len(x)):
		r_i = xyz_to_r(x[i],y[i],z[i])
		r.append(r_i)
	
	radlim = []
	for i in xrange(nbin+1):
		bin_i = i*binsize
		radlim.append(bin_i)
	
	mass = [0]*nbin
	for r_j in r:
		for i in xrange(nbin):
			if radlim[i] <= r_j < radlim[i+1]:
				mass[i] += bodymass
				break
	
	dens = []
	for i in xrange(nbin):
		volume_i = 4./3.*pi*(radlim[i+1]**3 - radlim[i]**3)
		dens_i = mass[i]/volume_i
		dens.append(dens_i)
	
	return dens

############################################################

def ang_dist(ra1, de1, ra2, de2):
	
	cosa = sin(de1)*sin(de2) + cos(de1)*cos(de2)*cos(ra1 - ra2)
	
	if cosa >= 1.:
		return 0.
	
	return acos(cosa)

############################################################

def rd2lb(ra, dec, deg=False):
	
	if np.shape(ra) == ():
		array = False
	else:
		array = True
	
	RMAT = np.array( [[-0.05487572, -0.87343729, -0.48383453],
					  [ 0.49410871, -0.44482923,  0.7469821 ],
					  [-0.86766654, -0.19807649,  0.45598456]], dtype=np.float32 )
	
	
	if deg:
		ra = np.array(ra) * np.pi / 180.
		dec = np.array(dec) * np.pi / 180.
	else:
		ra = np.array(ra)
		dec = np.array(dec)
	
	sind = np.sin(dec)
	sinr = np.sin(ra)
	cosd = np.cos(dec)
	cosr = np.cos(ra)
	
	x_j2 = cosd * cosr
	y_j2 = cosd * sinr
	z_j2 = sind
	
	xyz_j2 = np.array([x_j2,y_j2,z_j2])
	x, y, z = np.dot(RMAT, xyz_j2)
	
	l = np.arctan2(y, x)
	b = np.arcsin(z)
	if array:
		l[ l < 0. ] += 2.*np.pi
	else:
		if l < 0.: l += 2.*np.pi
	
	if deg:
		l *= 180. / np.pi
		b *= 180. / np.pi
	
	return l, b
	

############################################################

def lb2rd(l, b, deg=False):
	
	if np.shape(l) == ():
		array = False
	else:
		array = True
	
	RMAT = np.array( [[-0.05487572,  0.49410871, -0.86766654],
					  [-0.87343729, -0.44482923, -0.19807649],
					  [-0.48383453,  0.7469821,   0.45598456]], dtype=np.float32 )
	
	if deg:
		l = np.array(l) * np.pi / 180.
		b = np.array(b) * np.pi / 180.
	else:
		l = np.array(l)
		b = np.array(b)
	
	sind = np.sin(b)
	sinr = np.sin(l)
	cosd = np.cos(b)
	cosr = np.cos(l)
	
	x = cosd * cosr
	y = cosd * sinr
	z = sind
	
	xyz = np.array([x,y,z])
	x_j2, y_j2, z_j2 = np.dot(RMAT, xyz)
	
	ra = np.arctan2(y_j2, x_j2)
	dec = np.arcsin(z_j2)
	if array:
		ra[ ra < 0. ] += 2.*np.pi
	else:
		if ra < 0.: ra += 2.*np.pi
	
	if deg:
		ra *= 180. / np.pi
		dec *= 180. / np.pi
	
	return ra, dec

############################################################

def _create_edges(shape, binmap_r):
	edge_bindices = []
	for i in shape:
		edge_bindices.append((0,i-1))
	
	edges = set([])
	for i in xrange(np.product(shape)):
		candidate = binmap_r[i]
		for j,ind in enumerate(candidate):
			if ind in edge_bindices[j]:
				edges.add(i)
	
	return edges

def _create_binmap(shape, reverse=False):
	
	shape = np.array(shape).astype(int)
	
	binmap = {}
	x = np.arange(np.product(shape)).reshape(shape)
	if reverse:
		for b, i in np.ndenumerate(x):
			binmap[i] = b
	else:
		for b, i in np.ndenumerate(x):
			binmap[b] = i
	
	return binmap

def bindata(binarrays, binlens, lims, datarrays=[], errarrays=[], mincount=1, centered=True, verbose=False):
	"""\n
	Bins over the dimension of the values in the bin arrays over specified
	bin lengths and limits for each dimension. The number of objects in each bin
	is counted, and the values in the bin arrays and/or data arrays are averaged
	throughout each individual bin.
	
	Returns:
		centers - The centers of each bin dimensions from [binarrays]
			(returns average values if centered=False)
	    vals - The averages of values in the data arrays [datarrays]
	    errs - Propagated errors given in the error arrays [errarrays]
	    count - Single array containing the count inside each bin
	"""
	numdim = len(binarrays)
	numpar = len(datarrays)
	numerr = len(errarrays)
	some_par = numpar!=0
	some_err = numerr!=0
	collen = len(binarrays[0])
	binnums = [None] * numdim
	bindices = [None] * numdim
	centers = [None] * numdim
	
	for dimnum in xrange(numdim):
		xmin, xmax = map(float, lims[dimnum])
		binlen = binlens[dimnum]
		if xmin > xmax and binlen > 0:
			binlen = -binlen
		binnum = int( (xmax - xmin) / binlen )
		bins = np.linspace(xmin,xmax,binnum+1)
		centers_i = np.linspace(xmin - .5*binlen, xmax + .5*binlen, binnum+2)
		
		binner = binarrays[dimnum]
		bindex = np.digitize(binner, bins)
		
		binnums[dimnum] = binnum
		bindices[dimnum] = bindex
		centers[dimnum] = centers_i
	
	
	factors = np.zeros(numdim).astype(int)
	for i,num in enumerate(binnums):
		factors[i] = num + 2

	tot_binnum = int(np.product(factors))
	factors = factors.tolist()
	binmap = _create_binmap(factors)
	binmap_r = _create_binmap(factors, reverse=True)
	edges = _create_edges(factors, binmap_r)
	lostbins = len(edges)

	bindices = np.array(bindices)
	bindices = bindices.T
	indices = [None] * collen
	for i in xrange(collen):
		bindex = tuple(int(n) for n in bindices[i])
		indices[i] = binmap[bindex]
	
	indices = np.array(indices)
	
	
	
	binbin = []
	binpar = []
	binerr = []
	bincount = []
	delind = set([])
	lostcount = 0
	outcount = 0
	lostbin = 0
	for i in xrange(tot_binnum):
		if verbose and i % 1000 == 0 and i != 0:
			print 'Sorted %d of %d bins...' %(i, tot_binnum)
		
		bin_these = np.where( indices == i )
		count_i = len(bin_these[0])
		
		if i in edges: 
			outcount += count_i
			delind.add(i)
		elif count_i < mincount:
			lostcount += count_i
			lostbin += 1
			delind.add(i)
		
		bin_these = np.where( indices == i )
		
		binbin.append([])
		binpar.append([])
		binerr.append([])
		for array in binarrays:
			binbin[-1].append(np.array(array[bin_these]))
		for array in datarrays:
			binpar[-1].append(np.array(array[bin_these]))
		for array in errarrays:
			binerr[-1].append(np.array(array[bin_these]))
		
		bincount.append( count_i )

	bincount = np.array(bincount)
	safety = np.where( bincount > 1 )
	sqn1 = bincount * 1.
	sqn1[safety] = np.sqrt(bincount[safety] - 1.)

	binavg = [None] * tot_binnum
	paravg = [None] * tot_binnum
	erravg = [None] * tot_binnum

	where_count0 = list(np.where( bincount == 0 )[0])
	where_count1 = list(np.where( bincount == 1 )[0])
	where_countG = list(np.where( bincount >  1 )[0])

	nanpar = [np.nan] * numpar
	nanerr = [np.nan] * numerr
	bigerr = [9999.] * numerr

	for i in where_count0:
		centers_i = []
		for dimnum in xrange(numdim):
			centers_i.append(centers[dimnum][binmap_r[i][dimnum]])
		binavg[i] = centers_i
		
		paravg[i] = nanpar
		erravg[i] = nanerr

	for i in where_count1:
		if not centered:
			binavg[i] = list(np.mean(binbin[i], axis=1))
		else:
			centers_i = []
			for dimnum in xrange(numdim):
				centers_i.append(centers[dimnum][binmap_r[i][dimnum]])
			binavg[i] = centers_i
		
		if some_par:
			paravg[i] = list(np.mean(binpar[i], axis=1))
		erravg[i] = bigerr


	for i in where_countG:
		if not centered:
			binavg[i] = list(np.mean(binbin[i], axis=1))
		else:
			centers_i = []
			for dimnum in xrange(numdim):
				centers_i.append(centers[dimnum][binmap_r[i][dimnum]])
			binavg[i] = centers_i
		
		if some_par:
			paravg[i] = list(np.mean(binpar[i], axis=1))
		if some_err:
			erravg[i] = list(np.sqrt(np.sum(np.power(binerr[i],2), axis=1)) / bincount[i])
	
	delind = list(delind)
	centers = np.array(binavg).T
	vals = np.array(paravg).T
	errs = np.array(erravg).T
	
	centers = np.delete(centers, delind, axis=1)
	if some_par:
		vals = np.delete(vals, delind, axis=1)
	if some_err:
		errs = np.delete(errs, delind, axis=1)
	count = np.delete(bincount, delind)
	
	if verbose:
		if outcount != 0:
			print '\t%d stars lost outside of boundaries.' %outcount
		if lostcount != 0:
			print '\t%d stars lost in deletion of small bins.' %lostcount
		if lostbin != 0:
			print '\t%d bins deleted for having less than %d stars.' %(lostbin, mincount)
			print '\t%d bins remaining.' %len(count)
		else:
			print '\tAll %d bins have at least %d stars :)' %(len(count), mincount)
	
	return centers,vals,errs,count

############################################################

def writeASCII(filename,columns,names=None,delim=','):
	f1 = open(filename, 'w+')
	last = -len(delim)
	# Header
	if not names==None:
		assert(len(names)==len(columns))
		line = ''
		for name in names:
			line += name+delim
		f1.write(line[:last]+'\n')
	# Data
	lines = ['' for i in xrange(len(columns[0]))]
	for i in xrange(len(columns[0])):
		for col in columns:
			lines[i] += str(col[i])+delim
		lines[i] = lines[i][:last]+'\n'
	lines[-1] = lines[-1][:-1]
	for line in lines:
		f1.write(line)
	f1.close()