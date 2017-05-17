#! /usr/bin/python2
"""
This is still used for l-b histogram, but could be easily replaced 
with the function alan.bindata
"""


import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import sys
import pyfits as fits
from alan import rd2lb


#Defining variables
infile = 'coords_pearl2.4.fits'
rootname = 'theta-z-vr,vz' #Snapshot file root name
xlim = '7.8,9.8'
xbinlen = 0.1
ylim = '-2, 2'
ybinlen = 0.2
mincount = 0

keys = ['rootname=', 'infile=', 'outfile=', 'xlim=', 'ylim=', 'xbinlen=', 'ybinlen=', 'mincount=']
keycalled = {}
for key in keys: keycalled[key] = False

save = False
flags = []
args = sys.argv


#Reading arguments
#===========================================================
if len(args)>0 and 'mkbinfiles.py' in args[0]:
	args = args[1:]
else:
	args = []

if len(args) > 0:
	while len(args)>0 and args[0][0] == '-':
		flags.append(args[0][1:].lower())
		args = args[1:]
	
	for flag in flags:
		print 'Ignored flag: -' + flag
		print 'mkbinfiles.py takes no flags'
			
	i = 0
	while len(args)>0 \
	and not '=' in args[0] \
	and len(keys)>0:
		exe = "%s'%s'" %(keys[0],args[0])
		exec(exe)
		keys = keys[1:]
		args = args[1:]
		keycalled[key[i]] = True
		i += 1
	
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

table = fits.open(infile)[1].data
coords = {}
for name in table.names:
	coords[name] = table[name]
del table
coords['l'], coords['b'] = rd2lb(coords['ra'],coords['dec'])

try:
	mincount = int(mincount)
except:
	mincount = int(float(mincount))

xbinlen = float(xbinlen)
ybinlen = float(ybinlen)
xlim = xlim.split(',')
ylim = ylim.split(',')
for i in [0,1]:
	xlim[i] = float(xlim[i])
	ylim[i] = float(ylim[i])

params = rootname.split('-')

if not len(params) == 3:
	sys.exit('WARNING: %d parameters given in rootname. Please give 3 parameters.' %i)

xname,yname,zname = params

zname = zname.split(',')
zcount = len(zname)
zrange = range(zcount)

if xname == 'theta':
	if not keycalled['xlim=']:
		xlim = [11. ,-14.5]
	if not keycalled['xbinlen=']:
		xbinlen = -1.3

x = tuple(coords[xname])
y = tuple(coords[yname])
#x_err = tuple(coords[xname + '_err'])
#y_err = tuple(coords[yname + '_err'])


z = []
z_err = []
outputcount = 'False'
for q in zrange:
	qname = zname[q].lower()
	if qname == 'count':
		outputcount = q
		zcount -= 1
		zrange = range(zcount)
		continue
	z.append(tuple(coords[qname]))
	z_err.append(tuple(coords[qname + '_err']))


#Now the binning begins ... a function of x, y, z, xlim, ylim, xbinlen, ybinlen
#===========================================================
for q in zrange:
	if not len(x) == len(y) == len(z[q]):
		print 'WARNING: Lists are of unequal length.'

xmin,xmax = xlim
ymin,ymax = ylim
xwid = float(abs(xmax-xmin))
ywid = float(abs(ymax-ymin))
xcount = int(abs(xwid/xbinlen))
ycount = int(abs(ywid/ybinlen))


subsets = [0]*(xcount*ycount)
for i in range(len(subsets)):
	subsets[i] = []


lostcount = 0
for i in range(len(x)):
	if i%500000 == 0:
		print 'Searching for star', i
	
	x_i = x[i]
	y_i = y[i]
	
	if not (xmin <= x_i < xmax or xmin >= x_i > xmax):
		lostcount += 1
		continue
	
	if not (ymin <= y_i < ymax or ymin >= y_i > ymax):
		lostcount += 1
		continue
	
	xper = abs(x_i - xmin)/xwid
	yper = abs(y_i - ymin)/ywid
	
	binnum = ycount*int(xper*xcount) + int(yper*ycount)
	
	subsets[binnum].append(i)

print '%d out of %d stars were out of range' %(lostcount, len(x))


xout = []
yout = []
zout = []
count = []
for q in zrange:
	zout.append([])
setcount = 0
setcountlen = len(subsets)
binmadecount = 0
print 'Making up to %d bins...' %setcountlen
for subset in subsets:
	n = len(subset)
	if n < mincount:
		continue
	if n == 0:
		count.append(n)
		
		xindex = setcount/ycount
		yindex = setcount%ycount
		xmini = xmin + xbinlen*xindex
		xmaxi = xmin + xbinlen*(xindex+1)
		ymini = ymin + ybinlen*yindex
		ymaxi = ymin + ybinlen*(yindex+1)
		
		xout.append((xmini + xmaxi) / 2.)
		yout.append((ymini + ymaxi) / 2.)
		for q in zrange:
			zout[q].append(n)
		
		setcount += 1
		binmadecount += 1
		continue
	
	binmadecount += 1
	xtot = 0
	ytot = 0
	ztot = []
	xdenom = 0
	ydenom = 0
	zdenom = []
	for q in zrange:
		ztot.append(0)
		zdenom.append(0)
	
	for i in subset:
		#x_err_i2 = x_err[i]**2
		#y_err_i2 = y_err[i]**2
		#if x_err_i2 == 0.:
		#	x_err_i2 = 1e-17
		#if y_err_i2 == 0.:
		#	y_err_i2 = 1e-17
		xtot += x[i]#/x_err_i2
		#xdenom += 1.#/x_err_i2
		ytot += y[i]#/y_err_i2
		#ydenom += 1.#/y_err_i2
		for q in zrange:
			#z_err_i2 = z_err[q][i]**2
			#if z_err_i2 == 0.:
			#	z_err_i2 = 1e-17
			ztot[q] += z[q][i]#/z_err_i2
			#zdenom[q] += 1.#/z_err_i2
	
	floatn = float(n)
	xmean = xtot/floatn#xdenom
	ymean = ytot/floatn#ydenom
	zmean = []
	for q in zrange:
		zmean.append(ztot[q]/floatn)#zdenom[q])
	
	count.append(n)
	xout.append(xmean)
	yout.append(ymean)
	for q in zrange:
		zout[q].append(zmean[q])
	setcount += 1

print 'Made %d bins.' %binmadecount
#===========================================================

if not keycalled['outfile=']:
	outfile = []
	for qname in zname:
		outfile.append('%s-%s-%s.asc' %(xname, yname, qname))

else:
	outfile = outfile.split(',,')


while not len(outfile) == len(zname):
	print 'Unable to save. %s filenames given, while there are %s files to save (length of the last parameter in rootname).' %(len(outfile),len(zname))
	userinput = raw_input('Enter %d filename(s) (or q to quit) => ' %len(zname))
	if userinput.lower() == 'q':
		sys.exit(':(')
	else:
		outfile = userinput.split(',')

else:
	if not outputcount == 'False':
		zcount += 1
		zrange = range(zcount)
		zout = zout[:outputcount] + [count] + zout[outputcount:]
	
	for q in zrange:
		lines = []
		
		for i in range(len(xout)):
			line = '{0} {1} {2}'.format(xout[i], yout[i], zout[q][i])
			lines.append(line)
		
		f = open(outfile[q], 'w+')
		
		for line in lines:
			f.write('%s\n' %line)
		
		print 'Successfully wrote to', outfile[q]
