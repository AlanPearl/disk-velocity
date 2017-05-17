#! /usr/bin/python2

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import pyfits as fits
import numpy as np
from time import time
from math import pi, sin, cos, sqrt, atan2
from alan import writefits, ang_dist, read_argv

# Define the alan functions here


# Defining variables
start_time = time()
infile = 'qagcoords_dedup.fits'
outfile = 'qagcoords_ltpm30.fits'
pmcut = '999999'
overwrite = False
stats = False
searchlim = '2.5'
searchi = '.15'
mess_int = '10000'


dist_vs = []
match_vs = []

argdict = read_argv(sys.argv)
vars().update(argdict)

pmcut = float(pmcut)**2
searchi = float(searchi)
mess_int = int(mess_int)
searchlim = float(searchlim) * .9999
searchlimrad = searchlim * pi / 180.






#Reading input files
#===========================================================

print 'Reading input files...'

# Fits keys should be 'ra', 'dec', 'pmra', and 'pmde' for:
# Right Ascension, Declination, RA Proper Motion, and DEC Proper Motion

t1 = fits.open(infile)[1].data
names = t1.names
arraydict = {}

for name in names:
	arraydict[name.lower()] = t1[name]


#Parameters

ra = 'ra'
de = 'dec'
pmra = 'pmra'
pmde = 'pmde'
while not arraydict.has_key(ra):
	print '\nCannot determine the QSO RA column.'
	print 'Please input its name (or "names" for options; or "quit" to quit):'
	ra = raw_input('> ')
	if ra.lower() == 'names':
		print
		print arraydict.keys()
	elif ra.lower() == 'quit':
		sys.exit()

while not arraydict.has_key(de):
	print '\nCannot determine the QSO DEC column.'
	print 'Please input its name (or "names" for options; or "quit" to quit):'
	de = raw_input('> ')
	if de.lower() == 'names':
		print
		print arraydict.keys()
	elif de.lower() == 'quit':
		sys.exit()

while not arraydict.has_key(pmra):
	print '\nCannot determine the QSO Proper Motion in RA (PMRA) column.'
	print 'Please input its name (or "names" for options; or "quit" to quit):'
	pmra = raw_input('> ')
	if pmra.lower() == 'names':
		print
		print arraydict.keys()
	elif pmra.lower() == 'quit':
		sys.exit()

while not arraydict.has_key(pmde):
	print '\nCannot determine the QSO Proper Motion in DEC (PMDE) column.'
	print 'Please input its name (or "names" for options; or "quit" to quit):'
	pmde = raw_input('> ')
	if pmde.lower() == 'names':
		print
		print arraydict.keys()
	elif pmde.lower() == 'quit':
		sys.exit()

ra = arraydict[ra].tolist()
de = arraydict[de].tolist()
pmra = arraydict[pmra].tolist()
pmde = arraydict[pmde].tolist()

del t1


numqso = len(ra)
qsorange = xrange(numqso)

#pmmaxlist = [(0., None)] * 25
#pmminlist = [(10., None)] * 25
#pmavglist = [(0., None)] * 25

remove_index = set([])
dup_index = set([])

for i in qsorange:
	
	ra[i] = float(ra[i])
	de[i] = float(de[i])
	pmra[i] = float(pmra[i])
	pmde[i] = float(pmde[i])
	
	pmsq = pmra[i]*pmra[i] + pmde[i]*pmde[i]
	if pmsq > pmcut:
		remove_index.add(i)
	

if len(remove_index) != 0:
	
	print 'Deleting %d objects with proper motions > %.2f mas/yr' %(
													len(remove_index), pmcut**.5)
	
	for i in sorted(remove_index, reverse=True):
		del ra[i], de[i], pmra[i], pmde[i]
	
	for arrayname in arraydict.keys():
		x = np.delete(arraydict[arrayname], list(remove_index))
		arraydict[arrayname] = x
	
	remove_index = set([])

#ra = tuple(ra) ; de = tuple(de) ; pmra = tuple(pmra) ; pmde = tuple(pmde)

numqso = len(ra)
numqso1 = numqso - 1
qsorange = xrange(numqso)


ra_rad = [None] * numqso
de_rad = [None] * numqso

for i in qsorange:
	ra_rad[i] = ra[i] * pi/180.
	de_rad[i] = de[i] * pi/180.

ra_sort = [None] * numqso 
de_sort = [None] * numqso
topset = set([])
bottomset = set([])
for i in qsorange:
	ra_sort[i] = (ra[i], i)
	de_sort[i] = (de[i], i)
	
	if de[i] > 87.:
		topset.add(i)
	elif de[i] < -87.:
		bottomset.add(i)

ra_sort.sort()
de_sort.sort()

ra_sort_under = [None] * numqso
ra_sort_over  = [None] * numqso
for i in qsorange:
	ra_sort_under[i] = (ra_sort[i][0]-360., ra_sort[i][1])
	ra_sort_over [i] = (ra_sort[i][0]+360., ra_sort[i][1])
	
#===========================================================


#Calculate zero point proper motions for each qso, using nearby quasars
#===========================================================
#===================== MAIN CODE ===========================
#===========================================================

pmra_fit = [0.] * numqso
pmde_fit = [0.] * numqso
pmra_se = [0.] * numqso
pmde_se = [0.] * numqso
if stats:
	numused_list = [None] * numqso
	searcha_list = [None] * numqso
	angwght_list = [None] * numqso
	magwght_list = [None] * numqso
	anglsep_list = [None] * numqso
	maxdist_list = [None] * numqso
	avgdist_list = [None] * numqso
	maxnrby_list = [None] * numqso
	avgnrby_list = [None] * numqso
dupdict = {}

fitfill = np.float64(0.)
errfill = np.float64(100.)

dced = []
over = False
under = False
j = numqso/2 ; k = numqso/2
for i in qsorange:
	if i % mess_int == 0:
		print 'Starting correction %d / %d' %(i,numqso)
	
	if i in dup_index:
		pmra_fit[i] = 12345.
		pmde_fit[i] = 12345.
		pmra_se[i] = 12345.
		pmde_se[i] = 12345.
		if stats:
			numused_list[i] = 12345.
			searcha_list[i] = 12345.
			angwght_list[i] = 12345.
			magwght_list[i] = 12345.
			anglsep_list[i] = 12345.
			maxdist_list[i] = 12345.
			avgdist_list[i] = 12345.
			maxnrby_list[i] = 12345.
			avgnrby_list[i] = 12345.
		continue
	
	if under:
		j = 0
	
	xweight = 0.
	yweight = 0.
	xmagw = 0.
	ymagw = 0.
	angwght = -1.
	magwght = 0.
	angles = []
	anglsep = 3.5
	numused = 0
	pmra_used = []
	pmde_used = []
	pm_used = []
	dist_used = []
	dist_using = []
	lastdist = 0.
	dist_calced = 0
	
	pmra_i = pmra[i]
	pmde_i = pmde[i]
	ra_i = ra[i]
	de_i = de[i]
	ra_rad_i = ra_rad[i]
	de_rad_i = de_rad[i]
	
	cosdec = cos(de_rad_i)
	try:
		rascale = searchi/cosdec
	except:
		rascale = 360.
	
	ura = ra_i+rascale; lra = ra_i-rascale
	ude = de_i+searchi; lde = de_i-searchi
	searcha = searchi
	
	if rascale > 180.:
		searchallra = True
	else:
		searchallra = False
	
	needmore = False
	
	over = False
	under = False
	minunder = False
	maxover = False
	
	searchbottom = False
	searchtop = False
	topdone = False
	bottomdone = False
	
	withinra = set([])
	withinde = set([])
	close_enough = set([])
	using_list = []
	while True: #{LOOP 1}
		
		if needmore:
    # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #   NOT FIRST TIME (SEARCHING PAST THE 2x2 SQUARE)    #		{OUTSIDE}
 # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
			
			ura += rascale; ude += searchi 
			lra -= rascale; lde -= searchi # Make bounds wider
			searcha += searchi
			
			if not searchallra and rascale > 180.:
				searchallra = True
			
			
			
############ Add DE INDICES within new specified bounds and get new MINIMUM {OUTSIDE}
			
			#k = min_de				# Should already be set
			if not searchbottom:
				de_k = de_sort[k][0]#
				while de_k >= lde: #{LOOP 5}
					withinde.add(de_sort[k][1])###
					
					k -= 1
					if k == -1:
						searchbottom = True
						break #{LOOP 5}
					
					de_k = de_sort[k][0]#
				
				min_de = k
			
############ Add DE INDICES within new specified bounds and get new MAXIMUM {OUTSIDE}
			
			k = max_de
			
			if not searchtop:
				de_k = de_sort[k][0]#
				while de_k <= ude: #{LOOP 5}
					withinde.add(de_sort[k][1])###
					
					k += 1
					if k == numqso:
						searchtop = True
						break #{LOOP 5}
					
					de_k = de_sort[k][0]#
				
				max_de = k
				
			k = min_de # Leave k at lower bound
			if k == -1:
				k = 0
			
############ Add RA INDICES within new specified bounds and get new MINIMUM {OUTSIDE}
			
			#j = min_ra				# Should already be set
			#under = minunder		# Should already be set
			
			if not searchallra:
				
				ra_j = ra_sort[j][0]#
				while ra_j >= lra:
					withinra.add(ra_sort[j][1])###
					
					j -= 1
					if j == -1:
						under = True
						j = numqso1
					
					if under:
						ra_j = ra_sort_under[j][0]#
					else:
						ra_j = ra_sort[j][0]#
				
				min_ra = j
				minunder = under
				
############ Add RA INDICES within new specified bounds and get new MAXIMUM {OUTSIDE}
				
				j = max_ra
				over = maxover
				
				ra_j = ra_sort[j][0]#
				while ra_j <= ura:
					withinra.add(ra_sort[j][1])###
					
					j += 1
					if j == numqso:
						over = True
						j = 0
					
					if over:
						ra_j = ra_sort_over[j][0]#
					else:
						ra_j = ra_sort[j][0]#
				
				max_ra = j
				maxover = over
				
				j = min_ra # Leave j at lower bound
			
		
		else:
			
			needmore = True
			
    # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #       FIRST TIME (SEARCHING INSIDE 2x2 SQUARE)      #		{INSIDE}
 # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
			
############ Find the LOWEST sorted DE index, min_de         {INSIDE 2x2 square}
			de_k = de_sort[k][0]#
			
			#### If initially ABOVE lower bound
			if de_k >= lde:
				while de_k >= lde: #{LOOP 3}
					k -= 1
					if k == -1:
						searchbottom = True
						break #{LOOP 3}
					
					de_k = de_sort[k][0]#
				
			
			#### If initially BELOW lower bound
			else:
				while de_k < lde: #{LOOP 3.2}
					k += 1
					if k == numqso:
						searchtop = True
						break #{LOOP 3.2}
					
					de_k = de_sort[k][0]#
				
				k -= 1
				
			min_de = k # Highest index below the lower bound
			k += 1
			
############ Add DE INDICES within specified bounds          {INSIDE 2x2 square}
			if not searchtop:
				de_k = de_sort[k][0]#
				while de_k <= ude: #{LOOP 4}
					withinde.add(de_sort[k][1])###
					
					k += 1
					if k == numqso:
						searchtop = True
						break #{LOOP 4}
					
					de_k = de_sort[k][0]
			
			max_de = k # Lowest index above the upper bound
			if not searchbottom:
				k = min_de # Leave k at lower bound
			else:
				k = min_de + 1
			
#									DE ABOVE
############################################################
#									RA BELOW
			
			if not searchallra:
################ Find the LOWEST sorted RA index, min_ra     {INSIDE 2x2 square}
				if lra < 0:
					j = 0
				ra_j = ra_sort[j][0]#
				
				#### If initially ABOVE lower bound
				if ra_j > lra:
					while ra_j > lra:
						j -= 1
						if j == -1:
							under = True
							j = numqso1
						
						if under:
							ra_j = ra_sort_under[j][0]#
						else:
							ra_j = ra_sort[j][0]#
					
				
				#### If initially BELOW lower bound
				else:
					while ra_j < lra:
						j += 1
						if j == numqso:
							over = True
							j = 0
						
						if over:
							ra_j = ra_sort_over[j][0]#
						else:
							ra_j = ra_sort[j][0]#
					
					j -= 1
				
				min_ra = j # Highest index below the lower bound
				minunder = under
				
################ Add RA INDICES within the specified bounds  {INSIDE 2x2 square}
				if under:
					ra_j = ra_sort_under[j][0]#
				else:
					ra_j = ra_sort[j][0]#
				
				while ra_j <= ura:
					withinra.add(ra_sort[j][1])###
					
					j += 1
					if j == numqso:
						if under:
							under = False
						else:
							over = True
						j = 0
					
					if under:
						ra_j = ra_sort_under[j][0]#
					elif over:
						ra_j = ra_sort_over[j][0]#
					else:
						ra_j = ra_sort[j][0]#
				
				
				max_ra = j # Lowest index above the upper bound
				maxover = over
				
				j = min_ra # Leave j at lower bound
				under = minunder
			
		
	############################################################
	###########  START CALCULATING DISTANCES  ##################
	############################################################
		
		if searchallra:
			if searchtop and not topdone:
				close_enough1 = withinde | topset
			elif searchbottom and not bottomdone:
				close_enough1 = withinde | bottomset
			else:
				close_enough1 = withinde
		
		else:
			if searchtop and not topdone:
				close_enough1 = (withinde & withinra) | topset
			elif searchbottom and not bottomdone:
				close_enough1 = (withinde & withinra) | bottomset
			else:
				close_enough1 = withinde & withinra
		
		close_ish = close_enough1 - close_enough
		close_enough = close_enough1 | close_enough
		
		if searchtop and not topdone:
			close_ish = close_ish
		
		
		for i1 in close_ish: #{LOOP 6}#i1 is index of close quasar
			
			if i1 == i:
				continue #{LOOP 6}
			
			dist_i1 = ang_dist(ra_rad_i, de_rad_i, ra_rad[i1], de_rad[i1])
			dist_calced += 1
			
			
			
			if dist_i1 < .0001 \
			and abs(pmra_i - pmra[i1]) < .5 \
			and abs(pmde_i - pmde[i1]) < .5:
				
				dup_index.add(i1)
				if not i in dupdict:
					dupdict[i] = [i1]
				else:
					dupdict[i].append(i1)
				
				continue #{LOOP 6}
			else:
				dist_using.append((dist_i1, i1))
				#dist_vs.append(dist_i1*180./pi)
				#match_vs.append(sqrt((pmra[i1]-pmra[i])**2 + (pmde[i1]-pmde[i])**2))
		
		
######## Use the closest quasars to get expected proper motion
		
		# Don't calculate anything. Just delete duplicates/outliers
		####################################################
		#needmore = False
		#dist_using = []
		dist_using.sort(reverse=True)
		####################################################
		
		searchstop = min( [searchlimrad, searcha * pi * .005555] )
		#searchbefore = (searcha - 1.01 ) * pi / 180.
		
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
		for _ in xrange(len(dist_using)): #{LOOP 2}
			
			dist_i1, i1 = dist_using.pop()
						
			#if dist_i1 < searchbefore or dist_i1 < lastdist:
			#	print '\n#%d' %i
			#	print 'SEARCH WARNING: Current Search Radius of %d' %searcha
			#	print 'And last distance used was:', lastdist
			#	print 'Center Position: RA: %f\t DEC: %f' %(ra_i, de_i)
			#	print 'Found Position:  RA: %f\t DEC: %f' %(ra[i1], de[i1])
			#	print 'Is a distance of %f away.' %( dist_i1 * 180./pi )
			#	print '\tPREVIOUS CONSTRAINTS: (why wasn\'t it already found?)'
			#	print '%f\t < RA  < %f' %(ra_i - rascale*(searcha-searchi), ra_i + rascale*(searcha-searchi))
			#	print '%f\t < DEC < %f\n' %(lde+1., ude-1.)
			#else:
			#	print 'Good:', dist_i1 * 180./pi, '>', searcha-1
			
			if dist_i1 > searchstop:
				dist_using.append((dist_i1, i1))
				break
			lastdist = dist_i1
			
			numused += 1
			floatnum = float(numused)
			floatnum1 = float(numused-1)
			
			pmra_i1 = pmra[i1]
			pmde_i1 = pmde[i1]
			#pm_i1 = sqrt( pmra_i1*pmra_i1 + pmde_i1*pmde_i1 )
			pmra_used.append(pmra_i1)
			pmde_used.append(pmde_i1)
			using_list.append(i1)
			
			if stats:
				pm_used.append(sqrt(pmra_i1**2 + pmde_i1**2))
				dist_used.append(dist_i1 * 180./pi)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				
				xmag = (ra[i1] - ra_i) * cosdec
				ymag =  de[i1] - de_i
				if xmag > 180.*cosdec:
					xmag -= 360.*cosdec
				elif xmag < -180.*cosdec:
					xmag += 360.*cosdec
				
				ang = atan2(xmag, ymag)
				if ang < 0.:
					ang += 2*pi
				
				angles.append(ang)
				angles.sort()
				separations = [2.*pi + angles[0] - angles[-1]]
				for anglenum in xrange(len(angles)-1):
					separation = angles[anglenum+1] - angles[anglenum]
					separations.append(separation)
				
				anglsep = max(separations)
				
				xw = cos(ang)
				yw = sin(ang)
				
				xweight = (xweight*floatnum1 + xw)/floatnum
				yweight = (yweight*floatnum1 + yw)/floatnum
				xmagw += xmag
				ymagw += ymag
				xmaga = xmagw / floatnum
				ymaga = ymagw / floatnum
				
				
				angwght = sqrt( xweight*xweight + yweight*yweight )
				magwght = sqrt( xmaga*xmaga + ymaga*ymaga )
				
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
# BREAK 2 (100 Quasars Found)
			if numused == 100:
				needmore = False
				break #{LOOP 2}
			
# BREAK 3 (Search Radius Limit of 2.5 Degrees)
		if searcha > searchlim:
			needmore = False
		
		
		if not needmore:
			
			'''#Sigma clipped average
			
			radev = sqrt( radevsq )
			dedev = sqrt( dedevsq )
			raoutu = pmra_fit_i + 2.5 * radev
			raoutl = pmra_fit_i - 2.5 * radev
			deoutu = pmde_fit_i + 2.5 * dedev
			deoutl = pmde_fit_i - 2.5 * dedev
			
			num = -1
			update_mean = False
			for i1 in using_list:
				num += 1
				if  raoutl < pmra[i1] < raoutu \
				and deoutl < pmde[i1] < deoutu:
					continue
				else:
					del using_list[num]
					update_mean = True
					numused -= 1
			
			if update_mean:
				rasum = 0
				desum = 0
				for i1 in using_list:
					rasum += pmra[i1]
					desum += pmde[i1]
				
				pmra_fit_i = rasum / float(numused)
				pmde_fit_i = desum / float(numused)
			'''
			
			
			if stats:
				numused_list[i] = numused
				searcha_list[i] = searcha
				anglsep_list[i] = anglsep*180./pi
				angwght_list[i] = angwght
				magwght_list[i] = magwght
				if numused == 0:
					avgnrby_list[i] = -1.
					maxnrby_list[i] = -1.
					avgdist_list[i] = -1.
					maxdist_list[i] = -1.
				else:
					avgnrby_list[i] = sum(pm_used) / float(numused)
					maxnrby_list[i] = max(pm_used)
					avgdist_list[i] = sum(dist_used) / float(numused)
					maxdist_list[i] = max(dist_used)
			
			if numused < 4:
				pmra_fit[i] = fitfill
				pmde_fit[i] = fitfill
				pmra_se[i] = errfill
				pmde_se[i] = errfill
			
			else:
				a = np.array(pmra_used)
				d = np.array(pmde_used)
				pmra_fit[i] = a.mean()
				pmde_fit[i] = d.mean()
				pmra_se[i] = a.std() / sqrt(numused-1)
				pmde_se[i] = d.std() / sqrt(numused-1)
			
			break #{LOOP 1}

#===========================================================
#===========================================================
#===========================================================

#for i in qsorange:
#	if i in remove_index:
#		continue
#	pmsq = pmra[i]*pmra[i] + pmde[i]*pmde[i]
#	if not pmsq > 900.:
#		remove_index.add(i)


for i in dup_index:
	remove_index.add(i)

#for i in qsorange:
#	if not i in dupdict:
#		print 'Could not find the exact match for:', i

#Make fixed (corrected) lists by subtracting each individual fit

pmra_fix = [None] * numqso
pmde_fix = [None] * numqso
for i in qsorange:
	if i in remove_index:
		pmra_fix[i] = 123456789.
		pmde_fix[i] = 123456789.
	else:
		pmra_fix[i] = pmra[i] - pmra_fit[i]
		pmde_fix[i] = pmde[i] - pmde_fit[i]
	#pmra_fix[i] = -pmra_fit[i]
	#pmde_fix[i] = -pmde_fit[i]

#Write outfile
############################################################





outdict = {	'pmra_pearl':	pmra_fix,
			'pmde_pearl':	pmde_fix,
			'pmra_se'	:	pmra_se,
			'pmde_se'	:	pmde_se,
		  }

if stats:
	outdict['numused'] = numused_list #Number of nearby QSO used in calculation
	outdict['searcha'] = searcha_list #Search radius ended on
	outdict['anglsep'] = anglsep_list #Largest angle with no nearby QSO found
	outdict['angwght'] = angwght_list #Tendency towards a direction
	outdict['magwght'] = magwght_list #Tendency towards a direction including mag
	outdict['maxdist'] = maxdist_list #Largest angle separation used in nearby QSO
	outdict['avgdist'] = avgdist_list #Average angle separation used in nearby QSO
	outdict['maxnrby'] = maxnrby_list #Largest proper motion used in nearby QSO
	outdict['avgnrby'] = avgnrby_list #Average proper motion used in nearby QSO


for arrayname in outdict.keys():
	outdict[arrayname] = np.array(outdict[arrayname])

for arrayname in arraydict.keys():
	if not arrayname in outdict:
		outdict[arrayname] = arraydict[arrayname]


if len(dupdict) != 0:
	
	ddlen = [None]*len(dupdict)
	i=-1

	for list_i in dupdict.values():
		i += 1
		ddlen[i] = len(list_i)
	
	print 'You will need to run this code again after duplicates are deleted.'
	print 'Found %d duplicates.' %len(dup_index)
	print 'This represents %d quasars, each with up to %d duplicates.\n' %(
														len(dupdict), max(ddlen) )
	
	
	dup_errors = set()
	print 'Checking duplicates...'
	for i in dupdict:
		for i1 in dupdict[i]:
			for param,name in [(ra,'ra'), (de,'de'), 
								 (pmra,'pmra'), (pmde,'pmde')]:
				if not abs(param[i] - param[i1]) < .001:
					print name, 'for', i, 'and', i1,
					print 'do not match: %f != %f' %(param[i], param[i1])
					dup_errors.add((i,i1))
	if len(dup_errors) != 0:
		print 'Incorrectly identified %d duplicates.' %len(dup_errors)
	
collen = len(outdict[outdict.keys()[0]])
if not collen == numqso:
	print 'WARNING: Column length of %d does not match original column length (%d).' %(
							collen, numqso)
maxkeylen = 0
for param in outdict:
	if len(param) > maxkeylen:
		maxkeylen = len(param)


if len(remove_index) != 0:
	print 'Eliminating duplicates and outliers from the columns...'

for param in outdict:
	if not len(outdict[param]) == collen:
		print 'WARNING: Uneven column lengths! len(%s) != len(%s)' %(
												param, outdict.keys()[0])


for arrayname in outdict.keys():
	x = np.delete(outdict[arrayname], list(remove_index))
	outdict[arrayname] = x


print 'Checking data types...'
for param in outdict:
	array = outdict[param]
	dtype = type(array[0])
	num = 0
	for i in xrange(len(array)):
		if not type(array[i]) is dtype:
			num += 1
			print 'WARNING: %s[%d] is %s, should be %s' %(param, i, 
												type(array[i]), dtype)
			if num == 100:
				print 'and so on...\n'
				break


print 'Writing to %s...' %outfile

writefits(outfile, outdict, overwrite=overwrite)


# run_time = time() - start_time
# rtmin = int(run_time // 60)
# rtsec = run_time %  60
# print '\t*** Using searchi=%.2f, that took %d min %.1f sec ***' %(
# 									searchi, rtmin, rtsec)