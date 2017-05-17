#! /usr/bin/python2
"""

Creates a proper motion correction table over the specified interval, given
an `infile` containing extragalactic data.

"""

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
import pyfits as fits
from math import pi, sin, cos, sqrt, atan2
from alan import writefits, ang_dist, read_argv
from numpy import delete as np_delete, array as np_array, float64 as np_float64
import numpy as np

# Define the alan functions here

# Defining variables
infile = 'qagcoords_ltpm30_dedup.fits'
outfile = 'pearl_corr.fits'

# Default:
# 360 dec bins over range (-15,75); 1440 ra bins over range (0,360)
shape = '360,1440'
ralim = '0,360'
declim = '-15,75'
searchi = '0.5' # Search radius increment

vars().update(read_argv(sys.argv))

searchi = float(searchi)

shape = map(float, shape.split(','))
ralim = map(float, ralim.split(','))
declim = map(float, declim.split(','))

binnum_r = int(shape[0])
binnum_d = int(shape[1])

binsize = [(ralim[1]-ralim[0])/shape[0], (declim[1]-declim[0])/shape[1]]

ra0 = np.linspace(ralim[0],ralim[1],binnum_r,endpoint=False) + binsize[0]/2.
de0 = np.linspace(declim[0],declim[1],binnum_d,endpoint=False) + binsize[1]/2.

ra0,de0 = np.meshgrid(ra0,de0)
ra0 = ra0.flatten().tolist()
de0 = de0.flatten().tolist()



#Reading input files
#===========================================================

print 'Reading input files...'

# Fits keys should be 'ra', 'dec', 'pmra', and 'pmde' for:
# Right Ascension, Declination, RA Proper Motion, and DEC Proper Motion


t1 = fits.open(infile)[1].data

ra = t1['ra'].tolist()
de = t1['dec'].tolist()
pmra = t1['pmra'].tolist()
pmde = t1['pmde'].tolist()

del t1

arraydict = {}
arraydict['ra'] = ra0
arraydict['dec'] = de0


numqso = len(ra)
qsorange = xrange(numqso)


remove_index = set([])
dup_index = set([])


numqso = len(ra)
numqso1 = numqso - 1
qsorange = xrange(numqso)

numstar = len(ra0)
numstar1 = numstar - 1
starrange = xrange(numstar)

ra_rad = [None] * numqso
de_rad = [None] * numqso
ra0_rad = [None] * numstar
de0_rad = [None] * numstar

for i in qsorange:
	ra_rad[i] = ra[i] * pi/180.
	de_rad[i] = de[i] * pi/180.
for i in starrange:
	ra0_rad[i] = ra0[i] * pi/180.
	de0_rad[i] = de0[i] * pi/180.

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


#Calculate zero point proper motions for each star, using nearby quasars
#===========================================================
#================= MAIN FUNCTION ===========================
#===========================================================

pmra_fit = [0.] * numstar
pmde_fit = [0.] * numstar
se_pma = [0.] * numstar
se_pmd = [0.] * numstar

dupdict = {}
dced = []
fitfill = np_float64(0.)
errfill = np_float64(100.)

over = False
under = False
j = numqso/2 ; k = numqso/2
for i in starrange:
	if i%10000 == 0:
		#if i != 0:
		#	print 'Starting correction %d; on average, calced %.3f distances' %(
		#										i, sum(dced)/10000.)
		#	dced = []
		#else:
			print 'Starting correction %d / %d' %(i,numstar)
	
	#if i in dup_index:
	#	continue

	if under:
		j = 0


	numused = 0
	pmra_used = []
	pmde_used = []
	dist_using = []
	lastdist = 0.
	dist_calced = 0

	ra_i = ra0[i]
	de_i = de0[i]
	ra_rad_i = ra0_rad[i]
	de_rad_i = de0_rad[i]

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
	hitbottom = False
	hittop = False
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
			if not hitbottom:
				de_k = de_sort[k][0]#
				while de_k >= lde: #{LOOP 5}
					withinde.add(de_sort[k][1])###
				
					k -= 1
					if k == -1:
						hitbottom = True
						break #{LOOP 5}
				
					de_k = de_sort[k][0]#
			
				min_de = k
		
############ Add DE INDICES within new specified bounds and get new MAXIMUM {OUTSIDE}
		
			k = max_de
		
			if not hittop:
				de_k = de_sort[k][0]#
				while de_k <= ude: #{LOOP 5}
					withinde.add(de_sort[k][1])###
				
					k += 1
					if k == numqso:
						hittop = True
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
						hitbottom = True
						break #{LOOP 3}
				
					de_k = de_sort[k][0]#
			
		
			#### If initially BELOW lower bound
			else:
				while de_k < lde: #{LOOP 3.2}
					k += 1
					if k == numqso:
						hittop = True
						break #{LOOP 3.2}
				
					de_k = de_sort[k][0]#
			
				k -= 1
			
			min_de = k # Highest index below the lower bound
			k += 1
		
############ Add DE INDICES within specified bounds          {INSIDE 2x2 square}
			if not hittop:
				de_k = de_sort[k][0]#
				while de_k <= ude: #{LOOP 4}
					withinde.add(de_sort[k][1])###
				
					k += 1
					if k == numqso:
						hittop = True
						break #{LOOP 4}
				
					de_k = de_sort[k][0]
		
			max_de = k # Lowest index above the upper bound
			if not hitbottom:
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
		
		if lde < -87.:
			searchbottom = True
		elif ude > 87.:
			searchtop = True
		
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
			
			dist_i1 = ang_dist(ra_rad_i, de_rad_i, ra_rad[i1], de_rad[i1])
			dist_calced += 1
			
			dist_using.append((dist_i1, i1))
	
	
######## Use the closest quasars to get expected proper motion
	
		dist_using.sort()
		searchstop = searcha * pi / 180.
		searchbefore = (searcha - 1.01 ) * pi / 180.
	
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
		for _ in xrange(len(dist_using)): #{LOOP 2}
			
			dist_i1, i1 = dist_using[0]
			# if dist_i1 < searchbefore or dist_i1 < lastdist:
			# 	print '\n#%d' %i
			# 	print 'SEARCH WARNING: Current Search Radius of %.3f' %searcha
			# 	print 'Center Position: RA: %f\t DEC: %f' %(ra_i, de_i)
			# 	print 'Found Position:  RA: %f\t DEC: %f' %(ra[i1], de[i1])
			# 	print 'Is a distance of %f away.' %( dist_i1 * 180./pi )
			# 	print '(Last one found was %f away)' %( lastdist * 180./pi )
			# 	print '\tPREVIOUS CONSTRAINTS: (why wasn\'t it already found?)'
			# 	print '%f\t < RA  < %f' %(ra_i - rascale*(searcha-1), ra_i + rascale*(searcha-1))
			# 	print '%f\t < DEC < %f\n' %(lde+1., ude-1.)
		
# BREAK 1 (Search Radius Limit of 2.5 Degrees)
			if dist_i1 > 0.04363323129985824: # =2.5 degrees
				needmore = False
				break
			if dist_i1 > searchstop:
				break
			lastdist = dist_i1
		
			del dist_using[0]
		
			numused += 1
			floatnum = float(numused)
			floatnum1 = float(numused-1)
		
			pmra_i1 = pmra[i1]
			pmde_i1 = pmde[i1]
			pm_i1 = sqrt( pmra_i1*pmra_i1 + pmde_i1*pmde_i1 )
			pmra_used.append(pmra_i1)
			pmde_used.append(pmde_i1)
			
			using_list.append(i1)
			
# BREAK 2 (100 Quasars Found)
			if numused == 100:
				needmore = False
				break #{LOOP 2}
	
		if searcha > 2.4:
			needmore = False
	
		if not needmore:
			
			if numused > 3.:
				a = np_array(pmra_used)
				d = np_array(pmde_used)
				pmra_fit[i] = a.mean()
				pmde_fit[i] = d.mean()
				se_pma[i] = a.std() / sqrt(numused-1)
				se_pmd[i] = d.std() / sqrt(numused-1)
			else:
				pmra_fit[i] = fitfill
				pmde_fit[i] = fitfill
				se_pma[i] = errfill
				se_pmd[i] = errfill
			dced.append(dist_calced)
			break #{LOOP 1}

#===========================================================
#==================== ALL DONE :D ==========================
#===========================================================

#Make fixed (corrected) lists by subtracting each individual fit
pmra_shift = [None] * numstar
pmde_shift = [None] * numstar
for i in starrange:
	if i in remove_index:
		continue
	pmra_shift[i] = -pmra_fit[i]
	pmde_shift[i] = -pmde_fit[i]

#Write outfile
############################################################
outdict = {	'pmra_shift':	pmra_shift,
			'pmde_shift':	pmde_shift,
			'pmra_se'	:	se_pma,
			'pmde_se'	:	se_pmd,
		  }


for arrayname in outdict.keys():
	outdict[arrayname] = np_array(outdict[arrayname])

for arrayname in arraydict.keys():
	if not arrayname in outdict:
		outdict[arrayname] = arraydict[arrayname]


collen = len(outdict[outdict.keys()[0]])
if not collen == numstar:
	print 'WARNING: Column length of %d does not match original column length (%d).' %(
							collen, numstar)
maxkeylen = 0
for param in outdict:
	if len(param) > maxkeylen:
		maxkeylen = len(param)

i_p = 0
if len(remove_index) != 0:
	print 'Eliminating some stars from the columns...'

for param in outdict:
	i_p += 1
	if not len(outdict[param]) == collen:
		print 'WARNING: Uneven column lengths! len(%s) != len(%s)' %(
												param, outdict.keys()[0])


for arrayname in outdict.keys():
	x = np_delete(outdict[arrayname], list(remove_index))
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


outnames = ['ra', 'dec', 'pmra_shift', 'pmde_shift', 'pmra_se', 'pmde_se']
outarray = []
for name in outnames:
	outarray.append(np.array(outdict[name], dtype=np.float32))

print 'Writing files...'

writefits(outfile, outarray=outarray, outnames=outnames)