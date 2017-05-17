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
from random import random
from math import pi, sin, cos, sqrt, atan2, acos
from alan import writefits, read_argv
from numpy import delete as np_delete, array as np_array

def ang_dist(ra1,de1,ra2,de2):
	
	cosa = sin(de1)*sin(de2) + cos(de1)*cos(de2)*cos(ra2-ra1)
	if cosa >= 1.:
		return 0.
	answer = acos(cosa)
	
	return answer


# Defining variables
infile = 'qagcoords_ltpm30.fits'
outfile = 'dist_vs_match_check.fits'
pmcut = '99999999'
# Only find neighbors to this percentage of the data (randomly selected)
datapercent = '50'
overwrite = False
# ALL neighbors within this radius will be found
search0 = '1.0'
# Maximum radius to search for neighbors
searchmax = '5'
mess_int = 5000


dist_vs = []
match_vs = []

argdict = read_argv(sys.argv)
vars().update(argdict)


pmcut = float(pmcut)**2
datapercent = float(datapercent)/100.
searchmax = float(searchmax)
search0 = float(search0)
Pmin = search0 / searchmax





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


numqso = len(ra)
qsorange = xrange(numqso)


numqso = len(ra)
numqso1 = numqso - 1
qsorange = xrange(numqso)


ra_rad = [0.] * numqso
de_rad = [0.] * numqso

for i in qsorange:
	ra_rad[i] = ra[i] * pi/180.
	de_rad[i] = de[i] * pi/180.

ra_sort = [0.] * numqso 
de_sort = [0.] * numqso
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

ra_sort_under = [0.] * numqso
ra_sort_over  = [0.] * numqso
for i in qsorange:
	ra_sort_under[i] = (ra_sort[i][0]-360., ra_sort[i][1])
	ra_sort_over [i] = (ra_sort[i][0]+360., ra_sort[i][1])

#===========================================================

over = False
under = False
j = numqso/2 ; k = numqso/2
rad2deg = 180./pi

mask = np.random.rand(numqso) <= datapercent
use_qsos = np.array(range(numqso))[mask].tolist()
done = 0

#Calculate zero point proper motions for each qso, using nearby quasars
#===========================================================
#===================== MAIN CODE ===========================
#===========================================================
for i in sorted(use_qsos):
	if done % mess_int == 0:
		print 'Searching for neighbors... %d / %d' %(done,len(use_qsos))
	done += 1
	if under:
		j = 0
	
	pmra_i = pmra[i]
	pmde_i = pmde[i]
	ra_i = ra[i]
	de_i = de[i]
	ra_rad_i = ra_rad[i]
	de_rad_i = de_rad[i]
	
	# Randomly throw out more distant neighbors which are more common
	Prob = random()
	if Prob < Pmin:
		searchi = searchmax
	else:
		searchi = search0 / Prob
	
	cosdec = cos(de_rad_i)
	try:
		rascale = searchi/cosdec
	except ZeroDivisionError:
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
			
			if dist_i1 < 0.09: # Distance less than 5 degrees
				dist_vs.append(dist_i1*rad2deg)
				match_vs.append(sqrt(pow((pmra[i1]-pmra_i),2) + pow((pmde[i1]-pmde_i),2)))
		
		
			
		
		needmore = False
		break #{LOOP 1}

#===========================================================
#===========================================================
#===========================================================

#Write outfile
############################################################


print 'Writing to %s...' %outfile

dist_vs = np.array(dist_vs, dtype=np.float32)
match_vs = np.array(match_vs, dtype=np.float32)

#writefits(outfile, outdict)
writefits(outfile, {'dist':dist_vs, 'match':match_vs}, overwrite=overwrite)