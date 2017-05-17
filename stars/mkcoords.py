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
import numpy as np
from alan import writefits, read_argv
from warnings import filterwarnings ; filterwarnings('ignore') # !!! NO WARNINGS !!!

def mkcoords(ra,dec,dist,dist_err,rv,rv_err,pmra,pmde,pmra_se,pmde_se,epma,epmd):
	RMAT = np.array( [[-0.05487572, -0.87343729, -0.48383453],
					  [ 0.49410871, -0.44482923,  0.7469821 ],
					  [-0.86766654, -0.19807649,  0.45598456]], dtype=np.float32 )
	
	ra_rad = ra.astype(np.float32) * (np.pi/180.)
	dec_rad = dec.astype(np.float32) * (np.pi/180.)
	
	sina = np.sin(ra_rad)
	cosa = np.cos(ra_rad)
	sind = np.sin(dec_rad)
	cosd = np.cos(dec_rad)
	
	# _j2 denotes in reference to the J2000 earth axes:
	#						with x_j2 axis at (ra,dec) = (0,0)
	#						 and y_j2 axis at (ra,dec) = (90,0)
	#						 and z_j2 axis at (ra,dec) = (0,90)
	x_j2 = dist.astype(np.float32) * cosd * cosa
	y_j2 = dist * cosd * sina
	z_j2 = dist * sind
	
	xyz_j2 = np.array([ x_j2, y_j2, z_j2 ])
	
	x_helio, y, z = np.dot(RMAT, xyz_j2)
	x = x_helio - 8. # Sun is located at (-8.0, 0.0, 0.0)_GC
	del x_helio
	R = np.sqrt(x**2 + y**2)
	
	vra =  4.741067035842384 * pmra * dist # / cosd
	vdec = 4.741067035842384 * pmde * dist
	
	vx_j2 = rv * cosd * cosa   -   vdec * sind * cosa   -   vra * sina # * cosd
	vy_j2 = rv * cosd * sina   -   vdec * sind * sina   +   vra * cosa # * cosd
	vz_j2 = rv * sind          +   vdec * cosd
	del vra, vdec
	
	vel_j2 = np.array([ vx_j2, vy_j2, vz_j2 ])
	vx, vy, vz = np.dot(RMAT, vel_j2)
	
	vx += 10.1 # Sun velocity is (10.1, 224.0, 6.7)_GC
	vy += 224.
	vz += 6.7
	
	theta = np.arctan2(y, x)
	
	theta[ theta < 0. ] += 2.*np.pi
	sint = np.sin(theta)
	cost = np.cos(theta)
	
	vR     = vx * cost + vy * sint
	vtheta = vy * cost - vx * sint
	
	
	x_j2_err = (x_j2 / dist) * (dist_err)
	y_j2_err = (y_j2 / dist) * (dist_err)
	z_j2_err = (z_j2 / dist) * (dist_err)
	
	vra_err = np.sqrt( ((4.741067035842384*dist)**2)*(epma**2) +
		((4.741067035842384*pmra)**2)*(dist_err**2) )
	vdec_err = np.sqrt( ((4.741067035842384*dist)**2)*(epmd**2) +
		((4.741067035842384*pmde)**2)*(dist_err**2) )
	
	vx_j2_err = np.sqrt( ((cosd*cosa)**2)*(rv_err**2) + 
		(sina**2)*(vra_err**2) + ((sind*cosa)**2)*(vdec_err**2) )
	vy_j2_err = np.sqrt( ((cosd*sina)**2)*(rv_err**2) + 
		(cosa**2)*(vra_err**2) + ((sind*sina)**2)*(vdec_err**2) )
	vz_j2_err = np.sqrt( (sind**2)*(rv_err**2) + (cosd**2)*(vdec_err**2) )
	del vra_err, vdec_err
	
	
	vra_se = (4.741067035842384*dist)*(pmra_se)
	vdec_se = (4.741067035842384*dist)*(pmde_se)
	
	vx_j2_se = np.sqrt( (sina**2)*(vra_se**2) + ((sind*cosa)**2)*(vdec_se**2) )
	vy_j2_se = np.sqrt( (cosa**2)*(vra_se**2) + ((sind*sina)**2)*(vdec_se**2) )
	vz_j2_se = (cosd)*(vdec_se)
	del vra_se, vdec_se
	
	xyz_err = np.array([ x_j2_err, y_j2_err, z_j2_err ])
	x_err, y_err, z_err = np.sqrt( np.dot(RMAT**2, xyz_err**2) )
	del xyz_err, x_j2_err, y_j2_err, z_j2_err
	
	vel_err = np.array([ vx_j2_err, vy_j2_err, vz_j2_err ])
	vx_err, vy_err, vz_err = np.sqrt( np.dot(RMAT**2, vel_err**2) )
	del vel_err, vx_j2_err, vy_j2_err, vz_j2_err
	
	vel_se = np.array([ vx_j2_se, vy_j2_se, vz_j2_se ])
	vx_se, vy_se, vz_se = np.sqrt( np.dot(RMAT**2, vel_se**2) )
	del vel_se, vx_j2_se, vy_j2_se, vz_j2_se
	
	
	R_err = np.sqrt( (x/R)**2*(x_err)**2 + (y/R)**2*(y_err)**2 )
	theta_err = np.sqrt( (y/(x**2+y**2))**2*(x_err)**2 + 
		(x/(x**2+y**2))**2*(y_err)**2 )
	
	vR_err = np.sqrt( (cost)**2*(vx_err)**2 + (sint)**2*(vy_err)**2 + 
		(vy*cost - vx*sint)**2*(theta_err)**2 )
	vtheta_err = np.sqrt( (sint)**2*(vx_err)**2 + (cost)**2*(vy_err)**2 + 
		(-vx*cost - vy*sint)**2*(theta_err)**2 )
	
	
	vR_se = np.sqrt( (cost)**2*(vx_se)**2 + (sint)**2*(vy_se)**2 )
	vtheta_se = np.sqrt( (sint)**2*(vx_se)**2 + (cost)**2*(vy_se)**2 )
	
	theta = 180.*(theta/np.pi - 1.)
	theta_err *= 180./np.pi
	
	return R, z, theta, R_err, z_err, theta_err, vR, vz, vtheta, vR_err, vz_err, vtheta_err, vR_se, vz_se, vtheta_se


def combinePearlVickers(pmra, pmde, pmra_se, pmde_se, pm_se_lim, pmra_pearl, pmde_pearl, pmra_vick, pmde_vick):
	use_vick_ra = pmra_se > pm_se_lim
	use_vick_de = pmde_se > pm_se_lim
	use_pearl_ra = np.logical_not(use_vick_ra)
	use_pearl_de = np.logical_not(use_vick_de)
	# Use Pearl correction where se is below the limit
	pmra[use_pearl_ra] = pmra_pearl[use_pearl_ra]
	pmde[use_pearl_de] = pmde_pearl[use_pearl_de]
	# Use Vickers correction elsewhere
	pmra[use_vick_ra] = pmra_vick[use_vick_ra]
	pmde[use_vick_de] = pmde_vick[use_vick_de]
	# Estimate that the Vickers correction has systematic error exactly at the limit
	pmra_se[use_vick_ra] = pm_se_lim
	pmde_se[use_vick_de] = pm_se_lim

def goodData(pos_errlim, vel_errlim, poslim, vellim, snrlim, R, z, theta, vR, vz, vtheta, R_err, z_err, theta_err, vR_err, vz_err, vtheta_err, snr, subclass=None, verbose=True, elimsys=False, vel_syslim=[8,10,10], vR_se=None, vz_se=None, vtheta_se=None, Fonly=False, LAMOSTonly=False, RAVEonly=False):
	Rlim, zlim, thetalim = poslim
	vRlim,vzlim,vthetalim= vellim
	good_pos = np.array([ Rlim[0] <= R, R <= Rlim[1], zlim[0] <= z , 
				z <= zlim[1], thetalim[0] <= theta, theta <= thetalim[1]])
	good_pos = np.all(good_pos, axis=0)
	good_vel = np.array([ vRlim[0] < vR, vR < vRlim[1], 
						vthetalim[0] < vtheta, vtheta < vthetalim[1], 
						vzlim[0] < vz, vz < vzlim[1] ])
	good_vel = np.all(good_vel, axis=0)
	good_snr = snr > snrlim
	good_vel_err = np.array([ vR_err < vel_errlim[0], 
						      vz_err < vel_errlim[1], 
						  vtheta_err < vel_errlim[2] ])
	good_vel_err = np.all(good_vel_err, axis=0)
	good_pos_err = np.array([ R_err		< pos_errlim[0], 
						z_err		< pos_errlim[1], 
						theta_err	< pos_errlim[2] ])
	good_pos_err = np.all(good_pos_err, axis=0)
	GOOD = np.array([good_pos, good_vel, good_pos_err, good_vel_err, good_snr])
	
	if Fonly:
		print 'F only'
		GOOD = np.concatenate([GOOD,[subclass.startswith('F')]])
	elif LAMOSTonly:
		print 'LAMOST only'
		GOOD = np.concatenate([GOOD,[subclass!='RAVE']])
	elif RAVEonly:
		print 'RAVE only'
		GOOD = np.concatenate([GOOD,[subclass=='RAVE']])
	if elimsys:
		print 'Eliminating systematic error'
		GOOD = np.concatenate([GOOD,[vR_se<vel_syslim[0],
						vz_se<vel_syslim[1],vtheta_se<vel_syslim[2]]])
	GOOD = np.all(GOOD, axis=0)

	
	if verbose:
		catalog_len = R.size
		print '\tCatalog length:', catalog_len
		print '\tNumber with bad positions:', catalog_len - np.where(good_pos)[0].size
		print '\tNumber with bad velocities:', catalog_len - np.where(good_vel)[0].size
		print '\tNumber with bad pos err:', catalog_len - np.where(good_pos_err)[0].size
		print '\tNumber with bad vel err:', catalog_len - np.where(good_vel_err)[0].size
		print '\tNew column length:', np.where(GOOD)[0].size
	
	return GOOD

if __name__ == '__main__':
	verbose = True
	overwrite = False
	# Maximum pmra_se/pmde_se value allowed to use Pearl correction
	# If above this value, then try Vickers correction instead (if not -pearl_corr_only)
	pm_se_lim = 2.0
	pos_errlim = [.2,.2,1.3] # Maximum allowed R_err, z_err, theta_err
	vel_errlim = [40,40,40] # Maximum allowed vR_err, vz_err, vtheta_err
	vel_syslim = [8,10,10] # Maximum allowed vR_se, vz_se, vtheta_se (if -elimsys)
	infile = 'obscoords_LAMOST-RAVE.fits'
	outfile = 'coords_LAMOST-RAVE.fits'
	Rlim = '8.0,10.0'
	zlim = [-2.0,2.0]
	thetalim = [-15.0,11.0]
	vRlim = [-150.0,150.0]
	vzlim = [-150.0,150.0]
	vthetalim = [-400.0,-100.0]
	snrlim = 5.0 # Minimum allowed signal to noise ratio
	elimsys = False
	LAMOSTonly = False
	RAVEonly = False
	Fonly = False
	pearl_corr_only = False # Throw out stars instead of using the Vickers correction
	
	vars().update(read_argv(sys.argv))
	if pearl_corr_only: pm_se_lim = 100.0
	Rlim = map(float, Rlim.split(','))
	
	poslim = [Rlim,zlim,thetalim]
	vellim = [vRlim, vzlim, vthetalim]
	
	input_names = ['ra', 'dec', 'dist', 'dist_err', 'rv', 'rv_err', 'pmra', 
	'pmde', 'epma', 'epmd', 'pmra_vick', 'pmde_vick', 'pmra_pearl', 'pmde_pearl', 'pmra_se', 'pmde_se', 'snr', 'M_K50', 'k0', 'j0', 'subclass']
	output_names = ['R', 'z', 'theta', 'R_err', 'z_err', 'theta_err', 'vR', 'vz', 
	'vtheta', 'vR_err', 'vz_err', 'vtheta_err', 'vR_se', 'vz_se', 'vtheta_se', 'M_K50', 
	'k0', 'j0', 'ra', 'dec']
	
	# Read file
	print 'Opening file...'
	catalog = pyfits.open(infile)[1].data
	# Ensure column names match what is expected
	if not (set(input_names).issubset(set(catalog.names))):
		print "ERROR:"
		print "Expected names:", input_names
		print "Given names:", catalog.names
		assert(False)
	# Set columns to script variables
	for name in input_names:
		vars()[name] = catalog[name]
	del catalog
	
	# Decide which Pearl correction values to throw out and use Vickers instead
	print 'Comparing Pearl and Vickers corrections...'
	combinePearlVickers(pmra, pmde, pmra_se, pmde_se, pm_se_lim, 
	pmra_pearl, pmde_pearl, pmra_vick, pmde_vick)
	
	# Transform coordinates
	print 'Beginning coordinate transformations...'
	(R, z, theta, R_err, z_err, theta_err, vR, vz, vtheta, 	vR_err, vz_err, 
	vtheta_err, vR_se, vz_se, vtheta_se) = mkcoords(ra, dec, dist, dist_err, 
	rv, rv_err, pmra, pmde, pmra_se, pmde_se, epma, epmd)
	
	# Make selection cuts
	print 'Beginning data cuts...'
	good = goodData(pos_errlim, vel_errlim, poslim, vellim, snrlim,
					R, z, theta, vR, vz, vtheta, R_err, z_err, theta_err,
					vR_err, vz_err, vtheta_err, snr, subclass, verbose, elimsys, vel_syslim, vR_se, vz_se, vtheta_se, Fonly, LAMOSTonly, RAVEonly)
	
	outdict = {}
	for name in output_names:
		outdict[name] = vars()[name][good]
		#del vars()[name]
	
	
	print 'Writing to %s...' %outfile
	writefits(outfile, outdict, output_names, overwrite=overwrite)
	print 'All done.'