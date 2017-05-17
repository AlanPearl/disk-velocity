#! /usr/bin/python2
"""

### combine_gal_qso.py
### Alan Pearl

"""

import sys,os
if '__file__' in vars():
	path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../module'))
else:
	path = os.path.abspath(os.path.join(os.getcwd(), '../module'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from alan import writefits, prse, rd2lb, read_argv
import numpy as np


qsofile = 'dr3dr4q1q2_qso_PPMXL_pms_only_v2.csv'
galfile = 'dr3dr4q1q2_galaxy_PPMXL_pms_only_recentered.csv'
outfile = 'qagcoords.fits'
overwrite = False

vars().update(read_argv(sys.argv))

qso = prse(qsofile, delim=',', dictin=True)
gal = prse(galfile, delim=',', dictin=True)

jmag = qso['jmag'] + gal['jmag']

ra = qso['ra'] + gal['ra']
dec = qso['de'] + gal['de']
pmra = qso['pmr_mas'] + gal['pmr_mas']
pmde = qso['pmd_mas'] + gal['pmd_mas']

pmra_vick = qso['pmr_corr_mas'] + gal['pmr_corr_mas']
pmde_vick = qso['pmd_corr_mas'] + gal['pmd_corr_mas']

ra = np.array(ra, dtype=np.float32)
dec = np.array(dec, dtype=np.float32)
pmra = np.array(pmra, dtype=np.float32)
pmde = np.array(pmde, dtype=np.float32)
pmra_vick = np.array(pmra_vick, dtype=np.float32)
pmde_vick = np.array(pmde_vick, dtype=np.float32)
l,b = rd2lb(ra, dec, deg=True)

outdict = {'ra':ra, 'dec':dec, 'pmra':pmra, 'pmde':pmde, 'l':l, 'b':b,  
			'pmra_vick':pmra_vick, 'pmde_vick':pmde_vick, 'jmag':jmag}

writefits(outfile, outdict, overwrite=overwrite)
