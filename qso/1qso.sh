#! /bin/bash

echo ">>> $ ./combine_gal_qso.py -overwrite qsofile=../data/dr3dr4q1q2_qso_PPMXL_pms_only_v2.csv galfile=../data/dr3dr4q1q2_galaxy_PPMXL_pms_only_recentered.csv outfile=../data/qagcoords.fits"
# First combine all quasars and galaxies into one file
# from Carlin's catalogs of QSOs and galaxies from LAMOST
./combine_gal_qso.py -overwrite qsofile=../data/dr3dr4q1q2_qso_PPMXL_pms_only_v2.csv galfile=../data/dr3dr4q1q2_galaxy_PPMXL_pms_only_recentered.csv outfile=../data/qagcoords.fits

echo ">>> $ ./nearestpmfit.py -overwrite infile=../data/qagcoords.fits outfile=../data/qagcoords_dedup.fits"
# Delete all duplicates
./nearestpmfit.py -overwrite infile=../data/qagcoords.fits outfile=../data/qagcoords_dedup.fits

echo ">>> $ ./nearestpmfit.py -overwrite infile=../data/qagcoords_dedup.fits outfile=../data/qagcoords_ltpm30.fits pmcut=30"
# Then delete all objects with > 30 mas/yr magnitude of proper motion and test correction
./nearestpmfit.py -overwrite infile=../data/qagcoords_dedup.fits outfile=../data/qagcoords_ltpm30.fits pmcut=30

echo ">>> $ ./mk_pearl_corr.py -overwrite infile=../data/qagcoords_ltpm30.fits outfile=../data/pearl_corr.fits shape=1440,360 ralim=0,360 declim=-15,75"
# Make the actual correction table using bins of size 0.25 x 0.25 degrees
# over RA range (0,360) and DEC range (-15,75)
# Therefore, the correction table contains 360 x 1440 = 518400 rows (each one is a bin)
./mk_pearl_corr.py -overwrite infile=../data/qagcoords_ltpm30.fits outfile=../data/pearl_corr.fits shape=1440,360 ralim=0,360 declim=-15,75


################################################################################

# Example of how to add pearl correction to file (also refer to ../stars/1stars.sh)
#../stars/add_pearl_corr.py corrfile=pearl_corr.fits infile=qagcoords_ltpm30.fits outfile=qag_pearlcorr.fits shape=1440,360 ralim=0,360 declim=-15,75

################################################################################

echo ">>> $ ./dist_vs_match.py -overwrite infile=../data/qagcoords_ltpm30.fits outfile=../data/dist_vs_match.fits datapercent=40 search0=0.8 searchmax=5"
# Make data file to find neighbor correlation in proper motion (Figure 5)
./dist_vs_match.py -overwrite infile=../data/qagcoords_ltpm30.fits outfile=../data/dist_vs_match.fits datapercent=40 search0=0.8 searchmax=5
# For greater accuracy: adjust datapercent (up to 100) and search0 (up to 5).
# Warning: It could take several hours to run/save using datapercent=100 and search0=5