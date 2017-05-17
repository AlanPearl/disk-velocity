#! /bin/bash

echo ">>> $ ./mkLAMOST.py -overwrite infile=../data/dr3_dr4q1q2_stellar_gr0_ebv_2MASS-dist_shorter_dupecombine_PPMXL_UCAC4_1p5arcsec_vickers_pmcorr.fits outfile=../data/LAMOST.fits"
# Take error weighted mean of distance calculations, 
# and format header names to be compatible with mkcoords.py
./mkLAMOST.py -overwrite infile=../data/dr3_dr4q1q2_stellar_gr0_ebv_2MASS-dist_shorter_dupecombine_PPMXL_UCAC4_1p5arcsec_vickers_pmcorr.fits outfile=../data/LAMOST.fits

echo ">>> $ ./add_pearl_corr.py -overwrite infile=outfile=../data/LAMOST.fits corrfile=../data/pearl_corr.fits"
# Add columns to the star data with Pearl-corrected proper motions
./add_pearl_corr.py -overwrite infile=outfile=../data/LAMOST.fits corrfile=../data/pearl_corr.fits

echo ">>> $ ./mkcoords.py -overwrite -Fonly -elimsys -pearl_corr_only infile=../data/LAMOST.fits outfile=../data/coords.fits"
# Transform the coordinates/velocities into Galactocentric Cylindric
./mkcoords.py -overwrite -Fonly -elimsys -pearl_corr_only infile=../data/LAMOST.fits outfile=../data/coords.fits

echo ">>> $ ./spatialbins.py ../data/coords.fits ../data/spatial-bins.csv"
# Create spatial bins of the data, each 0.2 kpc x 0.2 kpc x 1.3 degrees
./spatialbins.py ../data/coords.fits ../data/spatial-bins.csv

echo ">>> $ ./metabins.py infile=../data/spatial-bins.csv binnames=r,theta outfile=../data/r-theta-metabins.csv"
# Average over z
./metabins.py infile=../data/spatial-bins.csv binnames=r,theta outfile=../data/r-theta-metabins.csv

echo ">>> $ ./metabins.py infile=../data/spatial-bins.csv binnames=r,z outfile=../data/r-z-metabins.csv"
# Average over theta
./metabins.py infile=../data/spatial-bins.csv binnames=r,z outfile=../data/r-z-metabins.csv

echo ">>> $ ./metabins.py infile=../data/spatial-bins.csv binnames=theta,z outfile=../data/theta-z-metabins.csv"
# Average over R
./metabins.py infile=../data/spatial-bins.csv binnames=theta,z outfile=../data/theta-z-metabins.csv

echo ">>> $ ./metabins.py infile=../data/spatial-bins.csv binnames=r,theta zcut=and zgt=-.2 zlt=.2 outfile=../data/r-theta-metabins_thindisk.csv"
# Average over z (thin disk only)
./metabins.py infile=../data/spatial-bins.csv binnames=r,theta zcut=and zgt=-.2 zlt=.2 outfile=../data/r-theta-metabins_thindisk.csv

echo ">>> $ ./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=and zgt=-.2 zlt=.2 outfile=../data/r-metabins_thindisk.csv"
# Average over theta and z for rotation curve (thin disk only)
./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=and zgt=-.2 zlt=.2 outfile=../data/r-metabins_thindisk.csv

echo ">>> $ ./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=or zgt=.2 zlt=-.2 outfile=../data/r-metabins_thickdisk.csv"
# Average over theta and z for rotation curve (thick disk only)
./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=or zgt=.2 zlt=-.2 outfile=../data/r-metabins_thickdisk.csv
