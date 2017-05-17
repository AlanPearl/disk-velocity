#! /bin/bash

# NOTE that ./1stars.sh will need to be run again after this 
# because all the spatial bin and metabin csv files are overwritten
# Assume we have already made the file LAMOST.fits with Pearl corrected proper motions

# Same as in 1stars.sh except without the -Fonly and -pearl_corr_only flags (-LAMOSTonly flag would be necessary if I had already combined the RAVE data)
./mkcoords.py -overwrite -elimsys infile=../data/LAMOST.fits outfile=../data/coords.fits

./spatialbins.py ../data/coords.fits ../data/spatial-bins.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r,theta outfile=../data/r-theta-metabins.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r,z outfile=../data/r-z-metabins.csv

./metabins.py infile=../data/spatial-bins.csv binnames=theta,z outfile=../data/theta-z-metabins.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r,theta zcut=and zgt=-.2 zlt=.2 outfile=../data/r-theta-metabins_thindisk.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=and zgt=-.2 zlt=.2 outfile=../data/r-metabins_thindisk.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=or zgt=.2 zlt=-.2 outfile=../data/r-metabins_thickdisk.csv

# Make the two important figures in this directory for comparison to F star data

# Figure 10
./sideview.py -save=sideviewLAMOST.pdf

# Figure 11
./steptheta.py ../data/spatial-bins.csv -save=StepThetaLAMOST.pdf