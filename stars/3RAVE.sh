#! /bin/bash

# NOTE that ./1stars.sh will need to be run again after this 
# because all the spatial bin and metabin csv files are overwritten
# Assume we have already made the file LAMOST.fits with Pearl corrected proper motions

./combine_LAMOST-RAVE.py outfile=../data/LAMOST-RAVE.fits

# This time with the -RAVEonly flag (and change R limits to [7,9] instead of [8,10])
./mkcoords.py -overwrite -RAVEonly -elimsys infile=../data/LAMOST-RAVE.fits outfile=../data/coords.fits Rlim=7,9

./spatialbins.py ../data/coords.fits ../data/spatial-bins.csv lims=7.0:9.0,9.7:-12.4,-2.0:2.0

./metabins.py infile=../data/spatial-bins.csv binnames=r,theta outfile=../data/r-theta-metabins.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r,z outfile=../data/r-z-metabins.csv

./metabins.py infile=../data/spatial-bins.csv binnames=theta,z outfile=../data/theta-z-metabins.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r,theta zcut=and zgt=-.2 zlt=.2 outfile=../data/r-theta-metabins_thindisk.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=and zgt=-.2 zlt=.2 outfile=../data/r-metabins_thindisk.csv

./metabins.py infile=../data/spatial-bins.csv binnames=r zcut=or zgt=.2 zlt=-.2 outfile=../data/r-metabins_thickdisk.csv

# Make the two important figures for comparison
# This time, change R limits and theta velocity shift

# Figure 10
./sideview.py -save=sideviewRAVE.pdf rlim1=rlim2=7,9 vtheta_avg=-200

# Figure 11
./steptheta.py ../data/spatial-bins.csv -save=StepThetaRAVE.pdf ylim=7,9 vtheta_avg=-200 zlim=5.8,-4.6