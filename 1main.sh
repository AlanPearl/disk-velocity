#! /bin/bash

# First, go to https://drive.google.com/open?id=0B83JqdngCoUBVDk4a2t2dVE3OGM
# and download the file datafiles.tar.gz, placing it in the data/ directory

# Uncompress data files folder
cd data
tar -xzvf datafiles.tar.gz
cd ..

# Give permission to execute my code
chmod +x *.sh
chmod +x */*.sh
chmod +x */*.py

# Make proper motion correction from LAMOST QSO/galaxy data (takes about an hour)
cd qso
./1qso.sh
cd ..

# Implement proper motion correction, and find binned velocities of disk stars
cd stars
./1stars.sh
cd ..

# Create all figures used in the paper
cd Paper
./1figures.sh
cd ..