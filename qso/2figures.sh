#! /bin/bash
#
# Makes some additional figures of proper motion histograms

./plot.py corrname="'pearl'" save="'pmqag_hist_pearl.png'" gridline=True
./plot.py corrname="'vick'" save="'pmqag_hist_vick.png'" gridline=True
./plot.py corrname="'none'" save="'pmqag_hist_none.png'" gridline=True

./plot.py corrname="'pearl'" infile="'qagcoords_dedup.fits'" save="'pmqag_abs_hist_pearl.png'" xlim=[0,45] gridline=True plotcomps=False
./plot.py corrname="'vick'" infile="'qagcoords_dedup.fits'" save="'pmqag_abs_hist_vick.png'" xlim=[0,45] gridline=True plotcomps=False
./plot.py corrname="'none'" infile="'qagcoords_dedup.fits'" save="'pmqag_abs_hist_none.png'" xlim=[0,45] gridline=True plotcomps=False
