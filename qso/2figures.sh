#! /bin/bash
#
# Makes bin files, using the python script, mkbinfiles.py
#
#

# Proper Motion Histograms
./plot.py corrname="'pearl'" save="'pmqag_hist_pearl.png'" gridline=True
./plot.py corrname="'vick'" save="'pmqag_hist_vick.png'" gridline=True
./plot.py corrname="'none'" save="'pmqag_hist_none.png'" gridline=True

./plot.py corrname="'pearl'" infile="'qagcoords_dedup.fits'" save="'pmqag_abs_hist_pearl.png'" xlim=[0,45] gridline=True plotcomps=False
./plot.py corrname="'vick'" infile="'qagcoords_dedup.fits'" save="'pmqag_abs_hist_vick.png'" xlim=[0,45] gridline=True plotcomps=False
./plot.py corrname="'none'" infile="'qagcoords_dedup.fits'" save="'pmqag_abs_hist_none.png'" xlim=[0,45] gridline=True plotcomps=False

# QSO proper motion plot (hexbin): qsoplot2.py
./qsoplot2.py -save=qagplot2_rd_none_ltpm30 infile=qagcoords_ltpm30.fits ra_dec=True corrname= ylim=-15,75
./qsoplot2.py -save=qagplot2_rd_vick_ltpm30 infile=qagcoords_ltpm30.fits ra_dec=True corrname=vick ylim=-15,75
./qsoplot2.py -save=qagplot2_rd_pearl_ltpm30 infile=qagcoords_ltpm30.fits ra_dec=True corrname=pearl ylim=-15,75
./qsoplot2.py -save=qagplot2_lb_none_ltpm30 infile=qagcoords_ltpm30.fits corrname=
./qsoplot2.py -save=qagplot2_lb_vick_ltpm30 infile=qagcoords_ltpm30.fits corrname=vick
./qsoplot2.py -save=qagplot2_lb_pearl_ltpm30 infile=qagcoords_ltpm30.fits corrname=pearl

# Distance vs. match
./plot.py plotnum=2 save="'dist_vs_match_ltpm30.png'" infile="'dist_vs_match_ltpm30.fits'" xlim=[0,5] gridline=True

# Subsampling plots (combined)
./subsampling_comb.py -save=subsampling_comb.png
# (individual plots)
#./subsampling.py -save=subsampling_pearl_ltpm30.png corrname=pearl
#./subsampling.py -save=subsampling_vick.png corrname=vick
#./subsampling.py -save=subsampling_none.png corrname=none





## Pearl vs. Vickers vs. None over a testing parameter (fit_err, numused, etc.)
## howfar.py : Test against AVERAGE OF MAGNITUDE OF PM to test precision
## howfar2.py: Test against MAGNITUDE OF AVERAGE OF PM to test systematic error
## save outfiles to howfar/ directory
# ./howfar.py -save=howfar/jmag_hf1.png xparam=jmag minval=10 maxval=20 xlim=10,20 ylim=5,13
# ./howfar.py -save=howfar/fit_err_hf1.0.png xparam=fit_err minval=0 maxval=10 xlim=0,6 ylim=5.5,9
# ./howfar.py -save=howfar/fit_err_hf1.1.png xparam=fit_err minval=0 maxval=2 xlim=0,2 ylim=5.5,9
# ./howfar.py -save=howfar/dec_hf1.png xparam=dec minval=-10 maxval=45 xlim=-10,45 ylim=6,9 xmajtick=15 xmintick=1
# ./howfar.py -save=howfar/b_hf1.png xparam=b minval=-45 maxval=45 xlim=-45,45 ylim=6,9 xmajtick=15 xmintick=1
# 
# ./howfar2.py -save=howfar/jmag_hf2.png xparam=jmag minval=10 maxval=20 xlim=10,20 yim=0,6
# ./howfar2.py -save=howfar/fit_err_hf2.0.png xparam=fit_err minval=0 maxval=10 xlim=0,6 ylim=0,5
# ./howfar2.py -save=howfar/fit_err_hf2.1.png xparam=fit_err minval=0 maxval=3 xlim=0,3 ylim=0,5
# ./howfar2.py -save=howfar/dec_hf2.png xparam=dec minval=-10 maxval=45 xlim=-10,45 ylim=0,5 xmajtick=15 xmintick=1
# ./howfar2.py -save=howfar/b_hf2.png xparam=b minval=-45 maxval=45 xlim=-45,45 ylim=0,5 xmajtick=15 xmintick=1