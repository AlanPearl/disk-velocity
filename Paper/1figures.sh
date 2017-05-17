echo "Figure 1 >>> $ cd ../qso; ./plot.py corrname="'none'" infile="'../data/qagcoords_dedup.fits'" save="'../Paper/pmqag_abs_hist_none.pdf'" xlim=[0,45] ytit="'N'" gridline=True plotcomps=False; cd ../Paper"
# Figure 1
cd ../qso; ./plot.py corrname="'none'" infile="'../data/qagcoords_dedup.fits'" save="'../Paper/pmqag_abs_hist_none.pdf'" xlim=[0,45] ytit="'N'" gridline=True plotcomps=False; cd ../Paper

echo "Figure 2 >>> $ cd ../qso; ./qsoplot2.py -save=../Paper/qagplot2_rd_none_ltpm30.pdf infile=../data/qagcoords_ltpm30.fits ra_dec=True corrname= xlim=360,0 ylim=-15,75; cd ../Paper"
# Figure 2
cd ../qso; ./qsoplot2.py -save=../Paper/qagplot2_rd_none_ltpm30.pdf infile=../data/qagcoords_ltpm30.fits ra_dec=True corrname= xlim=360,0 ylim=-15,75; cd ../Paper

echo "Figure 3 >>> $ cd ../qso; ./qsoplot2.py -save=../Paper/qagplot2_rd_vick_ltpm30.pdf infile=../data/qagcoords_ltpm30.fits ra_dec=True corrname=vick xlim=360,0 ylim=-15,75; cd ../Paper"
# Figure 3
cd ../qso; ./qsoplot2.py -save=../Paper/qagplot2_rd_vick_ltpm30.pdf infile=../data/qagcoords_ltpm30.fits ra_dec=True corrname=vick xlim=360,0 ylim=-15,75; cd ../Paper

echo "Figure 4 >>> $ cd ../qso; ./qsoplot2.py -save=../Paper/qagplot2_rd_pearl_ltpm30.pdf infile=../data/qagcoords_ltpm30.fits ra_dec=True corrname=pearl xlim=360,0 ylim=-15,75; cd ../Paper"
# Figure 4
cd ../qso; ./qsoplot2.py -save=../Paper/qagplot2_rd_pearl_ltpm30.pdf infile=../data/qagcoords_ltpm30.fits ra_dec=True corrname=pearl xlim=360,0 ylim=-15,75; cd ../Paper

echo "Figure 5 >>> $ cd ../qso; ./plot.py plotnum=2 save="'../Paper/dist_vs_match_ltpm30.pdf'" infile="'../data/dist_vs_match.fits'" xlim=[0,5] gridline=True; cd ../Paper"
# Figure 5
cd ../qso; ./plot.py plotnum=2 save="'../Paper/dist_vs_match_ltpm30.pdf'" infile="'../data/dist_vs_match.fits'" xlim=[0,5] gridline=True; cd ../Paper

echo "Figure 6 >>> $ cd ../qso; ./subsampling_comb.py -save=../Paper/subsampling_comb.pdf; cd ../Paper"
# Figure 6
cd ../qso; ./subsampling_comb.py -save=../Paper/subsampling_comb.pdf; cd ../Paper

echo "Figure 7 >>> $ cd ../qso; ./pearl_corr_map.py infile=../data/pearl_corr.fits outfile=../Paper/pearl_corr_map.pdf; cd ../Paper"
# Figure 7
cd ../qso; ./pearl_corr_map.py infile=../data/pearl_corr.fits outfile=../Paper/pearl_corr_map.pdf; cd ../Paper

echo "Figure 8 >>> $ cd ../stars; ./countplot2.py -save=../Paper/l-b-countF8,10,10.pdf rootname=l-b-count infile=../data/coords.fits xlim=275,80 ylim=-70,90 xbinnum=360 ybinnum=180 xmintick=5 ymintick=5 xmajtick=30 ymajtick=30 pixel=800,400 cmax=50; cd ../Paper"
# Figure 8
cd ../stars; ./countplot2.py -save=../Paper/l-b-countF8,10,10.pdf rootname=l-b-count infile=../data/coords.fits xlim=275,80 ylim=-70,90 xbinnum=360 ybinnum=180 xmintick=5 ymintick=5 xmajtick=30 ymajtick=30 pixel=800,400 cmax=50; cd ../Paper

echo "Figure 9 >>> $ cd ../stars; ./HRplot.py -save=../Paper/HRplot.pdf infile=../data/coords.fits -JKcolor xlim=-1.5,1 ylim=5,0 ybinlen=.05; cd ../Paper"
# Figure 9
cd ../stars; ./HRplot.py -save=../Paper/HRplot.pdf infile=../data/coords.fits -JKcolor xlim=-1.5,1 ylim=5,0 ybinlen=.05; cd ../Paper

echo "Figure 10 >>> $ cd ../stars; ./sideview.py -save=../Paper/substructureF_sideview.pdf; cd ../Paper"
# Figure 10
cd ../stars; ./sideview.py -save=../Paper/substructureF_sideview.pdf; cd ../Paper

echo "Figure 11 >>> $ cd ../stars; ./steptheta.py ../data/spatialbins.csv -save=../Paper/StepThetaF.pdf; cd ../Paper"
# Figure 11
cd ../stars; ./steptheta.py ../data/spatial-bins.csv -save=../Paper/StepThetaF.pdf; cd ../Paper

echo "Figure 12 >>> $ cd ../stars; ./rotationcurve.py ../Paper/rotationcurve.pdf; cd ../Paper"
# Figure 12
cd ../stars; ./rotationcurve.py ../Paper/rotationcurve.pdf; cd ../Paper

#===========================================================
# INCLUDED FILES:

echo ">>> $ cd ../stars; cp pearl_corr.fits ../Paper; cd ../Paper"
# pearl_corr.fits
cd ../data; cp pearl_corr.fits ../Paper; cd ../Paper

echo ">>> $ cd ../stars; ./csv2fits.py infile=../data/spatial-bins.csv outfile=../Paper/spatial-bins.fits; cd ../Paper"
# spatial-bins.fits
cd ../stars; ./csv2fits.py infile=../data/spatial-bins.csv outfile=../Paper/spatial-bins.fits; cd ../Paper
