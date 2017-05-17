

======
PPMXL Correction
======

This is python code of the method presented in the paper "A Global Correction to PPMXL Proper Motions" from Vickers+2016.

The method recenter.recenter(ra, dec, j) returns proper motion corrections, in mas/yr (NOT deg/yr!) for the columns "pmra (Proper Motion in RA*cos(delta))" and "pmde (Proper Motion in Dec)" from the ppmxl.main database on the VO (TAP address: http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap)

Dependencies
------------

- numpy
- scipy

Tests
-----

Make sure the above dependencies are installed. In the PPMXL_Correction directory, run:

``$ python recenter.py``

This should create a new file named vc_qso_recentered.csv


``diff vc_qso_recentered.fit vc_qso_recentered_original.csv``

Should show the rubric file and the file created on your machine are identical.

Note that these test data are from Veron+2010 (located at: http://cdsarc.u-strasbg.fr/viz-bin/VizieR-3?-source=VII/258)


Another qualitative test would be to open the data in iPython and do something along the lines of:

```python
In [1]: d = np.genfromtxt('vc_qso_recentered.csv', delimiter=',', names=True)

In [2]: np.average(d['pmr_mas'])
Out[2]: -0.70694456500066327

In [3]: np.average(d['pmr_corr_mas'])
Out[3]: 0.033116644977046202

In [4]: np.average(d['pmd_mas'])
Out[4]: -2.2456262454584768

In [5]: np.average(d['pmd_corr_mas'])
Out[5]: -0.18019566627583183
```

Note on the PPMXL Database
-----
We are currently working to implement the correction across the PPMXL database.
http://dc.zah.uni-heidelberg.de/ppmxl/q/cone/info

Contact
-----
johnjvickersATshaoDOTacDOTcn
