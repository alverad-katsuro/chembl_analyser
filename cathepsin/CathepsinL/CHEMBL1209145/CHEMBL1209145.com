%nprocshared=4
%mem=512Mb
%chk=chembl1209145.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1209145 


0  1
C           0.76046         2.12344         2.30901
C           1.58271         3.40555         2.55624
N           2.30343         3.84564         1.33607
C           3.44165         2.97748         0.90296
C           4.78213         3.36252         1.57575
N           5.08069         4.80898         1.40838
C           5.33175         5.26668         0.00987
C           6.49100         4.56752        -0.72902
C           6.50877         5.05371        -2.19413
O           6.85546         6.46042        -2.25390
C           8.13595         6.98261        -2.13804
C           9.24592         6.13319        -2.15572
C          10.53612         6.64710        -2.13867
C          10.75597         8.02497        -2.19471
C          12.06454         8.52393        -2.24684
C          12.31905         9.78266        -2.81960
C          13.63790        10.20535        -2.85053
C          14.65916         9.49762        -2.31519
N          15.83636        10.14520        -2.44044
C          15.54119        11.08703        -2.99238
N          14.21121        11.30235        -3.33584
C          13.61855        12.37732        -4.12668
C          14.35809         8.26184        -1.74063
C          15.36046         7.46282        -1.07842
N          16.13931         6.83474        -0.52858
N          13.09270         7.77577        -1.74081
C           9.63996         8.87310        -2.15519
C           8.33338         8.37519        -2.06618
C           7.20813         9.39748        -1.86092
F           5.99504         8.89659        -1.56742
F           7.49500        10.22591        -0.83810
F           7.04811        10.18479        -2.94014
C           3.98456         5.59673         2.03076
C           2.62854         5.29864         1.34029
H           0.21220         1.85171         3.21330
H           1.39474         1.27909         2.03577
H           0.03410         2.28411         1.50808
H           2.28177         3.26068         3.38348
H           0.86932         4.18434         2.84233
H           3.56060         3.06747        -0.17924
H           3.23316         1.92750         1.11915
H           5.59222         2.75453         1.16667
H           4.72038         3.13842         2.64494
H           5.54512         6.33989         0.03137
H           4.43423         5.12318        -0.59446
H           6.34434         3.48634        -0.73592
H           7.44224         4.77619        -0.23472
H           5.50654         4.93909        -2.61751
H           7.19701         4.46442        -2.80428
H           9.13115         5.05441        -2.19103
H          11.36889         5.94783        -2.16764
H          11.54169        10.38363        -3.28542
H          16.32672        11.78128        -3.25231
H          14.38759        13.10993        -4.38822
H          12.83441        12.86579        -3.54467
H          13.18837        11.96520        -5.04505
H           9.77781         9.95065        -2.11520
H           4.20851         6.66416         1.95669
H           3.91302         5.33853         3.09150
H           2.66251         5.65420         0.30821
H           1.83625         5.86304         1.83982

