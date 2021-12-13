%nprocshared=4
%mem=512Mb
%chk=chembl3605414.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl3605414 


0  1
N           7.67679        -6.93123         5.16623
C           6.72429        -6.55970         4.66898
C           5.50192        -6.03800         4.03386
N           5.64086        -4.63020         3.64810
C           5.66370        -3.59910         4.52213
O           5.40192        -3.78367         5.69621
C           6.10908        -2.22882         3.97608
N           7.58645        -2.10950         3.91283
C           8.47844        -2.97593         3.37062
O           8.13950        -3.92067         2.68720
C           9.92672        -2.72072         3.68631
C          10.36211        -1.53174         4.28469
C          11.71085        -1.33530         4.57628
C          12.64716        -2.34412         4.32483
N          14.04076        -2.22481         4.65132
C          14.64192        -1.46912         5.60058
O          14.02403        -0.67888         6.28699
N          15.97511        -1.68764         5.79349
C          16.78358        -0.82124         6.67605
C          17.74379         0.05914         5.84841
C          18.65284        -0.82768         4.96793
C          17.80885        -1.79312         4.10396
C          16.83035        -2.59332         4.99100
C          12.21412        -3.52111         3.71137
F          13.09427        -4.50424         3.43966
C          10.87368        -3.70246         3.38029
C           5.58161        -1.14206         4.96639
C           5.83544         0.30498         4.48848
C           5.19093         0.53941         3.11064
C           5.71303        -0.49279         2.09244
C           5.46208        -1.93780         2.58559
H           4.66381        -6.14675         4.72728
H           5.30031        -6.62887         3.13718
H           5.97738        -4.45440         2.71329
H           7.96643        -1.36714         4.47769
H           9.66391        -0.72893         4.50495
H          12.01236        -0.38337         5.00151
H          14.61775        -2.90570         4.18316
H          17.36152        -1.45301         7.35638
H          16.14221        -0.17654         7.28278
H          18.35451         0.67159         6.51571
H          17.16487         0.73678         5.21560
H          19.32242        -1.40652         5.60921
H          19.27566        -0.20070         4.32574
H          18.46781        -2.48012         3.56771
H          17.24634        -1.22529         3.35860
H          17.40494        -3.21979         5.67925
H          16.23947        -3.26141         4.36122
H          10.56471        -4.63716         2.91775
H           6.05293        -1.27118         5.94457
H           4.50606        -1.27533         5.11467
H           6.90888         0.50352         4.43641
H           5.41277         1.00556         5.21273
H           5.41697         1.55092         2.76486
H           4.10435         0.45720         3.19496
H           6.78150        -0.33547         1.92363
H           5.21002        -0.34191         1.13441
H           5.84918        -2.63185         1.83459
H           4.38348        -2.10777         2.64630

