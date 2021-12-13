%nprocshared=4
%mem=512Mb
%chk=chembl4060936.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl4060936 


0  1
N           1.96396        -1.20479        -4.96957
C           3.01788        -1.88940        -5.07251
C           4.28831        -2.71464        -5.19659
C           4.59796        -3.63603        -4.27065
C           5.84131        -4.49834        -4.31831
C           5.48555        -5.98811        -4.52067
C           4.80297        -6.24088        -5.88266
C           4.42080        -7.69933        -6.08299
C           3.07690        -8.06798        -6.18377
C           2.72695        -9.40059        -6.39286
C           3.71842       -10.37387        -6.49788
C           5.05993       -10.01335        -6.39145
C           5.41035        -8.68019        -6.18701
N           6.59946        -4.30947        -3.07622
C           7.94860        -4.31291        -2.98766
O           8.65194        -4.52534        -3.95756
C           8.51229        -4.01614        -1.60243
C           7.97052        -2.67817        -1.02129
C           8.27691        -1.36762        -1.75471
C           9.05699        -1.28237        -2.91690
C           9.27547        -0.05023        -3.53161
C           8.73265         1.11103        -2.99016
C           7.96884         1.03938        -1.82894
C           7.74429        -0.19312        -1.21893
N           9.98124        -4.07773        -1.58498
C          10.72869        -4.21911        -0.47432
O          10.26766        -4.31351         0.64741
O          12.02713        -4.04214        -0.79168
C          13.02873        -3.98863         0.25617
C          13.78135        -5.30730         0.32486
C          15.01822        -5.43303        -0.31128
C          15.73938        -6.62222        -0.21966
C          15.22520        -7.69380         0.50702
C          13.98749        -7.57733         1.13630
C          13.26537        -6.38637         1.04547
H           4.92608        -2.51999        -6.05786
H           3.93415        -3.80300        -3.42449
H           6.46728        -4.17501        -5.15631
H           6.40169        -6.58206        -4.46431
H           4.83090        -6.32577        -3.71284
H           5.47534        -5.93348        -6.68719
H           3.90819        -5.62060        -5.96450
H           2.29440        -7.31677        -6.10669
H           1.67994        -9.68032        -6.47833
H           3.44556       -11.41315        -6.66423
H           5.83464       -10.77218        -6.47126
H           6.46053        -8.40789        -6.11199
H           6.06570        -4.13316        -2.24006
H           8.14735        -4.82529        -0.96462
H           6.88564        -2.75958        -0.92433
H           8.33862        -2.57804         0.00213
H           9.50624        -2.15965        -3.36385
H           9.87604         0.00341        -4.43640
H           8.90637         2.07116        -3.47024
H           7.54828         1.94374        -1.39749
H           7.14492        -0.23342        -0.31166
H          10.45105        -3.99941        -2.47469
H          12.56788        -3.74774         1.21817
H          13.71975        -3.17971         0.00729
H          15.42843        -4.59947        -0.87656
H          16.70432        -6.71206        -0.71204
H          15.78964        -8.62013         0.58292
H          13.58477        -8.41373         1.70256
H          12.30348        -6.29986         1.54310

