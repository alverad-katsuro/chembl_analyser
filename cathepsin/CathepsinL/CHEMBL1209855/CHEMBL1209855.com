%nprocshared=4
%mem=512Mb
%chk=chembl1209855.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1209855 


0  1
C           0.96748         0.40294         0.21783
O           2.33148         0.10051        -0.15904
C           2.99000         0.68547        -1.22909
C           2.54173         1.90506        -1.74277
C           3.18316         2.50448        -2.82100
C           4.32690         1.93121        -3.39039
C           4.95742         2.55167        -4.49048
C           6.26087         2.18130        -4.88079
C           6.79863         2.82572        -5.98535
N           7.97988         2.71568        -6.57883
C           7.83190         3.62873        -7.61656
N           6.86927         4.20907        -7.74071
C           6.12923         3.75539        -6.70455
C           4.84255         4.09092        -6.29022
C           4.03376         5.03247        -7.02099
N           3.37961         5.76432        -7.60742
N           4.28432         3.51241        -5.19478
C           4.77941         0.71648        -2.85478
C           4.11660         0.07041        -1.80449
C           4.66913        -1.30149        -1.38706
F           3.97054        -1.97493        -0.45667
F           4.74250        -2.13120        -2.44541
F           5.92372        -1.19763        -0.91131
H           0.67216        -0.26752         1.02924
H           0.30923         0.24615        -0.64139
H           0.89883         1.43764         0.56381
H           1.68262         2.40999        -1.31061
H           2.79057         3.45177        -3.18308
H           6.86652         1.46194        -4.34007
H           8.78606         2.16460        -6.31572
H           8.62535         3.83796        -8.31877
H           5.62929         0.20278        -3.29027
