%nprocshared=4
%mem=512Mb
%chk=chembl1208978.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1208978 


0  1
C          -2.08163        -0.22391         4.82335
C          -1.15966         0.54043         3.89098
C           0.22818         0.41058         3.96302
C           1.01978         1.13893         3.07657
C           0.40130         1.97351         2.14755
C          -0.99146         2.05851         2.12822
C          -1.70250         2.96481         1.12064
C          -3.23815         2.86821         1.20818
O          -3.89749         3.73857         0.25368
C          -5.28374         3.83377         0.18589
C          -6.07321         2.96378         0.94527
C          -7.45634         3.03693         0.91552
C          -8.11562         4.00963         0.15533
C          -9.52652         4.06033         0.17231
C         -10.23161         5.16865        -0.34223
C         -11.61907         5.13334        -0.25626
C         -12.30606         4.08463         0.25630
N         -13.64067         4.28384         0.21504
C         -13.70582         5.31936        -0.23057
N         -12.55837         6.00819        -0.60752
C         -12.42962         7.43201        -0.90774
C         -11.57266         3.00388         0.73902
C         -12.20743         1.83165         1.28578
N         -12.69695         0.89054         1.71084
N         -10.21731         3.00914         0.71045
C          -7.32100         4.88048        -0.60711
C          -5.92038         4.79039        -0.62894
C          -5.17956         5.75983        -1.56182
F          -3.86033         5.54955        -1.71641
F          -5.69115         5.71965        -2.80750
F          -5.31333         7.03347        -1.14867
N          -1.76199         1.34909         2.99268
H          -1.51519        -0.84569         5.51996
H          -2.74785        -0.87299         4.25073
H          -2.69450         0.46850         5.40480
H           0.68883        -0.24789         4.69680
H           2.10296         1.05599         3.11103
H           1.00100         2.55054         1.44604
H          -1.38691         2.69556         0.10982
H          -1.39831         4.00043         1.28962
H          -3.55851         3.15281         2.21519
H          -3.53605         1.83477         1.00626
H          -5.63204         2.20371         1.58246
H          -8.00530         2.33324         1.53609
H          -9.74181         6.04535        -0.75157
H         -14.68650         5.77277        -0.24199
H         -11.37693         7.72548        -0.91643
H         -12.95702         8.01145        -0.14498
H         -12.87152         7.63763        -1.88626
H          -7.78048         5.62619        -1.24695

