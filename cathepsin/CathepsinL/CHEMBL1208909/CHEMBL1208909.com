%nprocshared=4
%mem=512Mb
%chk=chembl1208909.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1208909 


0  1
C           7.45197        -0.77586         1.55465
N           6.17414        -0.45304         0.92605
C           5.79025         0.74103         0.32505
N           4.89734         0.77691        -0.36816
C           4.50707        -0.51507        -0.39434
C           3.48434        -1.13920        -1.10578
C           2.59298        -0.39424        -1.95740
N           1.87131         0.19448        -2.61996
N           3.27824        -2.47858        -1.01124
C           4.03823        -3.24415        -0.17132
C           3.78868        -4.62830        -0.07444
C           3.05626        -5.28130        -1.07130
C           2.85733        -6.65582        -1.02155
C           3.31821        -7.41932         0.05285
C           3.09435        -8.92944        -0.03551
C           1.75885        -9.39177         0.58457
O           0.65961        -8.85282        -0.15198
C           3.97278        -6.76753         1.11984
C           4.44115        -7.43970         2.42231
F           3.93798        -8.65714         2.68046
F           5.78022        -7.56602         2.45160
F           4.10855        -6.70458         3.49820
C           4.23281        -5.39469         1.01208
C           5.07476        -2.64858         0.57798
C           5.27507        -1.28569         0.41065
H           8.01608         0.14271         1.73680
H           8.02998        -1.42662         0.89233
H           7.27565        -1.28718         2.50404
H           6.38121         1.64254         0.39518
H           2.68209        -4.74126        -1.93813
H           2.33581        -7.13842        -1.84560
H           3.92585        -9.47053         0.41889
H           3.10752        -9.23435        -1.08552
H           1.68965        -9.07119         1.62697
H           1.70675       -10.48359         0.55144
H           0.67126        -7.88259        -0.05951
H           4.74803        -4.91767         1.83948
H           5.74360        -3.21123         1.22094

