%nprocshared=4
%mem=512Mb
%chk=chembl2153161.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl2153161 


0  1
C           7.96768        -0.05784        -6.89212
C           7.35520         1.18486        -6.21339
C           5.91876         1.41221        -6.73806
C           7.40988         1.00615        -4.67344
S           6.77322         2.41113        -3.78115
O           7.41792         3.57568        -4.32332
O           5.33976         2.33796        -3.82285
C           7.33280         2.10165        -2.10596
C           6.80077         3.12404        -1.06126
N           5.32967         3.09522        -1.01857
C           4.59953         2.15390        -0.38616
O           5.14044         1.23377         0.19797
C           3.10157         2.32211        -0.42462
C           2.49232         3.40198        -1.07698
C           1.10436         3.54279        -1.06667
C           0.31131         2.60577        -0.41077
C           0.90450         1.52389         0.22993
C           2.29187         1.38088         0.22281
C           7.38453         4.52539        -1.28316
O           8.39632         4.66246        -1.94225
N           6.85623         5.60092        -0.61961
C           5.85064         5.57337         0.45990
N           7.06780         6.96135        -0.98895
C           8.47911         7.39081        -1.13798
C           6.21601         7.80718        -1.70000
N           5.48211         8.58912        -2.07996
H           7.97742         0.06265        -7.97789
H           7.39207        -0.95489        -6.65045
H           8.99738        -0.20876        -6.55995
H           7.96240         2.05070        -6.48954
H           5.23964         0.64298        -6.36313
H           5.89897         1.38786        -7.82980
H           5.54158         2.38970        -6.43060
H           8.45090         0.87185        -4.36759
H           6.83363         0.12837        -4.36874
H           8.42568         2.10331        -2.10743
H           6.99125         1.10023        -1.83169
H           7.16512         2.81955        -0.07676
H           4.85367         3.82290        -1.52772
H           3.08350         4.15158        -1.59633
H           0.64229         4.38907        -1.56975
H          -0.77024         2.72067        -0.39971
H           0.28649         0.78743         0.74088
H           2.73762         0.53170         0.73174
H           5.95793         4.66680         1.06022
H           4.85030         5.59960         0.01650
H           5.96816         6.44031         1.11651
H           9.12567         6.84011        -0.45041
H           8.58063         8.46176        -0.94255
H           8.80788         7.19372        -2.16411

