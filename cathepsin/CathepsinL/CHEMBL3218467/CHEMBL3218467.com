%nprocshared=4
%mem=512Mb
%chk=chembl3218467.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl3218467 


0  1
C           0.33752        -0.43022         2.72511
O           1.59886        -0.66125         2.05585
C           2.29637         0.33271         1.38417
C           1.95146         1.68384         1.50348
C           2.65369         2.65593         0.79303
C           3.70483         2.28470        -0.04079
C           4.06059         0.94424        -0.15421
C           3.36806        -0.03868         0.55838
C           3.77406        -1.49695         0.39780
C           2.95831        -2.19058        -0.71297
C           3.44624        -3.59923        -0.97882
O           3.96974        -3.89919        -2.03323
N           3.24739        -4.47731         0.03329
C           3.56567        -5.92202         0.05650
C           3.21568        -6.40367         1.47233
O           2.21160        -7.04453         1.72509
N           3.96628        -5.91098         2.48574
C           3.72461        -6.25571         3.89147
C           4.61064        -5.48836         4.85921
N           5.34567        -4.85178         5.66204
C           5.07387        -6.19436        -0.25264
C           5.40819        -7.70270        -0.30272
C           4.52301        -8.43775        -1.33239
C           3.02515        -8.19936        -1.05633
C           2.70003        -6.68563        -0.99600
H          -0.01879        -1.37640         3.14125
H          -0.39211        -0.05214         2.00282
H           0.47467         0.29139         3.53453
H           1.13394         1.99982         2.14488
H           2.37776         3.70230         0.88672
H           4.24812         3.04213        -0.60197
H           4.88157         0.66274        -0.80982
H           3.64205        -2.02948         1.34225
H           4.83727        -1.55699         0.15427
H           3.03575        -1.60884        -1.63422
H           1.90161        -2.21615        -0.43158
H           2.78750        -4.12923         0.86236
H           4.74789        -5.30960         2.27477
H           3.90286        -7.32558         4.02853
H           2.68103        -6.03963         4.13650
H           5.71433        -5.70107         0.48235
H           5.34301        -5.77330        -1.22195
H           5.28173        -8.14953         0.68576
H           6.45737        -7.82775        -0.58249
H           4.73335        -9.50942        -1.29616
H           4.76534        -8.09231        -2.34047
H           2.74708        -8.69118        -0.12237
H           2.43063        -8.66578        -1.84505
H           1.63770        -6.54949        -0.77938
H           2.86898        -6.26482        -1.98979
