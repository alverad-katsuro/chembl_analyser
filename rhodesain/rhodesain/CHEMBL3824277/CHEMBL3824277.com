%nprocshared=4
%mem=512Mb
%chk=chembl3824277.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl3824277 


0  1
N           2.81692         2.00926         1.91932
C           4.00819         1.62882         1.75753
C           5.44417         1.17023         1.56251
C           6.07252         0.34077         2.41730
C           5.48257        -0.17658         3.71837
C           4.60342        -1.44101         3.52754
N           3.89796        -2.06054         4.69684
C           3.85395        -1.49417         5.95426
O           4.41426        -0.46260         6.26259
C           3.11202        -2.25033         7.04790
N           3.31086        -3.32742         7.22650
C           3.75105        -4.20706         6.70222
C           4.59045        -5.23613         7.28309
C           4.62934        -5.43297         8.66887
C           5.39622        -6.46066         9.21588
C           6.12369        -7.30819         8.38410
C           6.09835        -7.11556         7.00531
C           5.34447        -6.07882         6.45830
C           3.15680        -4.35788         5.38502
C           2.44442        -5.53200         5.12712
C           1.83663        -5.75295         3.89801
C           1.91271        -4.78788         2.90241
C           2.57586        -3.53206         3.14123
C           3.24450        -3.33451         4.41466
H           5.93869         1.51807         0.65673
H           7.09334         0.02954         2.20906
H           4.89000         0.62583         4.16194
H           6.29560        -0.40447         4.41193
H           5.23332        -2.21572         3.08075
H           3.84291        -1.17041         2.79209
H           2.04580        -2.20408         6.80748
H           3.28514        -1.73170         7.99532
H           4.05673        -4.78908         9.33236
H           5.41924        -6.60631        10.29308
H           6.70976        -8.11896         8.81022
H           6.66261        -7.77845         6.35385
H           5.33462        -5.94603         5.37913
H           2.36450        -6.29514         5.89823
H           1.30154        -6.68124         3.71769
H           1.43610        -4.96365         1.94079
H           2.49179        -2.74085         2.40573

