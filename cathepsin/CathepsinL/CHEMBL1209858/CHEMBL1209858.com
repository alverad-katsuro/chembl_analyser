%nprocshared=4
%mem=512Mb
%chk=chembl1209858.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1209858 


0  1
C           1.50336        -0.46942        -1.35961
C           2.71923        -0.18408        -0.46066
C           3.28210        -1.48299         0.14236
O           4.40950        -1.13644         0.97913
C           5.21557        -2.05128         1.63907
C           4.87693        -3.40785         1.64811
C           5.67012        -4.33743         2.31183
C           6.85833        -3.94579         2.93921
C           7.68443        -4.86653         3.60942
C           7.14662        -6.07852         4.08743
C           7.99893        -6.90924         4.80074
N           7.77789        -8.05665         5.43504
C           9.04377        -8.30258         5.95473
N           9.90176        -7.59517         5.74827
C           9.30303        -6.62980         5.01717
C           9.80778        -5.44292         4.49111
C          11.20118        -5.09785         4.61460
N          12.31085        -4.83439         4.69248
N           9.00159        -4.56788         3.83191
C           7.20952        -2.59156         2.89287
C           6.38227        -1.62887         2.30267
C           6.81230        -0.16101         2.45814
F           5.93282         0.76713         2.04453
F           7.05006         0.13284         3.75044
F           7.96078         0.08641         1.80256
H           1.12218         0.46164        -1.78495
H           0.70135        -0.93836        -0.78484
H           1.77900        -1.13490        -2.18099
H           2.42671         0.49591         0.34312
H           3.49643         0.31281        -1.04650
H           3.61659        -2.15163        -0.65618
H           2.51824        -1.97685         0.74983
H           3.98482        -3.76396         1.14277
H           5.36765        -5.37933         2.26996
H           6.09826        -6.34476         3.98617
H           6.90616        -8.52937         5.62375
H           9.23976        -9.12063         6.63353
H           8.11493        -2.24926         3.38859
