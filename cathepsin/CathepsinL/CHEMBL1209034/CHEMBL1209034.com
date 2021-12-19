%nprocshared=4
%mem=512Mb
%chk=chembl1209034.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1209034 


0  1
C           1.49022        -0.02186        -1.87271
N           2.80519        -0.09472        -1.24129
C           3.59558         0.93849        -0.74862
N           4.48135         0.71054        -0.08245
C           4.41463        -0.63287         0.04174
C           5.20325        -1.52249         0.77077
C           6.35449        -1.07121         1.50977
N           7.28137        -0.71690         2.07769
N           4.93595        -2.85621         0.77878
C           3.91793        -3.36321         0.01873
C           3.69005        -4.75245        -0.01257
C           3.01524        -5.36166        -1.07635
C           2.88341        -6.74265        -1.13642
C           3.36858        -7.56022        -0.11146
O           3.29321        -8.94072        -0.24543
C           3.16792        -9.58878        -1.54056
C           3.35888       -11.10907        -1.40415
N           4.73031       -11.47844        -1.03292
C           5.05394       -11.75230         0.36151
C           6.58799       -11.84999         0.37042
C           6.96637       -12.11484        -1.09983
C           5.76702       -11.62488        -1.87579
O           5.71687       -11.45071        -3.07847
C           3.95269        -6.96111         1.02085
C           4.37690        -7.72290         2.28433
F           3.96863        -9.00036         2.36717
F           5.71348        -7.73419         2.43555
F           3.88088        -7.13332         3.38954
C           4.15656        -5.57504         1.02165
C           3.11094        -2.49977        -0.74741
C           3.39128        -1.14488        -0.67884
H           1.23667         1.02047        -2.08309
H           0.73864        -0.44702        -1.20133
H           1.50271        -0.58800        -2.80762
H           3.34644         1.97926        -0.89566
H           2.65533        -4.79046        -1.92547
H           2.40141        -7.16898        -2.01083
H           3.92408        -9.18705        -2.22189
H           2.17126        -9.38907        -1.94244
H           3.12484       -11.57800        -2.36419
H           2.66275       -11.49722        -0.65616
H           4.60507       -12.70503         0.65498
H           4.70233       -10.94901         1.01136
H           6.94768       -12.63447         1.03958
H           7.00855       -10.89287         0.68816
H           7.09802       -13.18549        -1.27313
H           7.87715       -11.58214        -1.38312
H           4.64946        -5.12967         1.88344
H           2.26192        -2.84095        -1.33041
