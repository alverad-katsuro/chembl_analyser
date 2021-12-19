%nprocshared=4
%mem=512Mb
%chk=chembl1209144.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1209144 


0  1
C           0.20004        -5.74982         0.54099
N          -0.15496        -4.37537         0.13274
C          -1.62469        -4.18656         0.06456
C          -1.95534        -2.75529        -0.41439
N          -1.34232        -1.73072         0.47519
C          -1.93866        -1.63643         1.83715
C          -3.43539        -1.24983         1.83050
C          -3.97248        -1.11529         3.26640
O          -3.82013        -2.39644         3.92297
C          -4.19496        -2.67823         5.22434
C          -4.74988        -1.68009         6.03066
C          -5.13581        -1.95475         7.33748
C          -5.05070        -3.25354         7.84995
C          -5.46227        -3.52412         9.16777
C          -5.81573        -4.83414         9.54699
C          -6.20237        -5.02731        10.86438
C          -6.20153        -4.04116        11.79045
N          -6.57819        -4.49984        13.00343
C          -6.79266        -5.59097        12.79818
N          -6.60620        -6.10916        11.52048
C          -6.98853        -7.41525        10.99069
C          -5.83094        -2.76155        11.37896
C          -5.72125        -1.66740        12.30510
N          -5.61022        -0.79991        13.04160
N          -5.49928        -2.50948        10.08301
C          -4.49683        -4.24875         7.03366
C          -4.02023        -3.97184         5.74803
C          -3.30073        -5.11746         5.01862
F          -2.71088        -4.81197         3.84821
F          -2.31584        -5.62732         5.78324
F          -4.13728        -6.13954         4.76369
C           0.13092        -1.92487         0.52707
C           0.50266        -3.35599         0.98903
H          -0.25821        -6.47102        -0.14122
H          -0.15449        -5.94586         1.55744
H           1.28537        -5.87950         0.50965
H          -2.07099        -4.36109         1.04802
H          -2.05464        -4.90445        -0.63918
H          -3.03708        -2.62517        -0.47656
H          -1.55411        -2.62777        -1.42385
H          -1.39746        -0.87151         2.40100
H          -1.81663        -2.58551         2.36500
H          -4.02722        -2.01043         1.31751
H          -3.57194        -0.30547         1.29906
H          -5.03268        -0.84624         3.24106
H          -3.40738        -0.35315         3.81056
H          -4.88475        -0.66903         5.65630
H          -5.55969        -1.14863         7.93070
H          -5.84995        -5.66271         8.84590
H          -7.18772        -6.18171        13.61170
H          -7.34664        -8.05321        11.80267
H          -6.12157        -7.88388        10.51732
H          -7.78462        -7.29298        10.25063
H          -4.34622        -5.25302         7.41993
H           0.58190        -1.19223         1.20184
H           0.54476        -1.76246        -0.47224
H           0.19869        -3.49699         2.03011
H           1.58796        -3.47639         0.93082
