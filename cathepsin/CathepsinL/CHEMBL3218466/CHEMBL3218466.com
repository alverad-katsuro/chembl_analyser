%nprocshared=4
%mem=512Mb
%chk=chembl3218466.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl3218466 


0  1
C          -1.23090         0.98261        -1.69759
O          -0.20392        -0.03562        -1.71814
C           0.62010        -0.26337        -2.80614
C           0.56473         0.50986        -3.97068
C           1.40870         0.22089        -5.04174
C           2.31343        -0.83937        -4.96301
C           3.21335        -1.16825        -6.14344
C           4.51367        -0.33783        -6.13654
C           5.44126        -0.74835        -7.26246
O           6.49192        -1.31611        -7.04024
N           4.98581        -0.47502        -8.51068
C           5.62838        -0.77703        -9.81379
C           4.64976        -0.26850       -10.89339
O           3.58146         0.25137       -10.62521
N           5.01397        -0.42730       -12.18297
C           4.21966        -0.01871       -13.34081
C           4.92036        -0.33390       -14.59474
N           5.49114        -0.59249       -15.54362
C           5.83122        -2.32277        -9.96001
C           6.79271        -2.75399       -11.09378
C           8.13313        -1.99613       -11.01839
C           7.92608        -0.46899       -11.04569
C           6.98534        -0.00480        -9.90829
C           2.37561        -1.60110        -3.79481
C           1.53220        -1.31407        -2.72327
H          -1.75636         0.93392        -0.74022
H          -1.93971         0.80298        -2.51074
H          -0.76975         1.96755        -1.81017
H          -0.13295         1.33732        -4.06148
H           1.35346         0.82531        -5.94401
H           2.66823        -0.98752        -7.07324
H           3.46508        -2.23175        -6.12801
H           5.02753        -0.47574        -5.18152
H           4.27509         0.72498        -6.22763
H           4.09000        -0.01601        -8.58611
H           5.90200        -0.86449       -12.37647
H           4.04039         1.05903       -13.29544
H           3.26157        -0.54544       -13.32752
H           4.86201        -2.80732       -10.10249
H           6.24045        -2.73508        -9.03740
H           6.33438        -2.61186       -12.07295
H           6.98574        -3.82554       -11.00264
H           8.76553        -2.29184       -11.85896
H           8.65824        -2.27558       -10.10143
H           7.54835        -0.15713       -12.01996
H           8.89384         0.02323       -10.92421
H           6.79014         1.06562       -10.01303
H           7.54195        -0.12124        -8.97746
H           3.07827        -2.42738        -3.71611
H           1.58181        -1.91712        -1.81984
