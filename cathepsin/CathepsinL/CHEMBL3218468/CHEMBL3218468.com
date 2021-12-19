%nprocshared=4
%mem=512Mb
%chk=chembl3218468.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl3218468 


0  1
C          -0.01802        -0.53717        -0.13995
O          -0.39905        -1.91395        -0.36188
C           0.51037        -2.96059        -0.42419
C           1.76918        -2.84175         0.17595
C           2.67512        -3.89900         0.14099
C           2.33689        -5.09301        -0.49233
C           3.31438        -6.25355        -0.51721
C           4.56472        -5.95982        -1.37235
C           5.54119        -7.11041        -1.28639
O           6.55826        -7.01711        -0.63110
N           5.12641        -8.25964        -1.88013
C           5.82276        -9.56603        -1.88877
C           4.91261       -10.56368        -2.63819
O           5.26661       -11.17864        -3.62508
N           3.63263       -10.66142        -2.19491
C           2.58592       -11.50969        -2.80021
C           1.76113       -10.73259        -3.81325
N           1.07691       -10.08792        -4.65365
C           6.03379       -10.08811        -0.43037
C           6.94733       -11.32945        -0.36521
C           8.31690       -11.01924        -0.99855
C           8.14192       -10.62140        -2.47610
C           7.19871        -9.40402        -2.61872
C           1.08135        -5.21973        -1.08500
C           0.16804        -4.16785        -1.06184
O          -1.04288        -4.41063        -1.69301
C          -1.85268        -3.38529        -2.31026
H          -0.91062         0.08923        -0.21577
H           0.70572        -0.23576        -0.90216
H           0.41486        -0.42681         0.85768
H           2.06025        -1.92983         0.68842
H           3.64671        -3.78691         0.61489
H           2.81676        -7.14449        -0.90657
H           3.62196        -6.48508         0.50498
H           5.05595        -5.05152        -1.01425
H           4.27959        -5.78961        -2.41362
H           4.26242        -8.23012        -2.39871
H           3.36053       -10.10267        -1.40117
H           1.92217       -11.87184        -2.01131
H           3.04344       -12.36886        -3.29855
H           5.07361       -10.31953         0.03639
H           6.49300        -9.31939         0.19120
H           6.47936       -12.17239        -0.87942
H           7.08825       -11.61947         0.67879
H           8.96570       -11.89502        -0.92711
H           8.80305       -10.20638        -0.45144
H           7.75787       -11.47219        -3.04168
H           9.11710       -10.36990        -2.90047
H           7.03643        -9.19165        -3.67866
H           7.73482        -8.54272        -2.21819
H           0.80472        -6.14860        -1.57838
H          -2.60447        -3.86484        -2.94175
H          -1.21557        -2.74031        -2.92124
H          -2.35422        -2.79675        -1.53792
