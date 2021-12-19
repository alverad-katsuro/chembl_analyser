%nprocshared=4
%mem=512Mb
%chk=chembl1836510.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl1836510 


0  1
N          -3.46152        -0.25288         1.81742
C          -2.37341         0.08282         1.71416
C          -1.00144         0.49408         1.57493
C          -0.25760         0.08700         0.46703
C           1.07062         0.48655         0.33105
C           1.66656         1.29510         1.30030
C           3.11035         1.73615         1.13009
S           3.17642         3.29523         0.18709
C           4.88916         3.74589         0.20734
N           5.85151         2.78854         0.22606
C           7.13348         3.22890         0.24846
S           8.42775         2.01415         0.25753
C           9.80013         2.77724         1.18711
C          10.96246         1.80745         1.31157
C          11.94322         1.75336         0.31888
C          13.02100         0.87753         0.44297
C          13.12294         0.05041         1.55958
C          12.14331         0.09558         2.54940
C          11.06535         0.97100         2.42525
N           7.46374         4.54608         0.24819
C           6.70768         5.36775         0.23014
O           7.12243         6.50800         0.23580
N           5.37926         5.09484         0.20321
C           0.91567         1.70064         2.40818
F           1.46290         2.48547         3.35441
C          -0.41234         1.30267         2.54638
H          -0.71203        -0.54182        -0.29656
H           1.64078         0.16831        -0.53892
H           3.68349         0.97567         0.59315
H           3.58410         1.89620         2.10234
H           9.43281         3.05990         2.17726
H          10.11726         3.68526         0.66730
H          11.87196         2.39695        -0.55488
H          13.78437         0.84172        -0.33056
H          13.96605        -0.62928         1.65868
H          12.22171        -0.55043         3.42056
H          10.30607         1.00113         3.20325
H           4.70648         5.84608         0.19042
H          -0.98635         1.62680         3.41200
