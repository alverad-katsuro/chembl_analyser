%nprocshared=4
%mem=512Mb
%chk=chembl3218474.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl3218474 


0  1
C           3.12684         1.13154         0.83278
C           2.56047         0.02859        -0.10045
C           1.02093         0.03822         0.09147
N           2.93388         0.40230        -1.48341
C           2.64397        -0.24701        -2.63269
O           1.89483        -1.19722        -2.73439
O           3.40152         0.27001        -3.62268
C           3.49575        -0.41729        -4.89793
C           2.33223        -0.08610        -5.82243
C           2.44248         0.97081        -6.72880
C           1.42539         1.22083        -7.65033
C           0.28743         0.41874        -7.66609
C           0.16465        -0.63081        -6.75729
C           1.18507        -0.88540        -5.84306
C           3.13981        -1.33177         0.32569
O           3.96857        -1.44780         1.20858
N           2.64490        -2.42630        -0.29420
C           2.95418        -3.81382         0.06947
C           1.89162        -4.37062         1.00293
N           1.01014        -4.83253         1.77731
H           2.72848         2.11375         0.56586
H           4.21658         1.18215         0.75961
H           2.86695         0.93565         1.87639
H           0.60546         1.01947        -0.15122
H           0.76155        -0.19220         1.12807
H           0.52804        -0.69907        -0.54479
H           3.60117         1.15336        -1.55767
H           3.56404        -1.49526        -4.71984
H           4.43094        -0.09940        -5.36682
H           3.33495         1.59211        -6.73350
H           1.52764         2.03396        -8.36369
H          -0.50128         0.60707        -8.39034
H          -0.72171        -1.25965        -6.77174
H           1.08955        -1.72381        -5.15785
H           1.95440        -2.29036        -1.01786
H           3.92916        -3.87035         0.56115
H           2.97889        -4.42138        -0.83842

