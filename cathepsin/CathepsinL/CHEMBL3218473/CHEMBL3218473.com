%nprocshared=4
%mem=512Mb
%chk=chembl3218473.chk 
# opt freq=noraman m062x/6-311G(d,p) scf=maxcycle=1000 maxdisk=100Gb 

chembl3218473 


0  1
N           5.71194        -2.32402         4.67413
C           6.26204        -2.19046         5.66098
C           6.99474        -1.98739         6.92023
N           7.81246        -0.77212         6.88357
C           7.35509         0.47096         7.14818
O           6.20297         0.64132         7.50017
C           8.35062         1.62764         6.95959
N           8.73813         1.78640         5.53594
C           9.41226         0.92460         4.73784
O           9.95737        -0.08356         5.14314
C           9.60410         1.43136         3.34433
C          11.01512         1.57029         2.85364
C          10.12443         0.50712         2.27601
C           9.60505         0.57011         0.87766
C          10.38682         1.09572        -0.15490
C           9.88980         1.15115        -1.45557
C           8.61239         0.67161        -1.73555
C           7.83276         0.13576        -0.71211
C           8.32748         0.08611         0.58992
C           9.58666         1.40402         7.88692
C          10.57257         2.59510         7.86420
C           9.86724         3.89685         8.29306
C           8.63853         4.16605         7.40236
C           7.67387         2.95706         7.41944
H           7.64940        -2.84648         7.08486
H           6.27845        -1.92610         7.74418
H           8.75645        -0.86552         6.53658
H           8.38132         2.62108         5.09747
H           8.89820         2.19193         3.01870
H          11.82116         1.24741         3.50620
H          11.24772         2.44196         2.24945
H          10.34197        -0.51714         2.56858
H          11.38902         1.46538         0.04853
H          10.50130         1.56661        -2.25316
H           8.22629         0.71214        -2.75072
H           6.83781        -0.24490        -0.93097
H           7.71287        -0.33581         1.38248
H          10.12625         0.49680         7.60889
H           9.24244         1.25191         8.91333
H          10.99387         2.71403         6.86298
H          11.40479         2.39254         8.54240
H          10.56377         4.73542         8.22605
H           9.55220         3.81584         9.33697
H           8.96163         4.37783         6.38041
H           8.11506         5.05375         7.76508
H           6.81673         3.17783         6.77759
H           7.28681         2.84089         8.43585
