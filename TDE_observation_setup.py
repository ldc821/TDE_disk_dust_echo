from math import pi
import constCGS

##### observation setup #####

Nphi = 20

lambobsmin, lambobsmax = 2, 4  # [um] observation wavelength
Nnuobs = 10

Nmu = 100                   

Ntobs = 100
tobsmin = 100
tobsmax = 200*constCGS.pc2cm/constCGS.C_LIGHT

thetaobs = pi/2 # observation angle