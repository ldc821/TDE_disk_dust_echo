##### model parameters #####
from math import pi 

# tidal disruption events
T_tde = 3e4             # [K], temperature of TDE
lumr = 100              # [au], luminous radius of the blackhole
tdur = 1e6              # [s], duration of the burst

# dust distribution
nH0 = 1                # [cm^{-3}], H number density
alpha = 0.5         # the exponent of the density profile
beta = 0
# beta = 1/(pi/4)
# beta = 1/(pi/6)
n02nH = 1.45e-15        # n0 over nH(r)
lamb0 = 2               # [um], critical wavelength

##### simulation setup #####

amin, amax = 0.01, 0.3      # [um] min/max grain radius
Na = 30

rmin, rmax = 0.4, 100.      # [pc] radial layers
Nr = 100

tmin, tmax = 0, tdur        # [s]  time since start of TDE events
Nt = 50

hnumin, hnumax = 0.1, 50    # [eV] source frequency
Nnu = 100

thetamin, thetamax = 0, pi/2
Ntheta = 20

tol = 0.1                   # fractional tolerance for sublimation radius

##### folder to store the outputs
name = "data"
folder = '{}_T{:.1e}_tdur{:.1e}_nH0{:.2f}_params{:.1f}{:.1f}'.format(name, T_tde, tdur, nH0, alpha, beta)