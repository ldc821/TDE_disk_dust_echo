import numpy as np
import os
from math import pi, log10, cos, sin, acos
import constCGS

##### import parameters

from TDE_parameters import *
from TDE_observation_setup import *

##### read data from file

# arrays
with open(os.path.join(folder, 'arrays_info.txt'), 'r') as file:
    file.readline()
    tarr_info = file.readline().strip('\n').split()
    rarr_info = file.readline().strip('\n').split()
    aarr_info = file.readline().strip('\n').split()
    thetaarr_info = file.readline().strip('\n').split()

# emissvity 
with open(os.path.join(folder, 'emissivity.txt'), 'r') as file:
    jdnuarr_shape = file.readline().strip('\n').split()
    jdnuarr = np.loadtxt(file)
jdnuarr_shape = [int(shape) for shape in jdnuarr_shape]
jdnuarr = np.reshape(jdnuarr, jdnuarr_shape)

##### discretization

# tarr 
tst, tend, Nt = float(tarr_info[2]), float(tarr_info[3]), int(tarr_info[4])
if tarr_info[1] == 'linear':
    tarr = np.linspace(tst, tend, Nt)
elif tarr_info[1] == 'logarithmic':
    tarr = np.logspace(log10(tst), log10(tend), Nt)

# rarr
rst, rend, Nr = float(rarr_info[2]), float(rarr_info[3]), int(rarr_info[4])
# if rarr_info[1] == 'linear':
#     rarr = np.linspace(rst, rend, Nr)
# elif rarr_info[1] == 'logarithmic':
#     rarr = np.logspace(log10(rst), log10(rend), Nr)
rarr = np.linspace(rst, rend, Nr)

# aarr
ast, aend, Na = float(aarr_info[2]), float(aarr_info[3]), int(aarr_info[4])
if aarr_info[1] == 'linear':
    aarr = np.linspace(ast, aend, Na)
elif aarr_info[1] == 'logarithmic':
    aarr = np.logspace(log10(ast), log10(aend), Na)

# thetaarr: polar angle
thetast, thetaend, Ntheta = float(thetaarr_info[2]), float(thetaarr_info[3]), int(thetaarr_info[4])
thetaarr = np.linspace(thetast, thetaend, Ntheta)

# phiarr: azimuthal angle
phiarr = np.linspace(0, 2*pi, Nphi) # rad

# observation frequency
lambobsarr = np.linspace(lambobsmin, lambobsmax, Nnuobs) # um
nuobsarr = constCGS.C_LIGHT/(lambobsarr/constCGS.cm2um) # in Hz

# observation time
tobsarr = np.logspace(log10(tobsmin), log10(tobsmax), Ntobs)
# tobsarr = np.linspace(tobsmin, tobsmax, Ntobs)

# specific luminosity
Ldnuarr = np.zeros((Nnuobs, Ntobs), dtype=float)

##### functions

# linear interpolation of the emissivity at a given distance
def jdnu_intp(t, i_r, i_theta, i_nuobs):
    i_floor = np.argmin(np.abs(t - tarr))
    if t > tarr[-1] + (tarr[-1] - tarr[-2]):
        return 0
    if tarr[i_floor] == t:
        return jdnuarr[i_floor, i_r, i_theta, i_nuobs]
    if (tarr[i_floor] > t and i_floor > 0) or i_floor == Nt-1:
        i_floor -= 1
    slope = (jdnuarr[i_floor+1, i_r, i_theta, i_nuobs] - jdnuarr[i_floor, i_r, i_theta, i_nuobs])/(tarr[i_floor+1] - tarr[i_floor])
    return max(jdnuarr[i_floor, i_r, i_theta, i_nuobs] + slope * (t - tarr[i_floor]), 0)

# calculate mu from tobs and t, Eq. 27
def mu(theta, phi, thetaobs):
    return cos(theta)*cos(phi)*sin(thetaobs) + sin(theta)*cos(thetaobs)

# solve tobs at a given r from mu
def mu2t(mu, rcm, tobs):
    return tobs - rcm/constCGS.C_LIGHT*(1-mu)

# specifc luminosity, Eq.28
def lum_dnu(tobs, i_nuobs):
    integral = 0
    for i_r in range(Nr-1):
        rcm = rarr[i_r]*constCGS.pc2cm
        drcm = (rarr[i_r+1] - rarr[i_r])*constCGS.pc2cm
        # light echo has already passed or light has not arrived
        # if rcm < constCGS.C_LIGHT*(tobs-tdur)/2 or rcm > constCGS.C_LIGHT*(tobs-tmin):
        #     continue 
        ang_integral = 0
        for i_theta in range(Ntheta-1):
            theta = thetaarr[i_theta]
            dtheta = thetaarr[i_theta+1] - theta
            for i_phi in range(Nphi-1):
                phi = phiarr[i_phi]
                dphi = phiarr[i_phi+1] - phi
                mu_pos = mu(theta, phi, thetaobs)
                mu_neg = mu(-theta, phi, thetaobs)
                if ((t_pos := mu2t(mu_pos, rcm, tobs)) > 0 and t_pos < tdur):
                    ang_integral += cos(theta) * dtheta * dphi * jdnu_intp(t_pos, i_r, i_theta, i_nuobs)
                if ((t_neg := mu2t(mu_neg, rcm, tobs)) > 0 and t_neg < tdur):
                    ang_integral += cos(theta) * dtheta * dphi * jdnu_intp(t_neg, i_r, i_theta, i_nuobs)
        integral += rcm**2 * ang_integral * drcm
    return 4 * pi * integral
            
##### calculate specific luminosity

print('Calculating...')
progress = 0

for i_tobs in range(Ntobs):
    tobs = tobsarr[i_tobs]
    for i_nuobs in range(Nnuobs):
        Ldnuarr[i_nuobs, i_tobs] = lum_dnu(tobs, i_nuobs)
    if i_tobs/Ntobs > progress:
        print('{:.2%}'.format(i_tobs/Ntobs))
        progress += 0.1
print('{:.2%}'.format(1))

##### save data

print('\nSaving data...')

Ldnuarr_shape = '{}\t{}\n'.format(Nnuobs, Ntobs)
with open(os.path.join(folder, 'specific_luminosity.txt'), 'w') as file:
    file.write(Ldnuarr_shape)
    np.savetxt(file, Ldnuarr)

print('All done!')