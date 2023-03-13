import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from math import log10
from TDE_parameters import *
# from PTF09ge_parameters import *
# from ASASSN14li_parameters import *
import os

##### read data from file

# arrays
with open(os.path.join(folder, 'arrays_info.txt'), 'r') as file:
    file.readline()
    tarr_info = file.readline().strip('\n').split()
    rarr_info = file.readline().strip('\n').split()
    aarr_info = file.readline().strip('\n').split()
    thetaarr_info = file.readline().strip('\n').split()

# tarr 
tst, tend, Nt = float(tarr_info[2]), float(tarr_info[3]), int(tarr_info[4])
if tarr_info[1] == 'linear':
    tarr = np.linspace(tst, tend, Nt)
elif tarr_info[1] == 'logarithmic':
    tarr = np.logspace(log10(tst), log10(tend), Nt)
tarr += (tarr[1] - tarr[0])/2.
# print(tarr)

# rarr
rst, rend, Nr = float(rarr_info[2]), float(rarr_info[3]), int(rarr_info[4])
if rarr_info[1] == 'linear':
    rarr = np.linspace(rst, rend, Nr)
elif rarr_info[1] == 'logarithmic':
    rarr = np.logspace(log10(rst), log10(rend), Nr)
# print(rarr)

# aarr
ast, aend, Na = float(aarr_info[2]), float(aarr_info[3]), int(aarr_info[4])
if aarr_info[1] == 'linear':
    aarr = np.linspace(ast, aend, Na)
elif aarr_info[1] == 'logarithmic':
    aarr = np.logspace(log10(ast), log10(aend), Na)
# print(aarr)

# thetaarr
thetast, thetaend, Ntheta = float(thetaarr_info[2]), float(thetaarr_info[3]), int(thetaarr_info[4])
thetaarr = np.linspace(thetast, thetaend, Ntheta)

# dust temperature
with open(os.path.join(folder, 'dust_temperature.txt'), 'r') as file:
    Tarr_shape = file.readline().strip('\n').split()
    Tarr_plt = np.loadtxt(file)
Tarr_shape = [int(shape) for shape in Tarr_shape]
Tarr_plt = np.reshape(Tarr_plt, Tarr_shape)


##### settings for the plot

# find the maximum temperature to set the ylimit
Tmax = np.max(Tarr_plt)

# the time slice to be plotted for each distance
pltt_st, pltt_end = 0, Nt
pltt = np.arange(pltt_st, pltt_end, int((pltt_end - pltt_st)/3))

# find the range of r to be plotted
# distances at which all dust has sublimated will not be in the range

# pltr_st, pltr_end: the indices of the smallest and biggest distances in plot
# pltr_nonsub: the index of the smallest distance where no dust is sublimated
try: 
    pltr_st = np.argwhere(Tarr_plt[pltt[0], :, -1] == 0).flatten()[-1]
except IndexError:
    pltr_st = 0
try:
    pltr_nonsub = np.argwhere(Tarr_plt[pltt[2], :, 0] != 0).flatten()[0]
except IndexError:
    pltr_nonsub = Nr
pltr_end = Nr

# use different interval for distances with and without sublimation
pltr = np.append(np.arange(pltr_st, pltr_nonsub, 5), np.arange(pltr_nonsub + 3, pltr_end, 10))
plttheta = range(0, Ntheta, 5)

# the number of lines to be plotted
NTplt = max(len(pltr), len(plttheta))

# set the sublimated dust temperature to nan to show a vertical line on plot
# Tarr_plt = Tarr_plt/1e3
Tarr_plt[Tarr_plt==0] = np.nan

# set up the colormap
colors = cm.viridis(np.linspace(0, 1, NTplt))[::-1]

##### plot
fig, axs = plt.subplots(1, 2, figsize=(20, 7))

# r as x_axis
i_plttheta = 2
for i_pltr, rcolor in zip(pltr, colors):
    axs[0].semilogx(aarr, Tarr_plt[pltt[0], i_pltr, i_plttheta, :], ':', color=rcolor, alpha=0.3, linewidth=2)
    vline1 = np.argwhere(np.isnan(Tarr_plt[pltt[0], i_pltr, i_plttheta, :]))
    if (len(vline1) > 0) and (len(vline1) < Na):
        vline1 = max(vline1)+1
        axs[0].vlines(x=aarr[vline1], ymin=Tarr_plt[pltt[0], i_pltr, i_plttheta, vline1], ymax=Tmax,
                   ls=':', color='k', alpha=0.3, lw=2)
    axs[0].semilogx(aarr, Tarr_plt[pltt[1], i_pltr, i_plttheta, :], '--', color=rcolor, alpha=0.6, linewidth=2.5)
    vline2 = np.argwhere(np.isnan(Tarr_plt[pltt[1], i_pltr, i_plttheta, :]))
    if (len(vline2) > 0) and (len(vline2) < Na):
        vline2 = max(vline2)+1
        axs[0].vlines(x=aarr[vline2], ymin=Tarr_plt[pltt[1], i_pltr, i_plttheta, vline2], ymax=Tmax,
                   ls='--', color='k', alpha=0.6, lw=1.5)
    axs[0].semilogx(aarr, Tarr_plt[pltt[2], i_pltr, i_plttheta, :], '-', color=rcolor, alpha=1, linewidth=3)
    vline3 = np.argwhere(np.isnan(Tarr_plt[pltt[2], i_pltr, i_plttheta, :]))
    if (len(vline3) > 0) and (len(vline3) < Na):
        vline3 = max(vline3)+1
        axs[0].vlines(x=aarr[vline3], ymin=Tarr_plt[pltt[2], i_pltr, i_plttheta, vline3], ymax=Tmax, 
                   ls='-', color='k', alpha=0.9, lw=1)
    axs[0].text(x=aarr[-1]*1.03, y=Tarr_plt[0, i_pltr, i_plttheta, -1]*0.99, s='{:<3.2f} pc'.format(rarr[i_pltr]), size=10)
axs[0].set_xlim(0.7*aarr[0], 1.8*aarr[-1])
axs[0].set_xlabel('a [$\\mu m$]', size=15)
# ax.set_ylim(ymax=Tmax/1e3)
axs[0].set_ylabel('T [K]', size=15)
axs[0].grid(True, 'both')
label_handle1 = axs[0].vlines(x=[], ymin=[], ymax=[], ls=':', color='k', alpha=0.3, lw=2,\
     label='t = {:.1e}s'.format(tarr[pltt[0]]))
label_handle2 = axs[0].vlines(x=[], ymin=[], ymax=[], ls='--', color='k', alpha=0.6, lw=1.5,\
     label='t = {:.1e}s'.format(tarr[pltt[1]]))
label_handle3 = axs[0].vlines(x=[], ymin=[], ymax=[], ls='-', color='k', alpha=0.9, lw=1,\
     label='t = {:.1e}s'.format(tarr[pltt[2]]))
axs[0].legend(handles=[label_handle1, label_handle2, label_handle3], fontsize=10, loc='lower left')
axs[0].annotate('{}_T{:.1e}_tdur{:.1e}_nH0{:.2f}_params{:2f}{:2f}'.format(name, T_tde, tdur, nH0, alpha, beta),\
        xy=(0.5, 0.01), xycoords='figure fraction',\
        xytext=(0.5, 0), textcoords= 'figure fraction', \
        horizontalalignment='center', verticalalignment='bottom', size=13)
axs[0].annotate("$n_H$" + " = {:.0f}".format(nH0) + "$cm^{-3}$",\
        xy=(0.15, 0.95), xycoords='axes fraction',\
        xytext=(0.15, 0.95), textcoords= 'axes fraction', \
        horizontalalignment='center', verticalalignment='bottom', size=12)
axs[0].set_title('Dust Temperature', size=15)

# theta as x_axis
i_pltr = 5
for i_plttheta, thetacolor in zip(plttheta, colors):
    axs[1].semilogx(aarr, Tarr_plt[pltt[0], i_pltr, i_plttheta, :], ':', color=thetacolor, alpha=0.3, linewidth=2)
    vline1 = np.argwhere(np.isnan(Tarr_plt[pltt[0], i_pltr, i_plttheta, :]))
    if (len(vline1) > 0) and (len(vline1) < Na):
        vline1 = max(vline1)+1
        axs[1].vlines(x=aarr[vline1], ymin=Tarr_plt[pltt[0], i_pltr, i_plttheta, vline1], ymax=Tmax,
                   ls=':', color='k', alpha=0.3, lw=2)
    axs[1].semilogx(aarr, Tarr_plt[pltt[1], i_pltr, i_plttheta, :], '--', color=thetacolor, alpha=0.6, linewidth=2.5)
    vline2 = np.argwhere(np.isnan(Tarr_plt[pltt[1], i_pltr, i_plttheta, :]))
    if (len(vline2) > 0) and (len(vline2) < Na):
        vline2 = max(vline2)+1
        axs[1].vlines(x=aarr[vline2], ymin=Tarr_plt[pltt[1], i_pltr, i_plttheta, vline2], ymax=Tmax,
                   ls='--', color='k', alpha=0.6, lw=1.5)
    axs[1].semilogx(aarr, Tarr_plt[pltt[2], i_pltr, i_plttheta, :], '-', color=thetacolor, alpha=1, linewidth=3)
    vline3 = np.argwhere(np.isnan(Tarr_plt[pltt[2], i_pltr, i_plttheta, :]))
    if (len(vline3) > 0) and (len(vline3) < Na):
        vline3 = max(vline3)+1
        axs[1].vlines(x=aarr[vline3], ymin=Tarr_plt[pltt[2], i_pltr, i_plttheta, vline3], ymax=Tmax, 
                   ls='-', color='k', alpha=0.9, lw=1)
    axs[1].text(x=aarr[-1]*1.03, y=Tarr_plt[0, i_pltr, i_plttheta, -1]*0.99, s='{:<3.2f} rad'.format(thetaarr[i_plttheta]), size=10)
axs[1].set_xlim(0.7*aarr[0], 1.8*aarr[-1])
axs[1].set_xlabel('a [$\\mu m$]', size=15)
# ax.set_ylim(ymax=Tmax/1e3)
axs[1].set_ylabel('T [K]', size=15)
axs[1].grid(True, 'both')
label_handle1 = axs[1].vlines(x=[], ymin=[], ymax=[], ls=':', color='k', alpha=0.3, lw=2,\
     label='t = {:.1e}s'.format(tarr[pltt[0]]))
label_handle2 = axs[1].vlines(x=[], ymin=[], ymax=[], ls='--', color='k', alpha=0.6, lw=1.5,\
     label='t = {:.1e}s'.format(tarr[pltt[1]]))
label_handle3 = axs[1].vlines(x=[], ymin=[], ymax=[], ls='-', color='k', alpha=0.9, lw=1,\
     label='t = {:.1e}s'.format(tarr[pltt[2]]))
axs[1].legend(handles=[label_handle1, label_handle2, label_handle3], fontsize=10, loc='lower left')
axs[1].annotate('{}_T{:.1e}_tdur{:.1e}_nH0{:.2f}_params{:2f}{:2f}'.format(name, T_tde, tdur, nH0, alpha, beta),\
        xy=(0.5, 0.01), xycoords='figure fraction',\
        xytext=(0.5, 0), textcoords= 'figure fraction', \
        horizontalalignment='center', verticalalignment='bottom', size=13)
axs[1].annotate("$n_H$" + " = {:.0f}".format(nH0) + "$cm^{-3}$",\
        xy=(0.15, 0.95), xycoords='axes fraction',\
        xytext=(0.15, 0.95), textcoords= 'axes fraction', \
        horizontalalignment='center', verticalalignment='bottom', size=12)
axs[1].set_title('Dust Temperature', size=15)

plt.savefig(os.path.join(folder, 'dust_tempearature.jpg'))
plt.show()