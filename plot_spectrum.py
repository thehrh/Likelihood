# Plot arrival time spectrum
# Set mass, dist, number of events and detector type

# imported from c-program
import spectrum

import os
import pickle
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import dblquad, quad
from numpy import transpose
from matplotlib.colors import LogNorm


spectrum_plot = spectrum.doubleArray(599000)
trigger_eff = spectrum.doubleArray(601)
for i in range(0,601):
    trigger_eff[i] = 1.0

#spectrum.generateDist(1.0,1.0, 160, spectrum_plot, trigger_eff, False)

spectrum.createSpectrum(spectrum_plot,1.0, 1.0,  160,True, True,pow(10,-5))

myArray = [[spectrum_plot[t*(599) +e] for t in range(0, 1000, 1)] for e in range(1, 590)]

X = np.arange(0, 10., 0.01)
Y = np.arange(0.1, 59., 0.1)
X, Y = np.meshgrid(X, Y)

Z = myArray

# Surface Plot of the arrival distribution - log scale
fig = plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.87)
#surf = ax.contourf(X,Y,Z, 8, cmap=plt.cm.jet)
surf = ax.contourf(X,Y,Z,levels=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],cmap=plt.cm.jet,norm = LogNorm())
#surf = ax.contourf(X,Y,Z,cmap=plt.cm.jet,norm = LogNorm())
ax.set_xlabel('time [s]', fontsize=19)
ax.set_ylabel('energy [MeV]', fontsize=19)
#ax.set_title('m = '+str(M)+' eV - ' + str(events) + ' events - D = '+str(D)+ ' Mpc \n'+str(det_type), fontsize=19)
ax.xaxis.set_tick_params(labelsize=19, width=2)
ax.yaxis.set_tick_params(labelsize=19, width=2)
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
# defining custom minor tick locations:
ax.xaxis.set_minor_locator(plt.FixedLocator([50,500,2000]))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both',reset=False,which='both',length=8,width=2)
cbar = fig.colorbar(surf, shrink=1, aspect=20, fraction=.12,pad=.02)
cbar.set_label('# of events',size=19)
# access to cbar tick labels:
cbar.ax.tick_params(labelsize=19)
#plt.xlim(0.0, 10.0)
#plt.ylim(1, 39)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(top=0.87)
surf = ax.contourf(X,Y,Z, 8, cmap=plt.cm.jet)
#surf = ax.contourf(X,Y,Z,levels=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],cmap=plt.cm.jet,norm = LogNorm())
#surf = ax.contourf(X,Y,Z,cmap=plt.cm.jet,norm = LogNorm())
ax.set_xlabel('time [s]', fontsize=19)
ax.set_ylabel('energy [MeV]', fontsize=19)
#ax.set_title('m = '+str(M)+' eV - ' + str(events) + ' events - D = '+str(D)+ ' Mpc \n'+str(det_type), fontsize=19)
ax.xaxis.set_tick_params(labelsize=19, width=2)
ax.yaxis.set_tick_params(labelsize=19, width=2)
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
# defining custom minor tick locations:
ax.xaxis.set_minor_locator(plt.FixedLocator([50,500,2000]))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='both',reset=False,which='both',length=8,width=2)
cbar = fig.colorbar(surf, shrink=1, aspect=20, fraction=.12,pad=.02)
cbar.set_label('# of events',size=19)
# access to cbar tick labels:
cbar.ax.tick_params(labelsize=19)
#plt.xlim(0.0, 10.0)
#plt.ylim(1, 39)
plt.show()

