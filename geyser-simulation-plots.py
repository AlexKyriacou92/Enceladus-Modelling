from enceladus_modelling import *

from numpy import sqrt
import sys
import os
import scipy
sys.path.append("./pymap3d/src/")
from matplotlib import pyplot as pl
import pymap3d as pm
enceladus_ellipsoid = pm.Ellipsoid(model = "enceladus")
#Velocity Probability Distribution
cmap = matplotlib.cm.get_cmap('gnuplot2')
day = 3600*24.
#Rotation Period of Enceladus
T_Enceladus = 1.37*day
#Rotation Rate of Enceladus
Omega_Enceladus = 2*pi / T_Enceladus

#LANDING POSITION SIMULATION FOR V <= 20 m/s

x_final = np.load('xposition_2.npy') #for v < 20 m/s
n_final = np.load('nposition_2.npy')
n_final_sum = np.load('nposition_sum_2.npy')

deg_min = 0
deg_max = 2 #

N_v2 = len(n_final)
N_a2 = len(n_final[0])
v_max = float(N_v2)
vspace = np.linspace(0, v_max, N_v2)
theta_dist = np.linspace(deg_min, deg_max, N_a2) #TODO: Check negative values
#vspace = np.linspace(0, 700, N_velocities+1)

topo = np.genfromtxt(path_to_directory + 'damascus-topography.txt')
X_topo, Y_topo = topo[:,0], topo[:,1]
X_center = 3235

n_list = []
x_list = []

x_max = 2000.
x_min = 0
dx = 10.
xbins = int((x_max - x_min)/dx)

n_flux = np.zeros(xbins)
z_hist = np.zeros(xbins)
#TODO: ax[0,1]: topography, a[1,0]: z_rate, a[1,1]: topography
fig, axs = pl.subplots(2,2, figsize=(18,8))
for r in range(N_radii):
    n_list_r = []
    for i in range(N_v2):
        for j in range(N_a2):
            n_list_r.append(n_final[r,i,j])
            if r == 0:
                x_list.append(x_final[i,j])
                n_list.append(n_final_sum[i,j])
    n_hist_r, x_bins_r = np.histogram(x_list, bins=xbins, range=(x_min, x_max), weights=n_list_r)
    n_flux_r = n_hist_r / (dx * 2 * pi * x_bins_r[1:])
    z_flux_r = n_flux_r * particle_volume(ice_grain_radii[r])*year

    n_flux += n_flux_r
    z_hist += z_flux_r

    label_r = '$r_{grain}$ = ' + str(round(ice_grain_radii[r]*1e5,2)) + ' $\mathrm{\mu m}$'
    color_r = cmap(ice_grain_radii[r]/11e-6)
    axs[0,0].plot(x_bins_r[1:], n_flux_r, label=label_r,c=color_r)
    #axs[0,1].plot(x_bins_r[1:], z_flux_r, label=label_r,c=color_r)

n_hist, x_bins = np.histogram(x_list, bins=xbins, range=(x_min, x_max), weights=n_list)

axs[0,0].set_title('Histogram, landing particles')
#axs[0,0].plot(x_bins[1:], n_hist, label='Sum over Radii',c='k')
axs[0,0].plot(x_bins[1:], n_flux, label='Sum over Radii',c='k')
axs[0, 1].plot(x_bins_r[1:], z_hist, label='Deposition Rate', c='k')
axs[0,1].set_ylabel('Depososition Rate [mm/yr]')

axs_2 = axs[0,1].twinx()
axs_2.plot(x_bins_r[1:], (z_hist*100e3) /1e3, c='b',label='Accumulated Ice')
axs_2.set_ylabel('Accumulated Ice [100 kyr]')
axs_2.set_yscale('log')
axs[0,0].set_yscale('log')
axs[0,0].set_ylabel('$\dfrac{dN}{dAdt}$ [$\mathrm{m^{-2}s^{-1}}$]')
#axs[0,1].set_ylabel('Ice Deposition [mm/yr]')

axs[1,0].plot(X_topo - X_center, Y_topo)
axs[1,0].set_xlabel('X [m]')
axs[1,0].set_ylabel('Y [m]')
axs[1,1].plot(X_topo - X_center, Y_topo)
axs[1,1].set_xlabel('X [m]')
axs[1,1].set_ylabel('Y [m]')

#ADD TEMPERATURE SINTERING
#Interplote temperature with distance data
from scipy import interpolate
sintering_times = np.genfromtxt(path_to_directory + 'sintering-time.txt')
T_sinter = [80, 100, 120, 140, 160, 180, 200, 220]
yr = 365.25*24*3600
t_sinter = sintering_times[:,0]/yr
r_sinter = sintering_times[:,1]

sinter = interpolate.interp1d(T_sinter, np.log10(t_sinter))
Temp_space = np.linspace(80, 200, 100)
Tsinter_space = sinter(Temp_space)

def get_sintering_time(T):
    sinter_log10 = sinter(T)
    t = 10**sinter_log10
    return t

surface_temp = np.genfromtxt(path_to_directory + 'T_Fracture_225.txt')


print(surface_temp)
x_temp = surface_temp[:,0] #Horizontal distance
T_temp = surface_temp[:,1] #Surface Temperature

f = interpolate.interp1d(T_temp, x_temp)
Temp_space = np.linspace(min(T_temp), max(T_temp), 100)
x_interp = f(Temp_space)

X_plot = X_topo - X_center
g = interpolate.interp1d(X_plot, Y_topo) #interpolate topographic X and Y coordinates (in metres)
x_space = np.linspace(min(X_plot), max(X_plot), 1000)
#print(x_space, min(X_topo))
y_space = g(x_space)
#pl.plot(x_space, y_space)

y_space = g(x_space)
#x_local = x_ref + x_space/1e3
x_local = x_space
Temp_list = [80, 90, 100, 120, 140, 160]

for i in range(len(Temp_list)):
    Ti = Temp_list[i]
    xi = f(Ti)
    print('xi', xi, Ti)

    yi = g(xi)
    print('local coords')
    print(xi, yi)
    tsinter = get_sintering_time(Ti)
    sinter_label = ', $ t_{sinter} = 10^{' + str(round(np.log10(tsinter), 1)) + '}$ yrs,  $10^{' + str(round(np.log10(tsinter * yr), 1)) + '}$ s'
    axs[1,1].scatter(xi, yi, label='T = ' + str(Ti) + ' K, ' + sinter_label)
    axs[1,0].scatter(xi, yi)

for m in range(2):
    axs[0, m].set_yscale('log')
    for n in range(2):
        axs[m,n].set_xlabel('X [m]')
        axs[m,n].legend()
        axs[m,n].grid()
        axs[m,n].set_xlim(x_min-10, x_max)

pl.show()