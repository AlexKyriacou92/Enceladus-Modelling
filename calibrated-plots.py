#Mercator Plots
#Global Plots

from enceladus_modelling import *

geyser_list = Geyser.load_geyser_list(f_tiger_stripe)
all_geysers = get_array('geyser-map.npy')

geyser_0 = Geyser.get_geyser(10, geyser_list)

heatmap_0 = geyser_0.get_heatmap(0, geyser_list, all_geysers)
print(heatmap_0)

#Generate Global Plots!

#geyser_0.plot_heatmap(0, geyser_list, all_geysers, 'lin')

#Create Calibrated Mass Depsoition
#Sum of all geysers -> Mrate

Mrate = 25 #kg/s
N_latitude = 180
N_longitude = 360

Nradius_spectrum = np.load(path_to_directory + 'N_particle_spectrum.npy', 'r')
radii = geyser_0.radius_list
N_collsions = Nradius_spectrum[:,2]
N_launched = Nradius_spectrum[:,1]
P_collsion = N_collsions/N_launched #probability of landing inside radius_range

N_collsions_norm = N_collsions/sum(N_launched)
N_launched_norm = N_launched/sum(N_launched)

#mrate_geyser_heatmap = open_memmap(path_to_directory + 'mass_rate_surface.npy', dtype='float32', mode='r+', shape = (N_geysers, N_radii, N_latitude, N_longitude))
'''
from scipy.interpolate import interp1d
r_list[0] = 0.6e-6
dr = 0.1e-6 #m
r_max = max(r_arr*1e-6)
r_min = min(r_arr*1e-6)
Nr = int((r_max - r_min)/dr)

r_space = np.linspace(r_min, r_max, Nr)
P = interp1d(radii, dP_dr)
'''
normalize_geysers = np.zeros((N_radii, N_latitude, N_longitude))

#Normalization
for geyser_i in range(N_geysers):
    for radius_j in range(N_radii):
        normalize_geysers[radius_j] += all_geysers[geyser_i][radius_j]

def sum2d(arr):
    norm = 0
    Narr = len(arr)
    for i in range(Narr):
        M = len(arr[i])
        for j in range(M):
            norm += arr[i][j]
    return norm

print('')

Normalization = 0
R_Norm = 0
U_Norm = 0
M_rate = 25 #Mass Eruption Rate kg/s
for j in range(N_radii):
    print(j, 'r = ', radii[j], 'm, P = ', P_collsion[j]*100, '%, norm =', sum_2d(normalize_geysers[j]) / N_geysers)
    Normalization += sum_2d(normalize_geysers[j])
    U_Norm += (4/3)*pi*rho_ice * radii[j]**3 * N_launched[j] / sum(N_launched)

N_particles = M_rate/U_Norm #Total number of particles launched per second
print('Particle launch rate', N_particles, 'particles/s')
print('Number of landing particles', N_particles * sum(N_collsions)/sum(N_launched), 'particles/s')
print('Number of launched particles', N_particles, 'particles/s')

for j in range(N_radii):
    mass_j = (4/3) * pi * rho_ice * radii[j]**3
    m_rate_launched_j = N_particles * N_launched[j]/sum(N_launched) * mass_j
    m_rate_landed_j = N_particles * N_collsions[j]/sum(N_launched) * mass_j
    print('r = ', radii[j]*1e6 ,'um, m =', mass_j, 'kg/s, m_rate_launched = ', m_rate_launched_j, 'kg/s, m_rate_landed = ', m_rate_landed_j, 'kg/s')

#Plot Mrate
mrate_geyser_heatmap = open_memmap(path_to_directory + 'mass_rate_surface.npy', dtype='float32', mode='w+', shape = (N_geysers, N_radii, N_latitude, N_longitude))
Zrate_geyser_heatmap = open_memmap(path_to_directory + 'Zrate_surface.npy', dtype='float32', mode='w+', shape = (N_geysers, N_radii, N_latitude, N_longitude))
Nrate_geyser_heatmap = open_memmap(path_to_directory + 'Nrate_surface.npy', dtype='float32', mode='w+', shape = (N_geysers, N_radii, N_latitude, N_longitude))

mrate_sum_all = open_memmap(path_to_directory + 'mass_rate_surface_sum.npy', dtype='float32', mode='w+', shape = (N_latitude, N_longitude)) # to map particle flux rate
Nrate_sum_all = open_memmap(path_to_directory + 'Nrate_surface_sum.npy', dtype='float32', mode='w+', shape = (N_latitude, N_longitude)) # to map depsoition
Zrate_sum_all = open_memmap(path_to_directory + 'Zrate_rate_surface_sum.npy', dtype='float32', mode='w+', shape = (N_latitude, N_longitude)) # to map total mass rate
Area_Array = np.load(path_to_directory + 'area_array.npy', 'r')

year = 365.25 * 3600 * 24




for geyser_i in range(N_geysers):
    for radius_j in range(N_radii):
        #mrate_sum_all += all_geysers[geyser_i][radius_j] * (((4 / 3) * pi * rho_ice * P(r_arr[j] * 1e-6) * (r_arr[j] ** 3)) * dr) / Area_Array
        #mass_j = (4 / 3) * pi * rho_ice * radii[j] ** 3
        #vol_j = ((4 / 3) * pi * radii[j] ** 3) * (N_particles * N_collsions[j] / sum(N_launched))
        #vol_j = particle_volume(radii[j]) * N_particles * (N_collsions[j]/sum(N_launched))
        N_particles_j = N_particles * N_collsions[radius_j] / sum(N_launched)
        vol_j = N_particles_j * particle_volume(radii[radius_j])
        mass_j = vol_j * rho_ice * N_particles_j
        #print(radii[j], radii[j]**3)

        #Mass Rate
        m_rate_landed_j = mass_j
        mrate_geyser_heatmap[geyser_i][radius_j] = m_rate_landed_j * all_geysers[geyser_i][radius_j] / Area_Array
        mrate_sum_all += m_rate_landed_j * all_geysers[geyser_i][radius_j] / Area_Array

        #Z Rate [mm/yr]
        Zrate_geyser_heatmap[geyser_i][radius_j] = vol_j * all_geysers[geyser_i][radius_j] / Area_Array * (year) #Note, is this m / yr or mm/yr -> check you calibration!!
        Zrate_sum_all +=  vol_j * all_geysers[geyser_i][radius_j] / Area_Array * (year) #TODO Check Deposition Units

        #N rate
        Nrate_geyser_heatmap[geyser_i][radius_j] = N_particles_j * all_geysers[geyser_i][radius_j] / Area_Array
        Nrate_sum_all += N_particles_j * all_geysers[geyser_i][radius_j] / Area_Array

#Plot Mass Deposition Rate
fig, ax, cbar = heatmap2d(mrate_sum_all, mode='log', alpha_value= 1.0)
ax.set_title('Mass Flux Rate [$\mathrm{kg s^{-1}}$]')
cbar.set_label('$\dot{M}$ [$\mathrm{kg m^{-2}s^{-1}}$')

show_plot(fig, True)

#Map surface deposition rate
img_1 = pl.imread(path_to_directory + 'Enceladus-Images/enceladus-mercator-compressed.jpg')
fig, ax, cbar = heatmap2d(Zrate_sum_all, mode='log', alpha_value= 0.5, vmin0=1e-6, vmax0=1e0, my_cmap='Blues')
ax.set_title('Deposition Rate [$\mathrm{mm yr^{-1}}$]')
cbar.set_label('$\dot{Z}$ [$\mathrm{mm yr^{-1}}$')
pl.savefig(path_to_directory + 'plots/mercator-zrate-map.png',bbox_inches='tight')
#pl.close(fig)
show_plot(fig)

fig0 = pl.figure(figsize=(15,15), dpi = 100)
ax, cbar = make_polar_plot(Zrate_sum_all, lat_max= -60, fig=fig0, Pixel_Radius=440)
pl.show()

#Collision FLux
#Map N flux rate
fig, ax, cbar = heatmap2d(Nrate_sum_all, mode='log', alpha_value= 0.5)
ax.set_title('Collision FLux [$\mathrm{m^{-2}s^{-1}}$]')

cbar.set_label('$\dot{N}$ [$\mathrm{particles m^{-2}s^{-1}}$')

show_plot(fig, True)

#print(Area_Array)
#sum2d(Area_Array)
A_enceladus = 4 * pi * (252*1e3)**2
print(A_enceladus, sum2d(Area_Array), A_enceladus/sum2d(Area_Array))

#Different Time Scenarios

T_active = [1e3, 1e4, 1e5, 1e6]
T_labels = ['1 kyr', '10 kyr', '100 kyr', '1 Myr']
for i in range(len(T_active)):
    # Map surface deposition rate
    z_deposition = T_active[i]*Zrate_sum_all/1e3

    fig, ax, cbar = heatmap2d(z_deposition, mode='log', alpha_value=1.0)
    ax.set_title('Surface Deposition [m], T = ' + T_labels[i])
    cbar.set_label('$\Delta Z $ [m]')
    pl.savefig('T_active=' + str(T_active) + 'yr.png')
    #pl.show()
    pl.close(fig)
    make_polar_plot(z_deposition, -60, lim_coords=[[-285, 285], [-285, 285]], fig=None, index=111, img=pl.imread(thirty_deg_img))
