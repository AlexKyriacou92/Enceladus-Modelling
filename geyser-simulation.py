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

day = 3600*24.
#Rotation Period of Enceladus
T_Enceladus = 1.37*day
#Rotation Rate of Enceladus
Omega_Enceladus = 2*pi / T_Enceladus

def p(v,r,r_c = R_critical,v_gas = V_gas):
    return (1 + r/r_c) * (r/r_c) * (v/(v_gas**2)) * (1 - (v/v_gas))**((r/r_c) - 1)

def get_magnitude(vec):
    return sqrt(sum(vec**2))
#3D COORDINATE SYSTEM

#STEP 1 SAMPLE FROM VELOCITY AND ANGLE DISTRIBUTION

#STEP 2: Simulate plume flights

#STEP 3: Save landing positions

#STEP 4: Scale by Particle rate

#STEP 5: FIND DELTA Z /DT

#STEP 6: Apply Timing Model

#STEP0:
#SELECT GEYSER YOU WANT TO SIMULATE

geyser_list = Geyser.load_geyser_list(f_tiger_stripe)
geyser_0 = Geyser.get_geyser(379, geyser_list)

N_velocities = 700 #from 0 m/s to 800 m/s
nbins = N_velocities
v_space_linear = np.linspace(0, V_gas, N_velocities*10)

N_angles = 60 #from -7.5 deg to +7.5 deg
#geyser_sample = open_memmap('geyser_velocity_dist.npy', mode='w+', dtype='int', shape=(N_radii, N_velocities))
#N_sample = int(N_particles/f_reduction)
cmap = matplotlib.cm.get_cmap('gnuplot2')

Nrate_geyser_heatmap = get_array('Nrate_surface.npy', path_to_directory, 'r+')
Nradius_spectrum = np.load(path_to_directory + 'N_particle_spectrum.npy', 'r')
#radii = geyser_0.radius_list
N_collisions = Nradius_spectrum[:,2]
N_launched = Nradius_spectrum[:,1]
f_ratio = 1e7 #COntrols number of particles simulated -> reduce for speed, increase for smoother distributions

zenith = geyser_0.zenith_d
z_width = geyser_0.opening_angle_d/2.
dtheta_dist = np.random.uniform(zenith-z_width, zenith+z_width, N_angles)

'''
for r in range(N_radii):
    #GET NUMBER OF PARTICLES
    heatmap_0 = geyser_0.get_heatmap(r, geyser_list, Nrate_geyser_heatmap)
    # GET LOCAL COORDINATE
    latitudes = np.linspace(-90, 90, 180)
    longitudes = np.linspace(-180, 180, 360) #LONGITUDE DEfINED FROM -180 to +180

    N_flux_local = 0  # Number of Ice Grains falling per unit area per unit time %m^2

    N_particles_tot = 68293003296349.93

    N_particles = N_launched[r]/sum(N_launched) * N_particles_tot
    print('ratio',N_launched[r]/N_collisions[r])

    N_sample = int(float(N_particles)/f_ratio)
    print(N_particles, N_sample)

    np.random.seed(1)
    # x_rand = v_escape*np.random.sample(N_sample)
    x_rand = V_gas * np.random.sample(10 * N_sample)
    #====================================================================================================
    r_grain = ice_grain_radii[r]
    velocity_probabilty = p(v_space_linear, r_grain)
    p_max = max(velocity_probabilty)

    np.random.seed(2)
    y_rand = p_max * np.random.sample(10*N_sample)
    print(r, r_grain)
    v_dist = []
    ii = 1
    #for i in range(N_sample):
    while ii < N_sample:
        i = np.random.randint(0, 10*N_sample - 1)
        y_trial =  p(x_rand[i], r_grain)
        if y_trial >= y_rand[i]:
            v_dist.append(x_rand[i])
            ii+=1
            #print(100*float(ii)/float(N_sample))

    N_survivors = len(v_dist)
    print('N Survivors', N_survivors, 'fraction = ', 100*float(N_survivors)/float(N_sample), '%')

    hist_r, v_bins = np.histogram(v_dist, bins=N_velocities, range=(0,700.))
    #print(len(hist_r), hist_r)
    hist_r *= int(f_ratio)
    #print(hist_r)
    geyser_sample[r][:] = hist_r
    #pl.hist(v_dist, bins= 100)
    #pl.show()
    '''

geyser_sample2 = np.load('geyser_velocity_dist.npy','r+')
#print(geyser_sample[6])
v_bins = np.linspace(0,700,N_velocities)
#pl.plot(v_bins, geyser_sample2[6])
#pl.show()
#STEP 1


#norm vector to surface:
def get_vec(latitude, longitude, h): #Degrees
    x,y,z = pm.geodetic2ecef(latitude, longitude, h, ell= enceladus_ellipsoid, deg=True)
    return np.array([x, y ,z])

def get_norm(latitude, longitude, h):
    vec_xyz = get_vec(latitude, longitude, h)
    norm = sqrt(sum(vec_xyz**2))
    return vec_xyz/norm

def get_antinorm(latitude, longitude, h):
    return -get_norm(latitude, longitude, h)

def coriolis_vector(v_vec, rot_vec = Omega_Enceladus*np.array([0,0,1])):
    return -2*np.cross(v_vec*rot_vec)


#CREATE 3D Acceleration Vector:
def acceleration_simple(r_vec):
    R_saturn = R_orbit * np.array([1, 0, 0])
    r_saturn = R_saturn + r_vec
    g_Enceladus3 = gravitational_vector(r_vec, M_Enceladus)
    g_Saturn3 = gravitational_vector(r_saturn, M_Saturn)
    #a_vec = g_Enceladus3 + g_Saturn3
    a_vec = g_Enceladus3
    return a_vec

def acceleration(r_vec, v_vec, rot_vec):
    R_saturn = R_orbit*np.array([1,0,0])
    r_saturn = R_saturn + r_vec
    g_Enceladus3 = gravitational_vector(r_vec, M_Enceladus)
    g_Saturn3 = gravitational_vector(r_saturn, M_Saturn)

    a_coriolis = coriolis_vector(v_vec, rot_vec)

    a_vec = g_Enceladus3 + g_Saturn3 + a_coriolis
    return a_vec


def gravitational_acceleration(M, R):
    return G_gravitational*M / R**2

day = 3600*24.
#Rotation Period of Enceladus
T_Enceladus = 1.37*day
#Rotation Rate of Enceladus
Omega_Enceladus = 2*pi / T_Enceladus
def coriolis_acceleration(v,theta, rot = Omega_Enceladus):
    return 2*rot*v*sin(theta)

def gravitational_force(M,m,R):
    return G_gravitational*M*m/R**2
def gravitational_vector(r_vec,M): #distance from center:
    R = sqrt(sum(r_vec**2))
    g3 = - r_vec * G_gravitational*M / R**3
    return g3
g_Saturn0 = gravitational_acceleration(M_Saturn, R_orbit)
g_Enceladus = gravitational_acceleration(M_Enceladus, R_enceladus*1e3)


def get_plume_velocity(geyser, velocity, dtheta = 0, dphi = 0):
    x,y,z = pm.geodetic2ecef(geyser.latitude_d, geyser.longitude_d, 0)
    R_vec = np.array([x,y,z])
    azimuth = geyser.azimuth_d + dphi
    elevation = 90 - (geyser.zenith_d + dtheta)
    #print('elevation', elevation, 'azimuth', azimuth)
    e, n, u = pm.aer2enu(azimuth, elevation, 1)

    lat2, lon2, h2 = pm.aer2geodetic(azimuth, elevation, 100, geyser.latitude_d, geyser.longitude_d, 0, ell=enceladus_ellipsoid)
    x2, y2, z2 = pm.geodetic2ecef(lat2, lon2, h2)
    #x2, y2, z2 = pm.aer2ecef(azimuth, elevation, 1, geyser.latitude_d, geyser.longitude_d, 0, ell =enceladus_ellipsoid)

    #x2, y2, z2 = pm.enu2ecef(e, n, u, geyser.latitude_d, geyser.longitude_d, h0=1, ell=enceladus_ellipsoid,deg=True) #TODO check that velocity is properly converted to ECEF
    p_vec = np.array([x2-x, y2-y, z2-z])
    #print('dx',x2-x,'dy',y2-y,'dz', z2-z)
    #print('Plume Vector [e,n,u],',e, n, u)
    #print('plume vector [x,y,z]', p_vec)
    p_vec /= sqrt(sum(p_vec**2))
    #print('plume vector norm [x,y,z]', p_vec)
    #print('')

    #print('velocity =', p_vec*velocity, 'm/s')
    return p_vec*velocity

topo = np.genfromtxt(path_to_directory + 'damascus-topography.txt')
X_topo, Y_topo = topo[:,0], topo[:,1]
X_center = 3235
Infl = 0

def collision_bool(x_particle, h_particle, X_surface, Y_surface, t,  t_wait = 100, tolerance= 50.): #wait 100 s
    dx_arr = x_particle - X_surface
    dy_arr = h_particle - Y_surface

    dr_arr = sqrt(dx_arr**2 + dy_arr**2)
    dr_min = min(dr_arr)
    #print('min dr', dr_min)
    if dr_min < tolerance and t > t_wait:
        return True
    else:
        return False

def particle_trajectory_3(geyser, r_start3, v_start3, h_start, tof, dt=1):  # Define Particle Trajectory given start position and velocity
    # Starting Position and Velocity
    r0 = r_start3  # Position
    v0 = v_start3  # Velocity
    a0 = acceleration_simple(r_start3)  # Acceleration

    Nt = int(tof / dt)  # Number of Time Steps
    tspace = np.linspace(0, tof, Nt)  # Timing Array

    r_list = []
    v_list = []
    a_list = []
    t_list = []

    r_list.append(r_start3)
    v_list.append(v_start3)
    a_list.append(acceleration_simple(r_start3))
    t_list.append(0)

    altitude = []
    arc = []
    altitude.append(h_start)
    displacement = 0
    arc.append(0)
    Ri0 = get_magnitude(r_start3)
    ti = 0
    #print(Nt)
    impact_bool = False
    range_bool = False
    for i in range(Nt):
        ai_0 = a_list[i - 1]
        vi_0 = v_list[i - 1]
        ri_0 = r_list[i - 1]

        dr = vi_0 * dt + 0.5 * ai_0 * dt ** 2
        displacement += dr
        ri = ri_0 + dr
        ai = acceleration_simple(ri)
        vi = vi_0 + ai_0 * dt
        Ri = get_magnitude(ri)

        r_list.append(ri)
        v_list.append(vi)
        a_list.append(ai)

        ti += dt
        t_list.append(ti)

        # CONVERT TO GEODETIC COORDS, calculate horizontal distance
        lat_i, lon_i, h_i = pm.ecef2geodetic(ri[0], ri[1], ri[2], ell=enceladus_ellipsoid)
        #print('geodetic', lat_i, lon_i, h_i)
        dlon = deg2rad(lon_i) - geyser.longitude_r
        dlat = deg2rad(lat_i) - geyser.latitude_r
        a = (sin(dlat / 2)) ** 2 + cos(geyser.latitude_r) * cos(deg2rad(lat_i)) * (sin(dlon / 2)) ** 2
        c = 2 * arctan2(sqrt(a), sqrt(1 - a))
        arc_distance = Ri0 * c
        altitude.append(h_i)
        arc.append(arc_distance)

        t_wait = 100.
        #print(ti, t_wait)
        if (collision_bool(arc_distance, h_i, X_topo - X_center, Y_topo, tspace[i], t_wait) ==True) or ((Ri -Ri0) < 0 and ti > t_wait):
            print(i, 't=', tspace[i], 's, collision! R = ', Ri)
            impact_bool = True
            break
        elif abs(arc_distance) > 5000:
            print(i, 't=', tspace[i], 's, out of range! R =', Ri)
            range_bool = True
            break

    r_list = np.array(r_list)  # 3D vector of particle positions (ECEF coordinates)
    v_list = np.array(v_list)  # 3D vector of particle velocity (ECEF coordinates)
    a_list = np.array(a_list)  # 3D vector of particle acceleration (ECEF coordinates)
    altitude = np.array(altitude)  # Altitude (height above surface)
    arc = np.array(arc)  # Arc Distance along Surface
    return r_list, v_list, a_list, altitude, arc, t_list, impact_bool


h_start = -190
lon0 = geyser_0.longitude_d
lat0 = geyser_0.latitude_d

coords_0 = pm.geodetic2ecef(lat0, lon0, h_start, ell= enceladus_ellipsoid,deg=True)
x0, y0, z0 = coords_0
r_vec_0 = np.array([x0, y0, z0])

vspace = np.linspace(0, 700, N_velocities+1)
Tend = int(1e4) #RUN SIMULATION FOR THIS MANY SECONDS
cmap = matplotlib.cm.get_cmap('gnuplot2')

N_v2 = 100

da = 0.1
deg_min = 0
deg_max = 7.5
N_a2 = int((deg_max - deg_min)/da) + 1
theta_dist = np.linspace(deg_min, deg_max, N_a2) #TODO: Check negative values

x_final = np.zeros((N_v2, N_a2))
n_final = np.zeros((N_radii, N_v2, N_a2))
n_final_sum = np.zeros((N_v2, N_a2))
for i in range(0, N_v2):
    vi = vspace[i]
    for j in range(N_a2):
        #N Paritcles
        n_ij = 0
        for k in range(N_radii):
            n_ijk = (da/15.)*float(geyser_sample2[k,i]) #TODO print the number of particles
            n_final[k,i,j] = n_ijk
            n_ij += n_ijk
        n_final_sum[i,j] = n_ij
        theta_j = theta_dist[j] - geyser_0.zenith_d
        v_vec_0 = get_plume_velocity(geyser_0, vi, theta_j, 0)
        #print('v = ', vi, 'm/s, theta =', theta_dist[j])

        r_list, v_list, a_list, altitude, arc, t_list, impact_bool = particle_trajectory_3(geyser_0, r_vec_0, v_vec_0, h_start, Tend)
        kf = len(r_list) - 1
        r_final = r_list[kf]
        v_final = v_list[kf]
        a_final = v_list[kf]
        altitude_final = altitude[kf]
        arc_final = arc[kf]
        t_final = t_list[kf]

        x_final[i,j] = arc_final
        print('v = ', vi, 'm/s, theta =', theta_j + geyser_0.zenith_d, 'x =', arc_final, 'm, tf =', t_final, 's impact?', impact_bool, 'N =', kf+1)
        print('')

#x_final = np.load('xposition_final.npy')
#print(x_final)
np.save('xposition_3.npy', x_final)
np.save('nposition_3.npy', n_final)
np.save('nposition_sum_3.npy', n_final_sum)

fig = pl.figure(figsize=(8,6), dpi = 100)
ax = fig.add_subplot(111)
for j in range(N_a2):
    color_j = cmap(theta_dist[j] / (deg_max + da))
    ax.plot(vspace[:N_v2], x_final[:,j], label= r'$\theta_{0}$ = ' + str(round(theta_dist[j],2)), c= color_j)
ax.set_xlabel('Starting Velocity $v_{0}$ [m/s]')
ax.set_ylabel(r'Landing Position $R \Delta \theta $ [m]')
ax.legend()
ax.grid()
pl.show()
#pl.close(fig)

'''
x_final = np.load('xposition_2.npy')
n_final = np.load('nposition_2.npy')
n_final_sum = np.load('nposition_sum_2.npy')
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
'''