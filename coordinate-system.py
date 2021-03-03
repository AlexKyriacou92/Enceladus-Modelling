import numpy as np
#import astropy as ast
from matplotlib import pylab as pl
import math
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D

from numpy import cos, sin, tan, sqrt, arctan
from math import pi
from numpy import deg2rad, rad2deg

import sys
import os

sys.path.append("./pymap3d/src/")
from matplotlib import pyplot as pl
import pymap3d as pm
from enceladus_modelling import *

def longitude_180(longitude_360):
    return (longitude_360 + 180)%360 - 180

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

R_orbit = 237948e3 #Semi-major axis to Saturn in metres, 237,948 km
M_Saturn = 5.6834e26 #Mass of Saturn [kg]

a_Enceladus = 256.6
b_Enceladus = 251.4
c_Enceladus = 248.3


geyser_list = Geyser.load_geyser_list(f_tiger_stripe)
geyser_379 = Geyser.get_geyser(379, geyser_list)
#change to 180 degrees
geyser_379.longitude_d = longitude_180(geyser_379.longitude_d)
geyser_379.longitude_r = deg2rad(geyser_379.longitude_d)

# The local coordinate origin (Zermatt, Switzerland)
lat0 = geyser_379.latitude_d # deg
lon0 = geyser_379.longitude_d  # deg
h0 = 0     # meters

# The point of interest
geyser_377 = Geyser.get_geyser(377, geyser_list)
lat1 = geyser_377.latitude_d  # deg
lon1 = geyser_377.longitude_d   # deg
h1 = 0
enceladus_ellipsoid = pm.Ellipsoid(model = "enceladus")

#coords0 = pm.geodetic2enu(lat1, lon1, h1, lat0, lon0, h0, ell=enceladus_ellipsoid)
coords_379 = pm.geodetic2ecef(lat0, lon0, 0, ell= enceladus_ellipsoid,deg=True)
coords_377 = pm.geodetic2ecef(lat1, lon1, 0, ell= enceladus_ellipsoid,deg=True)
print(coords_379, coords_377)
x0,y0,z0 = coords_379
x1, y1, z1 = coords_377

print(x0-x1,y0-y1,z0-z1)

zenith = geyser_379.zenith_d
elevation = 90 - zenith + 7.5
print(elevation)
azimuth = geyser_379.azimuth_d
slant = 1000 #1000 metres...

e,n,u = pm.aer2enu(azimuth, elevation, slant)
print(e,n,u)
print(sqrt(e**2 + n**2 + u**2))

X2, Y2, Z2 = pm.enu2ecef(e,n,u, geyser_379.latitude_d, geyser_379.longitude_d, h0 = 0, ell = enceladus_ellipsoid, deg =True)
print(X2, Y2, Z2)
print(X2-x0, Y2-y0, Z2-z0)
print(sqrt((X2-x0)**2 + (Y2 -y0)**2 + (Z2 - z0)**2))

vec3d = np.array([X2-x0, Y2-y0, Z2-z0])
norm3d = vec3d/sqrt(sum(vec3d**2))
print(norm3d)

r_vec = np.array([x0,y0,z0])
R = sqrt(sum(r_vec**2))
print(R/1e3)
print(R/(R_enceladus*1e3))

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
v_coriolis = v_escape
print('Gravity of Saturn:', g_Saturn0,'m/s^2, Enceladus:', g_Enceladus, 'm/s^2, coriolis: ', coriolis_acceleration(v_coriolis, deg2rad(90+80)), ' when v =', v_coriolis)

lon_Saturn = 0.
lat_Saturn = 0.
U_S, V_S, W_S = pm.geodetic2ecef(lat_Saturn, lon_Saturn, 0, ell= enceladus_ellipsoid,deg=True)

print('Saturn vector', U_S, V_S, W_S)
g_Saturn3 = g_Saturn0*np.array([1,0,1])

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

h_start = -190
coords_379 = pm.geodetic2ecef(lat0, lon0, h_start, ell= enceladus_ellipsoid,deg=True)
x0, y0, z0 = coords_379
r_vec_379 = np.array([x0, y0, z0])

#DEFINE INITIAL VEOLOCITY

print('g_acceleration at 379', acceleration_simple(r_vec_379))
norm_vector = get_norm(geyser_379.latitude_d, geyser_379.longitude_d,0)

def get_plume_velocity(geyser, velocity, dtheta = 0, dphi = 0):
    x,y,z = pm.geodetic2ecef(geyser.latitude_d, geyser.longitude_d, 0)
    R_vec = np.array([x,y,z])
    azimuth = geyser.azimuth_d + dphi
    elevation = 90 - (geyser_379.zenith_d + dtheta)
    #print('elevation', elevation, 'azimuth', azimuth)
    e, n, u = pm.aer2enu(azimuth, elevation, 1)

    lat2, lon2, h2 = pm.aer2geodetic(azimuth, elevation, 100, geyser.latitude_d, geyser.longitude_d, 0, ell=enceladus_ellipsoid)
    x2, y2, z2 = pm.geodetic2ecef(lat2, lon2, h2)
    #x2, y2, z2 = pm.aer2ecef(azimuth, elevation, 1, geyser.latitude_d, geyser.longitude_d, 0, ell =enceladus_ellipsoid)

    #x2, y2, z2 = pm.enu2ecef(e, n, u, geyser.latitude_d, geyser.longitude_d, h0=1, ell=enceladus_ellipsoid,deg=True) #TODO check that velocity is properly converted to ECEF
    p_vec = np.array([x2-x, y2-y, z2-z])
    print('dx',x2-x,'dy',y2-y,'dz', z2-z)
    print('Plume Vector [e,n,u],',e, n, u)
    #print('plume vector [x,y,z]', p_vec)
    p_vec /= sqrt(sum(p_vec**2))
    print('plume vector norm [x,y,z]', p_vec)
    print('')

    #print('velocity =', p_vec*velocity, 'm/s')
    return p_vec*velocity

print(norm_vector)
r_id = 4
r_grain_example = ice_grain_radii[r_id]
print(r_grain_example)
'''
theta_space = np.linspace(0,90,10)
for i in range(len(theta_space)):
    v_inital = get_plume_velocity(geyser_379, 50, dtheta=theta_space[i], dphi=0)
    print(theta_space[i], v_inital)
'''

speed = 20
dtheta_0 = 1
zenith_0 = geyser_379.zenith_d + dtheta_0
v_initial = get_plume_velocity(geyser_379, speed, dtheta_0)
#TODO: Plot velocity and acceleration

dt = 1
T_final = 650
Nt = int(T_final/dt)
tspace = np.linspace(0,T_final, Nt)
dt = tspace[1] - tspace[0]

r0 = r_vec_379
v0 = v_initial
a0 = acceleration_simple(r_vec_379)
Ri0 = sqrt(sum(r_vec_379**2))
r_list = []
R_list = []
v_list = []
a_list = []
r_list.append(r0)
R_list.append(Ri0)
v_list.append(v0)
a_list.append(a0)

def get_magnitude(vec):
    return sqrt(sum(vec**2))
t_counter =[]
t_counter.append(0)

v_abs = []
a_abs = []
v_abs.append(get_magnitude(v0))
a_abs.append(get_magnitude(a0))
horizontal = []
horizontal.append(0)

h_list = []
h_list.append(h_start)

fig = pl.figure(figsize=(20,10), dpi = 100)
ax = fig.add_subplot(221)
ax1 = ax.twinx()
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
displacement = 0

topo = np.genfromtxt(path_to_directory + 'damascus-topography.txt')
X_topo, Y_topo = topo[:,0], topo[:,1]
X_center = 3235
Infl = 0

def collision_bool(x_particle, h_particle, X_surface, Y_surface, t,  t_wait = 100, tolerance= 50.): #wait 100 s
    dx_arr = x_particle - X_surface
    dy_arr = h_particle - Y_surface

    dr_arr = sqrt(dx_arr**2 + dy_arr**2)
    dr_min = min(dr_arr)
    print('min dr', dr_min)
    if dr_min < tolerance and t > t_wait:
        return True
    else:
        return False

for i in range(Nt-1):
    print(i)
    ai_0 = a_list[i-1]
    vi_0 = v_list[i-1]
    ri_0 = r_list[i-1]

    dr = vi_0*dt + 0.5*ai_0*dt**2
    displacement += dr
    ri = ri_0 + dr
    ai = acceleration_simple(ri)
    vi = vi_0 + ai_0*dt

    r_list.append(ri)
    v_list.append(vi)
    a_list.append(ai)
    v_abs.append(get_magnitude(vi))
    a_abs.append(get_magnitude(ai))

    print(i, 'r=', sqrt(sum(r_list[i]**2))/1e3, 'v =', sqrt(sum(v_list[i]**2)), 'a =', sqrt(sum(a_list[i]**2)))
    Ri = sqrt(sum(r_list[i]**2))

    print(Ri, Ri/(R_enceladus*1e3))
    R_list.append(Ri)

    # CONVERT TO GEODETIC COORDS, calculate horizontal distance TODO: CHECK IF THIS IS ACCURATE
    lat_i, lon_i, h_i = pm.ecef2geodetic(ri[0], ri[1], ri[2], ell=enceladus_ellipsoid)
    print('geodetic', lat_i, lon_i, h_i)
    dlon = deg2rad(lon_i) - geyser_379.longitude_r
    dlat = deg2rad(lat_i) - geyser_379.latitude_r
    a = (sin(dlat / 2)) ** 2 + cos(geyser_379.latitude_r) * cos(deg2rad(lat_i)) * (sin(dlon / 2)) ** 2
    c = 2 * arctan2(sqrt(a), sqrt(1 - a))
    arc_distance = Ri0 * c
    h_list.append(h_i)
    horizontal.append(arc_distance)
    print('Ri0', Ri0,'distance', arc_distance,'height',h_i)
    t_counter.append(tspace[i+1])

    #TODO: Break Simulation once particle lands
    #THIS SHOULD HAPPEN ONCE THE PARTICLE REACHES MINIMUM VELOCITY
    #OR ONCE PARTICLE BREAKS WITH SURFACE
    if collision_bool(arc_distance, h_i, X_topo-X_center, Y_topo, tspace[i]) == True:
        print(i, 't=', tspace[i], 's, collision!')
        break

    '''
    if i%2 == 0:
        dv = get_magnitude(v_list[i]) - get_magnitude(v_list[i-2])

    if dv/dt == 0 and (h_i) < 0:
        print('dv =', dv)
        print('dv/dt = ', dv/dt, 'h_i = ', h_i - h_start)
        print('collison!')
        break

    if dv/dt == 0:
        print('inflections', tspace[i])

    if sqrt(((h_i - any(Y_topo))**2 + (arc_distance - any(X_topo - X_center))**2)) < 1:
        print('Ice Grain Collided')
        break
    '''
R_list = np.array(R_list)
dR = R_list - Ri0


#pl.plot(X_topo-X_center, Y_topo)
#pl.grid()
#pl.show()

suptitle_str = 'Geyser Particle, Launched at v = ' + str(round(speed,2)) + r' m/s, $\theta_{zenith} =$' + str(zenith_0) + ' TOF = ' + str(int(max(t_counter))) + ' s'
fig.suptitle(suptitle_str)

ax.plot(t_counter, np.array(R_list/Ri0), c='b')
ax1.plot(t_counter, h_list)
ax.axhline(1,c='k')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Radial Position R/R_encl')
ax1.set_ylabel('Altitude [m]')
ax.grid()

ax2.plot(t_counter, v_abs)
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Velocity [m/s]')
ax2.grid()

ax3.plot(t_counter, a_abs)
ax3.set_xlabel('Time [s]')

ax3.set_ylabel('Acceleration [m/s^2]')
ax3.grid()

ax4.plot(np.array(horizontal)/1e3, np.array(h_list)/1e3,label='Geyser Trajectory')
ax4.plot((X_topo-X_center)/1e3, Y_topo/1e3, label='Topography')
ax4.grid()
ax4.set_xlabel('Lateral Position [km]')
ax4.set_ylabel('Altitude [km]')
#ax.scatter(horizontal, R_list)
#ax.set_xlabel('Horizontal [m]')
#ax.set_ylabel('Altitude [m]')
print('')

Nf = len(t_counter)
lat_i, lon_i, h_i = pm.ecef2geodetic(r_list[Nf-1][0], r_list[Nf-1][1], r_list[Nf-1][2], ell=enceladus_ellipsoid)
print(lat_i, lon_i, h_i)
dlon = deg2rad(lon_i) - geyser_379.longitude_r
dlat = deg2rad(lat_i) - geyser_379.latitude_r
a = (sin(dlat / 2)) ** 2 + cos(geyser_379.latitude_r) * cos(deg2rad(lat_i)) * (sin(dlon / 2)) ** 2
c = 2 * arctan2(sqrt(a), sqrt(1 - a))
arc_distance = Ri0 * c

print(rad2deg(dlat), rad2deg(dlon))
print('Dispalcement', arc_distance, 't =', T_final)

pl.show()

def particle_trajectory_3(r_start3, v_start3, tof, dt = 1): #Define Particle Trajectory given start position and velocity
    #Starting Position and Velocity
    r0 = r_start3 #Position
    v0 = v_start3 #Velocity
    a0 = acceleration_simple(r_start3) #Acceleration

    Nt = int(tof/dt) #Number of Time Steps
    tspace = np.linspace(0, tof, Nt) #Timing Array

    r_list = []
    R_list = []
    v_list = []
    a_list = []
    t_list = []

    altitude = []
    arc = []

    for i in range(Nt):
        ai_0 = a_list[i - 1]
        vi_0 = v_list[i - 1]
        ri_0 = r_list[i - 1]

        dr = vi_0 * dt + 0.5 * ai_0 * dt ** 2
        displacement += dr
        ri = ri_0 + dr
        ai = acceleration_simple(ri)
        vi = vi_0 + ai_0 * dt

        r_list.append(ri)
        v_list.append(vi)
        a_list.append(ai)

        # CONVERT TO GEODETIC COORDS, calculate horizontal distance
        lat_i, lon_i, h_i = pm.ecef2geodetic(ri[0], ri[1], ri[2], ell=enceladus_ellipsoid)
        print('geodetic', lat_i, lon_i, h_i)
        dlon = deg2rad(lon_i) - geyser_379.longitude_r
        dlat = deg2rad(lat_i) - geyser_379.latitude_r
        a = (sin(dlat / 2)) ** 2 + cos(geyser_379.latitude_r) * cos(deg2rad(lat_i)) * (sin(dlon / 2)) ** 2
        c = 2 * arctan2(sqrt(a), sqrt(1 - a))
        arc_distance = Ri0 * c
        altitude.append(h_i)
        arc.append(arc_distance)

        if collision_bool(arc_distance, h_i, X_topo - X_center, Y_topo, tspace[i]) == True:
            print(i, 't=', tspace[i], 's, collision!')
            break

    r_list = np.array(r_list) #3D vector of particle positions (ECEF coordinates)
    v_list = np.array(v_list) #3D vector of particle velocity (ECEF coordinates)
    a_list = np.array(a_list) #3D vector of particle acceleration (ECEF coordinates)
    altitude = np.array(altitude) #Altitude (height above surface)
    arc = np.array(arc) #Arc Distance along Surface
    return r_list, v_list, a_list, altitude, arc



#==============================================
'''
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
'''


'''
fig = pl.figure(figsize=(10,10), dpi=100)
ax = fig.add_subplot(111,projection='3d')


coefs = (a_Enceladus, b_Enceladus, c_Enceladus)  # Coefficients in a0/c x**2 + a1/c y**2 + a2/c z**2 = 1
# Radii corresponding to the coefficients:
#rx, ry, rz = 1/np.sqrt(coefs)
rx, ry, rz = coefs
# Set of all spherical angles:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

# Cartesian coordinates that correspond to the spherical angles:
# (this is the equation of an ellipsoid):
x = rx * np.outer(np.cos(u), np.sin(v))
y = ry * np.outer(np.sin(u), np.sin(v))
z = rz * np.outer(np.ones_like(u), np.cos(v))

# Plot:
#ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='b', alpha=0.5)

for i in range(len(geyser_list)):
    geyser_i = geyser_list[i]
    if geyser_i.geyser_id != None:

        X, Y, Z =  pm.geodetic2ecef(geyser_i.latitude_d, geyser_i.longitude_d -pi, 0, ell= enceladus_ellipsoid)

        #E, N, U = pm.ecef2enu()
        print('Geyser Position', geyser_i.latitude_d, geyser_i.longitude_d - 180)
        print('Geyser Position [km] ECEF', X / 1e3, Y / 1e3, Z / 1e3)

        E, N, U = pm.ecef2enuv(X, Y, Z, geyser_i.latitude_r, geyser_i.longitude_r)
        print('E', E, 'N', N, 'Up', U)
        #East Coordinate
        plume_c = 1
        u = plume_c*sin(geyser_i.zenith_r) * sin(geyser_i.azimuth_r)
        v = plume_c*sin(geyser_i.zenith_r) * cos(geyser_i.azimuth_r)
        w = plume_c*cos(geyser_i.zenith_r)
        #E2 = E + u
        #N2 = N + v
        #U2 = U + w
        E2 = 0
        N2 = 0
        U2 = 10.

        #X2, Y2, Z2 = pm.enu2ecef(E2, N2, U2, geyser_i.latitude_r, geyser_i.longitude_r - pi, h0 = R_enceladus+10, ell=enceladus_ellipsoid)


        #print('Plume Direction', geyser_i.zenith_d, geyser_i.azimuth_d)
        #print('Plume Vector, E:', u, 'N:', v, "Up:", w)
        #print(E,N,U)
        #print('Shifted Coordinates (ECEF)', X2/1e3, Y2/1e3, Z2/1e3)
        #print((X2-X)/1e3, (Y-Y2)/1e3, (Z2-Z)/1e3)
        #dr = sqrt((X2-X)**2 + (Y2-Y)**2 + (Z2-Z)**2)/1e3
        #print(dr)
        #CREATE LOCAL ENU COORDINATES, for each geyser

        ax.scatter(X/1e3, Y/1e3, Z/1e3, c='r')
        #ax.scatter(X2/1e3, Y2/1e3, Z2, c='g')
        #ax.scatter(X_p/1e3, Y_p/1e3, Z_p/1e3, c='g')
        #a = Arrow3D([X/1e3, U/1e3], [Y/1e3, V/1e3], [Z/1e3, W/1e3], mutation_scale=1, lw=1, arrowstyle="-|>", color="k")
        #ax.add_artist(a)

        #ax.quiver(X/1e3, Y/1e3, Z/1e3, U, V, W, color='k')

ax.view_init(elev=-90, azim=0)
ax.set_xlim(-R_enceladus, R_enceladus)
ax.set_ylim(-R_enceladus, R_enceladus)
ax.set_zlim(-R_enceladus, R_enceladus)
range = 100
#ax.set_xlim(-range, range)
#ax.set_ylim(-range, range)
zrange = 20
#ax.set_zlim(-R_enceladus, -R_enceladus + zrange)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
pl.show()
# note: use np.dot for matrix multplication no *
'''