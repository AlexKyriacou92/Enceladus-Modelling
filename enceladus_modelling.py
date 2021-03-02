import numpy as np
import matplotlib
from matplotlib import pylab as pl
import scipy as sc
import h5py
from matplotlib import colors
from math import pi
import cmasher as cmr
#load all the geysers
from numpy.lib.format import open_memmap
import time
import matplotlib
from numpy import deg2rad, rad2deg, pi, sin, cos

from numpy import rad2deg, deg2rad, sin, cos, arccos
from scipy import ndimage

#Define Constants

path_to_directory = '/home/alex/Enceladus_Data/'
R_enceladus = 252.
# Circumfrence:
C_enceladus = 2 * pi * R_enceladus
C_30deg = 30 / 180 * C_enceladus
Prho = 2 * pi * R_enceladus * (30. / 180.)

#Other UNits
Pixel_Radius = 415
A_ratio = Pixel_Radius / C_30deg
year = 365.25*3600*24
mm = 1e3

C_enceladus = 2 * pi * R_enceladus
rho_ice = 920.

#Enceladus Properties:

g_Enceladus = 0.113 #m/s/s #Enceladus Surface Gravity
V_gas = 700. #m/s #Velocity of Gas Ejecta
R_critical = 0.2e-6 #Critical Radius of Ice Particle
v_escape = 239 #m/s #Escape Velocity Enceladus
M_Enceladus = 1.08e20 #kg, Mass of Enceladus

G_gravitational = 6.674e-11 #m^3kg^-1s^-2 #Gravitational Constant

N_latitude = 180
N_longitude = 360

radius_list_um = [
    0.60208946,
    1.097532,
    1.7218065,
    2.0006604,
    3.0458095,
    7.166058,
    10.2738695
]


Paticle_Spectra = np.load(path_to_directory + 'N_particle_spectrum.npy','r')
N_geysers = 199
N_radii = 7

#Select Geysers in ROI (Region of Interest)
f_hotspot = open(path_to_directory + 'geyser_CIRS_hotspot_weights.txt','r')
f_geysers = open(path_to_directory + 'geysers_positions.txt')
f_alexandria = open(path_to_directory + 'alexandria.txt')
f_baghdad = open(path_to_directory + 'baghdad.txt')
f_cairo = open(path_to_directory + 'cairo.txt')
f_damascus = open(path_to_directory + 'damascus.txt')
thirty_deg_img = path_to_directory + '30-deg-SP-25MB-convert.png.jpg'
f_tiger_stripe = open(path_to_directory + 'tiger_stripe_list.txt','r')

x_Damascus = [-180, -55]
y_Damascus = [30, -160]

x_Baghdad = [-150, 20]
y_Baghdad = [110, -130]

x_Cairo = [-95, 130]
y_Cairo = [150, -120]

x_Alexandria = [-10, 180]
y_Alexandria = [170, 0]

M_Damascus = (y_Damascus[1] - y_Damascus[0]) / (x_Damascus[1] - x_Damascus[0])
C_Damascus = y_Damascus[0] - M_Damascus * x_Damascus[0]

M_Cairo = (y_Cairo[1] - y_Cairo[0]) / (x_Cairo[1] - x_Cairo[0])
C_Cairo = y_Cairo[0] - M_Cairo * x_Cairo[0]

M_Baghdad = (y_Baghdad[1] - y_Baghdad[0]) / (x_Baghdad[1] - x_Baghdad[0])
C_Baghdad = y_Baghdad[0] - M_Baghdad * x_Baghdad[0]

M_Alexandria = (y_Alexandria[1] - y_Alexandria[0]) / (x_Alexandria[1] - x_Alexandria[0])
C_Alexandria = y_Alexandria[0] - M_Alexandria * x_Alexandria[0]

k = 0

#Define Functions
#def heatmap2d(arr: np.ndarray, mode = 'lin', limits = [[0,360],[-90,90]], my_cmap=cmr.arctic, alpha_value = 0.5):
def get_array(filename, prefix = path_to_directory, mode='r'):
    return np.load(prefix + filename, 'r')

def heatmap2d(arr, mode = 'lin', limits = [[0,360],[-90,90]], my_cmap=cmr.arctic, alpha_value = 0.5, vmin0 = None, vmax0 = None):
    arr_pl = arr

    # fig, ax = pl.subplots(figsize=(6,6))
    fig = pl.figure(figsize=(10, 5), dpi=100)
    ax = fig.add_subplot(111)
    # pl.subplot(111,projection="aitoff")
    #lon = np.linspace(0, 360, 360)
    #lat = np.linspace(-90, 90, 180)
    lat = np.linspace(-90, 90, 180) + 0.5
    lon = np.linspace(0, 359,360) + 0.5
    Lon, Lat = np.meshgrid(lon, lat)
    '''
    if bkgd_image != False:
        ax.imshow(bkgd_image)
    '''
    if vmin0 == None and vmax0 == None:
        if (mode == 'log'):
            print('log_plot')

            cmesh = ax.pcolormesh(Lon, Lat, arr_pl, cmap=my_cmap,  norm=colors.LogNorm(), alpha=alpha_value)
        elif (mode == 'lin'):
            print('lin_plot')
            cmesh = ax.pcolormesh(Lon, Lat, arr_pl, cmap=my_cmap, alpha=alpha_value)
        else:
            print('lin_plot')
            cmesh = ax.pcolormesh(Lon, Lat, arr_pl, cmap=my_cmap, alpha=alpha_value)
    else:
        if (mode == 'log'):
            print('log_plot')
            cmesh = ax.pcolormesh(Lon, Lat, arr_pl, cmap=my_cmap, vmin= vmin0, vmax=vmax0, norm=colors.LogNorm(), alpha=alpha_value)
        elif (mode == 'lin'):
            print('lin_plot')
            cmesh = ax.pcolormesh(Lon, Lat, arr_pl, cmap=my_cmap, vmin= vmin0, vmax=vmax0, alpha=alpha_value)
        else:
            print('lin_plot')
            cmesh = ax.pcolormesh(Lon, Lat, arr_pl, cmap=my_cmap, vmin= vmin0, vmax=vmax0, alpha=alpha_value)
    pc = fig.colorbar(cmesh)
    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')

    lat_min = limits[1][0]
    lat_max = limits[1][1]

    lon_min = limits[0][0]
    lon_max = limits[0][1]

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)

    #ax.grid()
    return fig, ax, pc

def r_dist(x0, y0):
    r_sq = x0 ** 2 + y0 ** 2
    r0 = np.sqrt(r_sq)
    return r0

def particle_mass(rho, r):
    return 4 * pi / 3. * rho * r ** 3

def particle_volume(r):
    return (4*pi/3) * r**3

def arc_length(lat0):
    R_encl = 252
    #return ((lat0 + 90.) / 360.) * 2 * pi * R_encl
    return ((lat0 + 90.) / 180.) * 2 * pi * R_encl

def Area(lat_d, lon_d, dlat_d, dlon_d):
    lon_r = deg2rad(lon_d)
    lat_r = deg2rad(lat_d)
    theta_r = pi / 2 - lat_r
    dlat_r = deg2rad(dlat_d)
    dlon_r = deg2rad(dlon_d)
    R_enceladus = 252 * 1e3

    dA = (R_enceladus ** 2) * sin(theta_r) * dlat_r * dlon_r
    return dA



def sum_2d(arr_2):
    N = len(arr_2)
    S = 0
    for i in range(N):
        S += sum(arr_2[i])
    return S

# Creates Geyser Scatter Plot, [lat, lon, depsoition rate]
# You can fill with M_rate, Z_dot_rate or N_flux
def geyser_scatter(arr2, limit = 0):
    lat = deg2rad(np.linspace(-90, 90, 180)) #Latitude radians
    lon = deg2rad(np.linspace(0, 360, 360)) #Longitude radians
    arr_t = arr2
    arr3 = []
    N_before = len(arr2) * len(arr2[0])
    for i in range(180):
        for j in range(360):
            deposit_lvl = arr_t[i, j]
            arr_ij = [lon[j], lat[i], deposit_lvl]
            if deposit_lvl >= limit:
                arr3.append(arr_ij)
    N_after = len(arr3)
    arr4 = np.array(arr3)
    return arr4


def get_coords(lat, lon):
    R = 252
    zenith_deg = 90 + lat
    zenith = deg2rad(90 + lat)
    azimuth = deg2rad(lon)

    prho = (zenith_deg / 180.) * 2 * pi * R
    X = -1 * prho * sin(azimuth)
    Y = prho * cos(azimuth)
    return X, Y


def get_coords_rot(lat, lon):
    R = 252
    zenith_deg = 90 + lat
    zenith = deg2rad(90 + lat)
    azimuth = deg2rad(lon)

    prho = (zenith_deg / 180.) * 2 * pi * R
    # print('Distance from South Pole', prho)
    X = -1 * prho * sin(azimuth + pi)
    Y = prho * cos(azimuth + pi)

    return X, Y


def get_coords_deposition(lat_d, lon_d, deposits, lat_max=-60.):  # lat_d : latitude in degrees
    R_enceladus = 252
    zenith_deg = 90 + lat_d
    zenith = deg2rad(90 + lat_d)
    azimuth = deg2rad(lon_d)

    prho = (zenith_deg / 180.) * 2 * pi * R_enceladus
    #X = -1 * prho * sin(azimuth + pi) #NOTE YOU HAVE TO ADD PI FACTOR HERE
    #Y = prho * cos(azimuth + pi)
    X = -1 * prho * sin(azimuth)
    Y = prho * cos(azimuth)
    X2, Y2, dep2 = [], [], []

    for i in range(len(X)):
        Ri = np.sqrt(X[i] ** 2 + Y[i] ** 2)
        if lat_d[i] <= -60:
            X2.append(X[i])
            Y2.append(Y[i])
            dep2.append(deposits[i])
    X2 = np.array(X2)
    Y2 = np.array(Y2)
    dep2 = np.array(dep2)
    return X2, Y2, dep2


def polar_plot(arr_2d, img = pl.imread(thirty_deg_img), lat_max=-60):
    arr4 = geyser_scatter(arr_2d)
    lat = rad2deg(arr4[:, 1])
    lon = rad2deg(arr4[:, 0])
    depositions = arr4[:, 2]
    X, Y, dep_geysers = get_coords_deposition(lat, lon, depositions, lat_max)

    rotated_img = ndimage.rotate(img, 180)

    fig = pl.figure(figsize=(10, 10), dpi=100)
    ax = fig.add_subplot(111)
    ax.imshow(img, extent=[-img.shape[1] / 2. / A_ratio, img.shape[1] / 2. / A_ratio, -img.shape[0] / 2. / A_ratio,
                           img.shape[0] / 2. / A_ratio])

    sc = ax.scatter(X, Y, s=100, c=dep_geysers, cmap=cmr.arctic, norm=matplotlib.colors.LogNorm(), edgecolors='none',
                    marker='H', alpha=0.5)
    # sc = ax.scatter(X, Y, s=100,  c = dep_geysers,cmap= cmr.arctic,edgecolors='none',marker='H', alpha=0.5)
    cbar = fig.colorbar(sc)

    lon_coords = [0, 90, 180, 270]
    for i in range(len(lon_coords)):
        X_text, Y_text = get_coords_rot(-60, lon_coords[i])
        ax.text(X_text - 15, Y_text, '$' + str(lon_coords[i]) + ' \degree $ W', fontsize=12, c='r')

    lat_coords = [-80, -70, -60]
    for i in range(len(lat_coords)):
        X_text, Y_text = get_coords_rot(lat_coords[i], 215)
        ax.text(X_text, Y_text, '$' + str(-lat_coords[i]) + ' \degree $ S', fontsize=12, c='r')

    ax.set_xlabel('Distance from South Pole X [km]')
    ax.set_ylabel('Distance from South Pole Y [km]')
    ax.set_xlim(-285, 285)
    ax.set_ylim(-285, 285)
    ax.legend()
    ax.grid()
    # pl.show()
    return fig, ax, cbar

#Pixel_Radius = 415
#A_ratio = Pixel_Radius / C_30deg
def make_polar_plot(arr_2d, lat_max, lim_coords = [[-285, 285],[-285, 285]], fig = None, index=111, img = pl.imread(thirty_deg_img), Pixel_Radius = 435, lower_limit = 0, my_alpha = 0.2, my_cmap = 'gist_heat'):
    rotated_img = ndimage.rotate(img, 180)
    A_ratio = Pixel_Radius / C_30deg
    if fig is None:
        fig = pl.figure()
    ax = fig.add_subplot(index)
    arr4 = geyser_scatter(arr_2d, limit=lower_limit)

    lat = rad2deg(arr4[:, 1])
    lon = rad2deg(arr4[:, 0])
    depositions = arr4[:, 2]


    #X, Y, dep_geysers = get_coords_deposition(lat, lon, depositions, lat_max) #GET COORDS DEPOSITON
    X, Y, dep_geysers = get_coords_deposition(lat, lon, depositions, lat_max)  # GET COORDS DEPOSITON
    ax.imshow(img, extent=[-img.shape[1] / 2. / A_ratio, img.shape[1] / 2. / A_ratio, -img.shape[0] / 2. / A_ratio,
                           img.shape[0] / 2. / A_ratio])
    #my_cmap = cmr.arctic
    sc = ax.scatter(X, Y, s=100, c=dep_geysers, cmap=my_cmap, norm=matplotlib.colors.LogNorm(), edgecolors='none',
                    marker='H', alpha=my_alpha)
    # sc = ax.scatter(X, Y, s=100,  c = dep_geysers,cmap= cmr.arctic,edgecolors='none',marker='H', alpha=0.5)
    cbar = fig.colorbar(sc)

    lon_coords = [0, 90, 180, 270]
    for i in range(len(lon_coords)):
        X_text, Y_text = get_coords_rot(-60, lon_coords[i])
        ax.text(X_text - 15, Y_text, '$' + str(lon_coords[i]) + ' \degree $ W', fontsize=12, c='r')

    lat_coords = [-80, -70, -60]
    for i in range(len(lat_coords)):
        X_text, Y_text = get_coords_rot(lat_coords[i], 215)
        ax.text(X_text, Y_text, '$' + str(-lat_coords[i]) + ' \degree $ S', fontsize=12, c='r')

    ax.set_xlabel('Distance from South Pole X [km]')
    ax.set_ylabel('Distance from South Pole Y [km]')

    ax.set_xlim(lim_coords[0][0], lim_coords[0][1])
    ax.set_ylim(lim_coords[1][0], lim_coords[1][1])

    ax.grid()
    return ax, cbar

def show_plot(fig, display_bool = True, time_display_plot = 30):
    if display_bool == True:
        timer = fig.canvas.new_timer(interval=time_display_plot * 1000)  # creating a timer object and setting an interval in millisceconds
        timer.add_callback(close_event)
        timer.start()
        pl.show()
        pl.close(fig)
    else:
        pl.close(fig)


from numpy import arctan2
def inverse_coords(X,Y):
    R = 252.
    zenith = 180*(np.sqrt(X**2 + Y**2))/(2*pi*R)
    Lat = -90 + zenith
    Lon = rad2deg(arctan2(Y,-1*X))
    return Lat, Lon

def inverse_coords_arr(X,Y):
    R = 252.
    zenith = 180*(np.sqrt(X**2 + Y**2))/(2*pi*R)
    Lat = -90 + zenith
    Lon = rad2deg(arctan2(Y,-1*X))
    return np.array([Lat, Lon])


def select_tiger_stripe(x,y,tiger_stripe1, tiger_stripe2):
    X_min1 = tiger_stripe1[0][0]
    X_max1 = tiger_stripe1[0][1]

    Y_min1 = tiger_stripe1[1][0]
    Y_max1 = tiger_stripe1[1][0]

    X_min2 = tiger_stripe2[0][0]
    X_max2 = tiger_stripe2[0][1]

    Y_min2 = tiger_stripe2[1][0]
    Y_max2 = tiger_stripe2[1][0]

    M1 = (Y_max1 - Y_min1) / (X_max1 - X_min1)
    C1 = Y_min1 - M1 * X_min1

    M2 = (Y_max2 - Y_min2) / (X_max2 - X_min2)
    C2 = Y_min2 - M2 * X_min2
    if (M1*x + C1 <= y <= M2*x + C2):
        return x, y
    else:
        return None


'''
if (M_Damascus * X + C_Damascus < Y < M_Baghdad * X + C_Baghdad):
    print(fname, X, Y, 'hello')
    ax.scatter(X, Y, c='b')



def get_heatmap(geyser_id_ii, radius_id, mode = 'lin', file = f_hotspot, geyser_arr = get_array('geyser-map.npy'), show_plot=True):
    j = radius_id

    k = 0
    for line in file:
        col = line.split()
        Ncol = len(col)


        geyser_id = int(col[0])
        lat_d = float(col[1])
        lon_d = float(col[2])
        zenith_d = float(col[3])
        azmiuth_d = float(col[4])
        opening_angle = float(col[5])
        X = float(col[6])
        Y = float(col[7])

        lat_r = rad2deg(lat_d)
        lon_r = rad2deg(lon_d)

        zenith_r = deg2rad(zenith_d)
        azmiuth_r = deg2rad(azmiuth_d)
        U = 10 * sin(zenith_r) * cos(azmiuth_r + pi)
        V = 10 * sin(zenith_r) * sin(azmiuth_r + pi)

        filename = path_to_directory + str(col[8])
        print('Ncol',Ncol)

        if Ncol < 10:
            ii = k
        elif Ncol == 10:
            tiger_stripe = str(col[9])
            ii = k
        elif Ncol > 10:
            ii = int(col[9])
            weight = float(col[10])
        print(geyser_id_ii, geyser_id)
        k += 1


        if int(geyser_id_ii) == int(geyser_id):
            geyser_all_radii = geyser_arr[ii][:]  # Deposition File for all Radii
            radius_label = str(round(radius_list_um[j], 2))
            ice_grain_radius = ice_grain_radii[j]
            geyser_deposition_array = geyser_all_radii[j]

            plot_limits = [[290, 310], [-85, -75]]
            # fig, ax, cbar = heatmap2d(geyser_deposition_array, 'lin', plot_limits)

            fig, ax, cbar = heatmap2d(geyser_deposition_array, mode='lin')
            ax.scatter(lon_d, lat_d, c='r', label=geyser_id)
            plot_info = 'Geyser ' + str(geyser_id) + ', $r_{grain} = ' + radius_label + '\, \mathrm{\mu m} $' + r' $ \theta_{Zen} = \,' + str(zenith_d) + '$, $\phi_{Az} = \,'  + str(azmiuth_d) + '$'
            ax.set_title(plot_info)
            ax.legend()

            if show_plot == True:
                pl.show()
            return fig, ax, cbar
        else:
            print('No matches')

'''




def close_event():
    pl.close() #timer calls this function after 3 seconds and closes the window

ice_grain_radii = np.array(radius_list_um)*1e-6 #List of Particle Radii (in metre, order of magntiude 10^-6)
ice_grain_mass = particle_mass(rho_ice, ice_grain_radii)

class Geyser:
    def __init__(self, id = 0, lat = -90, lon = 0, zen = 0, az = 0, sigma = 15, x = 0, y = 0, geyser_file = '', U = 0, V = 1, element = 0, weight = 1, tiger_stripe='unknown'):
        #super().__setattr__('_callbacks', {})
        self.geyser_id = id
        self.latitude_d = lat
        self.longitude_d = lon
        self.latitude_r = np.deg2rad(lat)
        self.longitude_r = np.deg2rad(lon)
        self.zenith_d = zen
        self.azimuth_d = az
        self.zenith_r = np.deg2rad(zen)
        self.azimuth_r = np.deg2rad(az)
        self.opening_angle_d = sigma
        self.opening_angle_r = np.deg2rad(sigma)
        self.X = x
        self.Y = y
        self.U = U
        self.V = V
        self.file = geyser_file
        self.element = element
        self.weight = weight
        self.radius_list = ice_grain_radii
        self.mass_list = ice_grain_mass
        self.tiger_stripe = tiger_stripe

    def load_geyser_list(self, fpath = f_tiger_stripe):
        geyser_list = []
        k = 0
        for line in fpath:
            col = line.split()
            Ncol = len(col)

            geyser_id = int(col[0])
            lat_d = float(col[1])
            lon_d = float(col[2])
            zenith_d = float(col[3])
            azimuth_d = float(col[4])
            opening_angle = float(col[5])
            X = float(col[6])
            Y = float(col[7])

            zenith_r = deg2rad(zenith_d)
            azimuth_r = deg2rad(azimuth_d)

            U = 10 * sin(zenith_r) * cos(azimuth_r + pi)
            V = 10 * sin(zenith_r) * sin(azimuth_r + pi)

            filename = path_to_directory + str(col[8])
            print('Ncol', Ncol)
            weight = 1
            tiger_stripe_id = str(col[9])
            geyser_1 = Geyser(geyser_id, lat_d, lon_d, zenith_d, azimuth_d, opening_angle, X, Y, filename, U, V, k, weight, tiger_stripe_id)
            geyser_list.append(geyser_1)
            k += 1
        return geyser_list
    def get_radius(radius_id):
        return ice_grain_radii[radius_id]

    def get_mass(mass_id):
        return ice_grain_mass[mass_id]

    def get_geyser(geyser_id, geyser_list):
        N_geyser_list = len(geyser_list)
        geyser_select = geyser_list[0]
        for i in range(N_geyser_list):
            geyser_i = geyser_list[i]
            geyser_id_i = geyser_i.geyser_id
            if geyser_id == geyser_id_i:
                geyser_select = geyser_list[i]
        return geyser_select

    def get_heatmap(self, radius_id, geyser_list, geyser_memmap):
        geyser_element = self.element
        heatmap_arr = geyser_memmap[geyser_element][radius_id]
        return heatmap_arr

    def plot_heatmap(self, radius_id, geyser_list, geyser_memmap, mode='lin', show_plot = True):
        plot_limits = [[290, 310], [-85, -75]]
        # fig, ax, cbar = heatmap2d(geyser_deposition_array, 'lin', plot_limits)
        heatmap_arr = self.get_heatmap(radius_id, geyser_list, geyser_memmap)
        fig, ax, cbar = heatmap2d(heatmap_arr)

        ax.scatter(self.longitude_d, self.latitude_d, c='r', label=str(self.geyser_id))
        plot_info = 'Geyser ' + str(
            self.geyser_id) + ', $r_{grain} = ' + str(self.radius_list[radius_id]) + '\, \mathrm{\mu m} $' + r' $ \theta_{Zen} = \,' + str(
            self.zenith_d) + '$, $\phi_{Az} = \,' + str(self.azimuth_d) + '$'
        ax.set_title(plot_info)
        ax.legend()

        if show_plot == True:
            pl.show()
        return fig, ax, cbar
    '''
    def make_polar_plot(self):
    '''
'''
class Geyser_Array:
    def __init__(self, items=None):
        if items is None: # watch out for the mutable default argument
            items = []
        self.items = items

    def add_geyser(self, Geyser): # use a better name
        self.items.append(Geyser) # don't create a needless intermediate, single-element list

    def get_attribute(self, attr, geyser_arr):
        geysers = [geyser_arr[i] for i in range(len(geyser_arr))]
        attribute_list = [Geyser.attr for Geyser in geysers]
        return attribute_list
'''