from enceladus_modelling import *
from numpy import sqrt, arcsin
from scipy.integrate import quad

geyser_list = Geyser.load_geyser_list(f_tiger_stripe)
geyser_379 = Geyser.get_geyser(379, geyser_list)

#Velocity Probability Distribution
def p(v,r,r_c = R_critical,v_gas = V_gas):
    return (1 + r/r_c) * (r/r_c) * (v/(v_gas**2)) * (1 - (v/v_gas))**((r/r_c) - 1)


#This calculates the distance travelled -> assuming constant acceleration
#Distance Ice Particle will fly from the geyser
def delta_x(v, theta, g = g_Enceladus): #Landing Position of Geyser Particle
    #return 2*(np.array(v)**2)*cos(theta)*sin(theta) / g
    return 2*(np.array(v)**2)*sin(2*theta) / g

def delta_y(v, theta, g=g_Enceladus):
    vy = np.array(v)*sin(theta)
    return np.array(vy)**2 / (2*g)

def time_of_flight(v, theta, g = g_Enceladus):
    vy = np.array(v)*sin(theta)
    return np.array(vy)/g
def count_population(x, xmin, xmax): #count number of elements of x that are within xmin and xmax
    pop = 0
    for i in range(len(x)):
        if x[i] >= xmin and x[i] <= xmax:
            pop += 1
    return pop
nbins = 700
N_vel = 1000
v_space_linear = np.linspace(0, V_gas, N_vel)


#f_reduction = 1e4 #Reduce Sample to Managable Number
#N_sample = int(N_particles/f_reduction)
cmap = matplotlib.cm.get_cmap('gnuplot2')

Nrate_geyser_heatmap = get_array('Nrate_surface.npy', path_to_directory, 'r+')
Nradius_spectrum = np.load(path_to_directory + 'N_particle_spectrum.npy', 'r')
#radii = geyser_0.radius_list
N_collisions = Nradius_spectrum[:,2]
N_launched = Nradius_spectrum[:,1]
#P_collsion = N_collsions/N_launched #probability of landing inside radius_rang

#Distribution of Velocity for a Geyser!

path_to_gphysics = path_to_directory + 'Geyser_Phyics/'
#geyser_velocities_launch = make_memmap
f_ratio = 1e9 #COntrols number of particles simulated -> reduce for speed, increase for smoother distributions

def geyser_position(v, theta, t, a = -g_Enceladus):
    vx = v*sin(theta)
    vy = v*cos(theta)
    x = vx*t
    y = vy*t + 0.5*a*t**2
    return x, y

def gravtional_acceleration(h, M = M_Enceladus, R = R_enceladus*1e3):
    R_local = R + h #The Surface Radius + the height above the surface
    g_local = G_gravitational * M / R_local**2
    return g_local

def get_angle(x_max, v_max, g = g_Enceladus):
    return arcsin(g * x_max / (2*v_max**2))/2.

def get_velocity(x_max, theta, g = g_Enceladus):
    return sqrt(x_max * g / (2*sin(2*theta)))

#1st Step: Create Velocity Distributions:
#for r in range(N_radii):

#FIND MAXIMUM ANGLE AND VELOCITY:
X_max = 2000. #Maximum of 2000 metres from the geyser

for r in range(N_radii):
    fig = pl.figure(figsize=(8, 8), dpi=160)
    ax = fig.add_subplot(211)
    ax2 = ax.twinx()

    ax_g = fig.add_subplot(212)
    ax_g2 = ax_g.twinx()

    #GET NUMBER OF PARTICLES
    heatmap_379 = geyser_379.get_heatmap(r, geyser_list, Nrate_geyser_heatmap)
    # GET LOCAL COORDINATE
    latitudes = np.linspace(-90, 90, 180)
    longitudes = np.linspace(0, 360, 360)

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
    v_dist = np.array(v_dist)
    N_survivors = len(v_dist)
    print('N Survivors', N_survivors, 'fraction = ', 100*float(N_survivors)/float(N_sample), '%')
    color_i = cmap(r_grain/12e-6)
    N_pop = count_population(np.array(v_dist), 0, v_escape)
    print('fraction of particles v < v_escape', round(float(N_pop)/float(N_survivors)*100,2),'%')

    f_landing = round(float(N_pop) / float(N_survivors) * 100, 2)
    hist_label = '$r_{grain} = ' + str(round(ice_grain_radii[r] * 1e6, 2))
    # + ' \, \mathrm{\mu m}$, $\dot{N}_{particles} =  ' + str(round(N_particles / 1e12, 2)) + r' \times 10^{12} \, \mathrm{s^{-1}}$, $f_{landed} = ' + str(f_landing) + '$ %'


    #ax.hist(v_dist, bins=nbins, color=color_i, label=hist_label, histtype='step', linestyle=('solid'))
    zenith = geyser_379.zenith_r
    opening_r = deg2rad(geyser_379.opening_angle_d)/2
    zenith_min = zenith - opening_r
    zenith_max = zenith + opening_r
    #theta_dist = np.random.uniform(zenith_min, zenith_min, N_sample-1)

    #ax.hist(v_dist, bins=nbins, color=color_i, label=hist_label, histtype='step', linestyle=('solid'))
    #y_dist = delta_y(v_dist, theta_dist)
    #ax.hist(y_dist, bins=nbins, color=color_i, label=hist_label, histtype='step', linestyle=('solid'))
    dt = 1e3
    N_angles = 2

    theta_space = np.linspace(zenith_min, zenith_max, N_angles)
    #SCAN OVER ANGLES
    #g_list = []
    for k in range(1, N_angles):
        #k_rand = np.random.randint(0, N_sample -1)
        #theta_k = theta_dist[k_rand]
        theta_k = theta_space[k]
        print('zenith angle', theta_k)
        Nt = 40
        #vy_dist = v_dist*cos(theta_k)
        #dt_max = max(vy_dist)/g_Enceladus

        #SCAN OVER VELOCITIES
        N_vel = 20
        #for l in range(N_vel):

        V_max = get_velocity(X_max, theta_k) #Calcuate Maximum Velocity given Outer Most Angle


        for l in range(1):
            velocity_median = np.median(v_dist)

            Theta_max = rad2deg(get_angle(X_max, velocity_median))
            print('Median Velocity: v_median =', velocity_median, 'm/s')
            print('Maximum Angle from Zenith given v_median', Theta_max)
            print('Maximum Velocity given theta = ', rad2deg(theta_k), 'v = ', V_max, 'm/s')
            #ll = np.random.randint(0, N_sample)
            #v_Launch = v_dist[ll]

            v_Launch = velocity_median
            #v_Launch = max(v_dist)
            vy_Launch = v_Launch*cos(theta_k)
            dt_max = 2*vy_Launch/g_Enceladus
            tspace = np.linspace(0, dt_max, Nt)
            print(k,l,'max time of flight', dt_max)

            vy_Launch2 = velocity_median*cos(Theta_max)
            dt_max2 = 2*vy_Launch2/g_Enceladus
            tspace2 = np.linspace(0, dt_max2, Nt)

            vy_Launch3 = V_max*cos(theta_k)
            dt_max3 = 2*vy_Launch2/g_Enceladus
            tspace3 = np.linspace(0, dt_max3, Nt)

            #for m in range(Nt)
            x_dt, y_dt = geyser_position(v_Launch, theta_k, tspace)
            x_dt2, y_dt2 = geyser_position(velocity_median, Theta_max, tspace2)
            x_dt3, y_dt3 = geyser_position(V_max, theta_k, tspace3)
            #t_landing

            #if m == 0:
            label1 = 'v =' + str(round(v_Launch,2)) + r' m/s, $\theta = ' + str(round(rad2deg(theta_k),2)) + '$, $\Delta t = $' + str(int(dt_max)) + ' s'
            label2 = 'v =' + str(round(velocity_median,2)) + r' m/s, $\theta = ' + str(round(rad2deg(Theta_max),2)) + '$, $\Delta t = $ ' + str(int(dt_max2)) + ' s'
            label3 = 'v =' + str(round(V_max,2)) + r' m/s, $\theta = ' + str(round(rad2deg(theta_k),2)) + '$, $\Delta t = $ ' + str(int(dt_max3)) + ' s'
            #else:
            #label1 = None
            #label2 = None
            #label3 = None


            ax.scatter(x_dt/1e3, y_dt/1e3, c = 'b', label=label1)
            ax.scatter(x_dt2/1e3, y_dt2/1e3, c='r',  label=label2)
            ##ax.scatter(x_dt3/1e3, y_dt3/1e3, c= 'g',  label=label3)
            ax2.scatter(x_dt/1e3, (y_dt / 1e3 + R_enceladus) / R_enceladus, c='b')

            ax_g.scatter(x_dt/1e3, gravtional_acceleration(y_dt), c= 'b',  label=label1)
            #ax_g.scatter(x_dt2 / 1e3, gravtional_acceleration(y_dt2), c='r',label=label2)
            #ax_g.scatter(x_dt3 / 1e3, gravtional_acceleration(y_dt3), c='g',label=label3)

            ax_g2.scatter(x_dt/1e3, gravtional_acceleration(y_dt)/g_Enceladus, c= 'b')
            #ax.scatter(x_dt, y_dt, c=cmap(dt/dt_max))
    print('')

    title_label = 'Ice Grain Trajectory (Geyser 379), $ \lambda = $' + str(-geyser_379.latitude_d) + ' S, $ \phi = $ ' + str(geyser_379.longitude_d) + '\n$r_{grain} = ' + str(round(r_grain*1e6,2)) + '\, \mathrm{\mu m}$ \n$v = ' + str(round(v_Launch,2)) + r'\, m/s $, $\theta_{zenith} = ' + str(round(rad2deg(theta_k),2))+ '$'
    ax.set_title(title_label)
    #ax.set_title('Geyser Particle Trajectory')
    ax.set_xlabel('Lateral Distance X [km]')
    ax.set_ylabel('Altitude Y [km]')
    ax_g.set_ylabel('Gravitational Acceleration g [m/s^2]')
    ax_g.set_xlabel('Lateral Distance X [km]')
    ax2.set_ylabel('$R/R_{Enceladus}$')
    ax_g2.set_ylabel('$g/g_{surface}$')

    ax_g.grid()
    ax.legend()
    ax.grid()
    pl.show()
#NEXT STEP: GET GRAVITATION

#THESIS PLOT============================================================================================================
'''
fig = pl.figure(figsize=(16,8), dpi = 160)
ax = fig.add_subplot(111)
for r in range(N_radii):
    #GET NUMBER OF PARTICLES
    heatmap_379 = geyser_379.get_heatmap(r, geyser_list, Nrate_geyser_heatmap)
    # GET LOCAL COORDINATE
    latitudes = np.linspace(-90, 90, 180)
    longitudes = np.linspace(0, 360, 360)

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
    color_i = cmap(r_grain/12e-6)
    N_pop = count_population(np.array(v_dist), 0, v_escape)
    print('fraction of particles v < v_escape', round(float(N_pop)/float(N_survivors)*100,2),'%')
    f_landing = round(float(N_pop)/float(N_survivors)*100,2)
    hist_label ='$r_{grain} = ' + str(round(ice_grain_radii[r] * 1e6, 2)) + ' \, \mathrm{\mu m}$, $\dot{N}_{particles} =  ' + str(round(N_particles / 1e12, 2)) + r' \times 10^{12} \, \mathrm{s^{-1}}$, $f_{landed} = ' + str(f_landing) + '$ %'
    ax.hist(v_dist,bins=nbins,  color=color_i, label=hist_label, histtype='step',linestyle=('solid'))
    print('')

    v_hist, v_bins = np.histogram(v_dist, bins = nbins, range(0, V_gas))

ax.set_title('Velocity Distribution of Plume Particles, Geyser 379 in Damascus Sulcus')
ax.set_xlabel('Ice grain launching velocity $v_{grain}$ [$\mathrm{m s^{-1}}$] ')
ax.axvline(v_escape,c='k', label='Encealadus Escape Velocity' + '\n$v_{escape} = ' + str(v_escape) + '\, \mathrm{m s^{-1}} $')
ax.set_ylabel(r'Particles per velocity bin ($\Delta v/\mathrm{bin} = 1 m/s$) $N_{particles} \times 10^{' + str(int(np.log10(f_ratio))) + '}$')
ax.legend()
#ax.set_yscale('log')
ax.grid()
pl.show()
'''
#END THESIS PLOT========================================================================================================


