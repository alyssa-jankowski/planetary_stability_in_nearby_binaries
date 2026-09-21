# Import statements
from turtle import pos
import rebound
import numpy as np
import matplotlib.pyplot as plt
import time
import astropy
import pandas as pd
import math
import sys

## Take user input on how many sims to run
# Check that the user gave at least one argument
if (len(sys.argv) < 3) or (len(sys.argv) > 3):
    print("Usage: python myscript.py <lower limit index> <upper>")
    sys.exit(1)
print("withc")
# Take the first CLI argument as ap (you can cast to float/int if needed)
nlow = int(sys.argv[1])
nhigh = int(sys.argv[2])

# Initialize all constant parameters and give sources
JUPITER_MASS_MSUN = 0.00095479191521124043
INITIAL_RADIUS_AU = 3.0

# MASSES SOURCED FROM TAKEDA (2023) - TABLE 2
A_MASS_MSUN = 1.66
B_MASS_MSUN = 1.55


# NEW HAB ZONE CALCULATION

# constants for both eq's
LSUN = 3.828*10**26 # watts
AU_METERS = 1.495979*10**11 # meters
F_EARTH = 1366 # watts*m^-2
F_INNER = F_EARTH * 1.7
F_OUTER = F_EARTH * 0.3

# STAR A HZ
 # Luminosity sourced from Takeda (2023), Table 2
L_A_LSUN = 251
L_STAR_A = L_A_LSUN*LSUN

INNER_HZ_A = np.sqrt(L_STAR_A/(4*math.pi*(F_INNER)))/AU_METERS
OUTER_HZ_A = np.sqrt(L_STAR_A/(4*math.pi*(F_OUTER)))/AU_METERS

# STAR B HZ
 # Luminosity sourced from Takeda (2023), Table 2
L_B_LSUN = 63
L_STAR_B = L_B_LSUN*LSUN

INNER_HZ_B = np.sqrt(L_STAR_B/(4*math.pi*(F_INNER)))/AU_METERS
OUTER_HZ_B = np.sqrt(L_STAR_B/(4*math.pi*(F_OUTER)))/AU_METERS


# VARIABLES THAT MUST BE CHANGED
 # Change to correct host name corresponding to the current host
host_name = 'GamLeoA'
incstar = 45

# Orbit posteriors 
 # Ranges chosen from Romanenko & Kiselev (2014)
a_mid = (554**2)**(1/3)
P_mid = 554 #years
P_err = 27 # years
a_err = 0.1
e_mid = 0.9
e_err = 0.03

# Choose num test particles wanted - This will be num_tests * num_bins - Default to 100 test particles - 10 tests, 10 bins
 # num_tests = number of test particles at each semi major axis step
num_tests = 1
 # num_bins = number of semi major axis steps
num_bins = 100

# Choose the primary around which the test particles will be - 0 for A, 1 for B
primary_choice = 0

# Choose max time and time steps for saving
tmax = 1e6
tstep = 1000 # this is the number output timesteps 

# Semi major axis ranges for each star
semi_major_axes_A = np.linspace(INNER_HZ_A, OUTER_HZ_A, num_bins)
semi_major_axes_B = np.linspace(INNER_HZ_B, OUTER_HZ_B, num_bins)


# PLANET B
# SOURCED FROM HAN 2010 - TABLE 2
GAMMA_LEOb_SMA = 1.19
GAMMA_LEOb_PER = 429./365.25
GAMMA_LEOb_ECCENTRICITY = 0.144
GAMMA_LEOb_ARG_OF_PERIAPSIS = 206.7
# SOURCED FROM TAKEDA (2023) - ABSTRACT - MSINi
GAMMA_LEOb_MASS_MSUN = 0.0102 ## 10.7 MJup

# PLANET C
#  SOURCED FROM HAN 2010 - SECTION 4, SECOND TO LAST PARAGRAPH
GAMMA_LEOc_SMA = 2.6
GAMMA_LEOc_PER = 1340 / 365.25
GAMMA_LEOc_ECCENTRICITY = 0.13
GAMMA_LEOc_MASS_MSUN = 0.00204 ## 2.14 MJup


# BEGIN LOOP FOR EACH POSTERIOR
for index in range(nlow, nhigh):
    # Break loop if there is a maximum desired posterior to run to
    # Initialize sim and set units
    sim = rebound.Simulation()
    sim.units = ('AU', 'years', 'Msun')
    sim.integrator = 'ias15'
    sim.dt = 0.5 / 20 
    sim.collision = "direct" 
    sim.collision_resolve = "merge"
    sim.collision_resolve_keep_sorted = 1
    sim.testparticle_type = 0 
    # Initialize dataframe to store values for each particle
    stored_vals = pd.DataFrame({'body':[],'time':[],'x':[],'y':[],'z': [],'vx':[],'vy':[],'vz':[],'a':[],'e':[],'inc':[],
                            'omega':[],'Omega':[],'f':[], 'starting_a':[], 'primary':[], 'hash':[]})

    # Add Oph A
    sim.add(m=A_MASS_MSUN, x=0, y=0, z=0, vx=0,vy=0,vz=0, hash = "GamLeoA")
    sim.move_to_hel()

    # Add Oph B
    sim.add(m = B_MASS_MSUN, P = np.random.normal(P_mid, P_err), e = np.random.normal(e_mid, e_err), inc = 0* np.pi / 180,
            Omega = 0, omega = 0, f = 0, hash = "GamLeoB")
    # sim.N_active = 2
    particle_index = 0
    # sim.move_to_com()
    # Store Vals for Gam Leo A
    stored_vals.loc[len(stored_vals)] = [particle_index, 0, sim.particles[particle_index].x, sim.particles[particle_index].y, sim.particles[particle_index].z, 
                                             sim.particles[particle_index].vx,sim.particles[particle_index].vy,sim.particles[particle_index].vz,
                                                 np.nan,np.nan,np.nan,np.nan,np.nan, np.nan, np.nan, np.nan, "GamLeoA"]
    particle_index = particle_index+1
    # Store Vals for Gam Leo B
    stored_vals.loc[len(stored_vals)] = [particle_index, 0, sim.particles[particle_index].x, sim.particles[particle_index].y, sim.particles[particle_index].z,
                                             sim.particles[particle_index].vx,sim.particles[particle_index].vy,sim.particles[particle_index].vz,
                                                sim.particles[particle_index].orbit(primary = sim.particles[0]).a, sim.particles[particle_index].orbit(primary = sim.particles[0]).e, 
                                                sim.particles[particle_index].orbit(primary = sim.particles[0]).inc, sim.particles[particle_index].orbit(primary = sim.particles[0]).omega,
                                                sim.particles[particle_index].orbit(primary = sim.particles[0]).Omega, sim.particles[particle_index].orbit(primary = sim.particles[0]).f, np.nan, 0, "GamLeoB"]
    particle_index = particle_index + 1

     #  Add Gamma Leo1b
    sim.add(primary = sim.particles[0], m = GAMMA_LEOb_MASS_MSUN, P= GAMMA_LEOb_PER, e = GAMMA_LEOb_ECCENTRICITY, 
        inc = np.random.rayleigh(0.01), Omega = np.random.uniform(0, 2*np.pi), omega = GAMMA_LEOb_ARG_OF_PERIAPSIS * np.pi / 180, 
        f = np.random.uniform(0, 2*np.pi), hash = "GammaLeo1b")

    # Store Vals for Gamma Leo1b
    stored_vals.loc[len(stored_vals)] = [particle_index, 0, sim.particles[particle_index].x, sim.particles[particle_index].y, sim.particles[particle_index].z,
                                         sim.particles[particle_index].vx,sim.particles[particle_index].vy,sim.particles[particle_index].vz,
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).a, sim.particles[particle_index].orbit(primary = sim.particles[0]).e, 
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).inc, sim.particles[particle_index].orbit(primary = sim.particles[0]).omega,
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).Omega, sim.particles[particle_index].orbit(primary = sim.particles[0]).f, 
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).a, 0, "GammaLeo1b"]
    particle_index = particle_index + 1

    #  Add Gamma Leo1c
    sim.add(primary = sim.particles[0], m = GAMMA_LEOc_MASS_MSUN, P = GAMMA_LEOc_PER, e = GAMMA_LEOc_ECCENTRICITY, 
            inc = np.random.rayleigh(0.01), Omega = np.random.uniform(0, 2*np.pi), omega = np.random.uniform(0, 2*np.pi), 
            f = np.random.uniform(0, 2*np.pi), hash = "GammaLeo1c")

    # Store Vals for Gamma Leo1c
    stored_vals.loc[len(stored_vals)] = [particle_index, 0, sim.particles[particle_index].x, sim.particles[particle_index].y, sim.particles[particle_index].z,
                                         sim.particles[particle_index].vx,sim.particles[particle_index].vy,sim.particles[particle_index].vz,
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).a, sim.particles[particle_index].orbit(primary = sim.particles[0]).e, 
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).inc, sim.particles[particle_index].orbit(primary = sim.particles[0]).omega,
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).Omega, sim.particles[particle_index].orbit(primary = sim.particles[0]).f, 
                                            sim.particles[particle_index].orbit(primary = sim.particles[0]).a, 0, "GammaLeo1c"]
    particle_index = particle_index + 1

    sim.N_active = 3
    sim.move_to_com()
    
    
    # Nested loops to create and add test particles to sim and save their initial values to the dataframe
    for i in range(num_bins):
        for j in range(num_tests):
            # Declare values
            test_e = np.random.uniform(0, 0.1)
            test_inc = incstar * np.pi / 180 + np.random.uniform(0, 0.01)
            test_Omega = np.random.uniform(0, 2*np.pi)
            test_omega = np.random.uniform(0, 2*np.pi)
            test_f = np.random.uniform(0, 2*np.pi)
            # Add to sim with correct corresponding starting semimajor axis depending on the chosen primary
            if primary_choice == 0:
                sim.add(primary = sim.particles[primary_choice], a = semi_major_axes_A[i], e = test_e, inc = test_inc, 
                        Omega = test_Omega, omega = test_omega, f = test_f, hash = "test%s" % (len(sim.particles)-1))
            elif primary_choice == 1:
                sim.add(primary = sim.particles[primary_choice], a = semi_major_axes_B[i], e = test_e, inc = test_inc, 
                        Omega = test_Omega, omega = test_omega, f = test_f, hash = "test%s" % (len(sim.particles)-1))
            particle = sim.particles[len(sim.particles)-1]

            # Save particle's initial values to the dataframe with correct values dependent on the chosen primary
            if primary_choice == 0:
                stored_vals.loc[len(stored_vals)] = [particle_index, 0, particle.x, particle.y, particle.z, particle.vx,particle.vy,particle.vz,
                                                particle.orbit(primary = sim.particles[primary_choice]).a, particle.orbit(primary = sim.particles[primary_choice]).e, 
                                                particle.orbit(primary = sim.particles[primary_choice]).inc, particle.orbit(primary = sim.particles[primary_choice]).omega,
                                                particle.orbit(primary = sim.particles[primary_choice]).Omega, particle.orbit(primary = sim.particles[primary_choice]).f, 
                                                semi_major_axes_A[i], primary_choice, "test%s" % (len(sim.particles)-1)]
            elif primary_choice == 1:
                stored_vals.loc[len(stored_vals)] = [particle_index, 0, particle.x, particle.y, particle.z, particle.vx,particle.vy,particle.vz,
                                                particle.orbit(primary = sim.particles[primary_choice]).a, particle.orbit(primary = sim.particles[primary_choice]).e, 
                                                particle.orbit(primary = sim.particles[primary_choice]).inc, particle.orbit(primary = sim.particles[primary_choice]).omega,
                                                particle.orbit(primary = sim.particles[primary_choice]).Omega, particle.orbit(primary = sim.particles[primary_choice]).f, 
                                                semi_major_axes_B[i], primary_choice, "test%s" % (len(sim.particles)-1)]
            particle_index = particle_index + 1 
    
    # Reset variables
    particle_index = 0
    print(sim.dt)
    sim.move_to_com()
    Eo = sim.energy()
    # Begin timer for running sim integration and start sim integration loop
    start_time = time.time()
    for t in np.linspace(0, tmax, tstep):
        i_time = time.time()
        # Integrate
        sim.integrate(t)
        # Save particles - 0 is star A, 1 is star B, everything else is test particles
        for particle_index in range(len(sim.particles)):
            particle = sim.particles[particle_index]
            if particle_index == 0:
                stored_vals.loc[len(stored_vals)] = [particle_index, t, particle.x, particle.y, particle.z, particle.vx,particle.vy,particle.vz,
                                                     np.nan,np.nan,np.nan,np.nan,np.nan, np.nan, np.nan, np.nan, "GamLeoA"]
            elif particle_index == 1:
                stored_vals.loc[len(stored_vals)] = [particle_index, t, particle.x, particle.y, particle.z, particle.vx,particle.vy,particle.vz,
                                                    particle.orbit(primary = sim.particles[0]).a, particle.orbit(primary = sim.particles[0]).e, 
                                                    particle.orbit(primary = sim.particles[0]).inc, particle.orbit(primary = sim.particles[0]).omega,
                                                    particle.orbit(primary = sim.particles[0]).Omega, particle.orbit(primary = sim.particles[0]).f, np.nan, 1, "GamLeoB"]
            elif particle_index == 2:
                stored_vals.loc[len(stored_vals)] = [particle_index, t, particle.x, particle.y, particle.z, particle.vx,particle.vy,particle.vz,
                                                    particle.orbit(primary = sim.particles[0]).a, particle.orbit(primary = sim.particles[0]).e, 
                                                    particle.orbit(primary = sim.particles[0]).inc, particle.orbit(primary = sim.particles[0]).omega,
                                                    particle.orbit(primary = sim.particles[0]).Omega, particle.orbit(primary = sim.particles[0]).f, np.nan, 0, "GammaLeoAb"]
       #     elif particle_index == 3:
       #         stored_vals.loc[len(stored_vals)] = [particle_index, t, particle.x, particle.y, particle.z, particle.vx,particle.vy,particle.vz,
        #                                            particle.orbit(primary = sim.particles[0]).a, particle.orbit(primary = sim.particles[0]).e, 
        #                                            particle.orbit(primary = sim.particles[0]).inc, particle.orbit(primary = sim.particles[0]).omega,
        #                                            particle.orbit(primary = sim.particles[0]).Omega, particle.orbit(primary = sim.particles[0]).f, np.nan, 0, "GammaLeoAc"]
            else:
                stored_vals.loc[len(stored_vals)] = [particle_index, t, particle.x, particle.y, particle.z, particle.vx,particle.vy,particle.vz,
                                                    particle.orbit(primary = sim.particles[primary_choice]).a, particle.orbit(primary = sim.particles[primary_choice]).e, 
                                                    particle.orbit(primary = sim.particles[primary_choice]).inc, particle.orbit(primary = sim.particles[primary_choice]).omega,
                                                    particle.orbit(primary = sim.particles[primary_choice]).Omega, particle.orbit(primary = sim.particles[primary_choice]).f, np.nan, primary_choice, "test%s" % (particle_index)]
                
        print(time.time() - i_time, time.strftime('%d_%H:%M:%S', time.localtime()), sim.t)
            
    # Print full execution time
    print("Time needed for execution: ", time.time() - start_time, "energy change: ", abs(sim.energy() - Eo) / Eo )
    
    # Save
    stored_vals.to_csv('results/gamleoA/'+ 'withCtest_id_num-' + str(index) + '_' + host_name + time.strftime('ias15_%Y_%m_%d_%H_', time.localtime()) + '_inclination-' + str(incstar) + 'deg.csv', index=False)
    