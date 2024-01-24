# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:27:13 2023

@author: Lyuben Davidov
@email: l.davidov@tue.nl
@email: ljubo.davidov@gmail.com
"""

import numpy as np
import matplotlib.pyplot as plt



# PARAMETERS FOR THE SIMULATION
_lambda   = 1.57e-6                             # wavelength [m]
n_g       = 3.85                                # group index
n_eff     = 3.0                                 # effective index of refraction
N_tr      = 8.7e+23                             # transparency carrier density [m^-3]
G         = 0.075                               # confinement factor
A_MQW     = 0.54e-12                            # optical mode area [m^2]
beta      = 0.1e+12                             # internal losses rate [s^-1]
I         = 40e-3                               # injection current [A]
alfa      = 1;                                  # linewidth enhancement factor
r1        = 0.32                                # reflectivity right mirror
r2        = 0.32                                # reflectivity left mirror
c0        = 3e+8                                # speed of light in vacuum [m/s]
omega0    = (2 * np.pi * c0) / _lambda          # angular frequency of light [rad/Hz]
v_g       = c0 / n_g                            # group velocity [m/s]
q         = 1.60217663e-19                      # electron charge [C]
gamma     = 1e+9                                # carrier relaxation rate [Hz]
s         = 0.1*1e+12                           # [Hz]
g         = 17.5*1e+20                          # differential gain coefficient at 1.525um wavelength
V         = 7.56*1e-17                          # active area volume
#N0        = g * G * ( (I)/(q*V) - gamma * N_tr)
g_abs     = 0.1*1e+12 
N0 = 0.025*g_abs**2



# PARAMETERS OF CHOICE
_len      = 5.0
N_seg     = 10
L_seg     = _len * _lambda
dt        = L_seg / v_g
delta     = 1e+12
timesteps = 150000
print("Δt is equal to: {:.3e}s".format(dt));
print("Total simulation time: {:.3e}s".format( dt*timesteps ))






'''SET INITIAL VALUES'''
# SET INITIAL FIELD IN SECTIONS (POSITIVES)
E_p     = np.full((N_seg, timesteps+1), 1e-3, dtype='complex128')

# SET INITIAL FIELD IN SECTIONS (NEGATIVES)
E_n     = np.full((N_seg, timesteps+1), 1e-3, dtype='complex128')

# SET INITIAL CARRIERS IN SECTIONS
N       = np.full((N_seg, timesteps+1), 1e+11)

# CREATE TIMESTEP ARRAY
time    = np.full((timesteps+1,), 0.0)








'''CALC PROCEDURE'''
for counter in range(timesteps):
    print(counter)
    
    # SET BOUNDARY CONDITIONS FOR NEGATIVE PROPAGATION
    E_n[0, counter+1] = np.sqrt(r1) * E_p[0, counter];
    
    # SET BOUNDARY CONDITIONS FOR POSITIVE PROPAGATION
    E_p[N_seg-1, counter+1] = np.sqrt(r2) * E_n[N_seg-1, counter];
    
    
    
    # CALCULATE NEW dt N
    for k in range(N_seg):
        P = np.abs(E_p[k,counter])**2 + np.abs(E_n[k,counter])**2
        dN_dt = N0 - gamma*N[k,counter] - s*N[k,counter]*P
        N[k, counter+1] = N[k,counter] + dt*dN_dt
    
    # CALCULATE NEW dt E+
    for i in range(N_seg-1):
        fac_p1 = (1/v_g) * (- E_p[i+1,counter] + E_p[i, counter])
        fac_p2 = E_p[i, counter]
        fac_p3 = dt * E_p[i, counter] * (-beta/2 + N[i, counter]*(1-1j*alfa)/2)
        E_p[i, counter+1] = fac_p1 + fac_p2 + fac_p3
        
        
    # CALCULATE NEW dt E- 
    for j in range(N_seg-1):
        A_inv = 1/(1/dt - v_g*dt + beta/2 - N[j+1,counter+1]*(1-1j*alfa)/2)
        E_n[j+1,counter+1] = A_inv * E_n[j+1,counter]/dt - A_inv * E_n[j,counter+1]/(dt*v_g)
        
        
    # FILL THE TIMESTEP ARRAY
    time[counter+1] = (counter + 1)*dt;
        
    
        
    


'''CREATE PLOTS'''



# PLOT OF THE INTENSITY
intensity = np.abs(E_n[0,:])**2;
int0 = plt.plot(time, intensity, label="Left-most segment");

intensity = np.abs(E_n[int(N_seg/2)-1,:])**2;
int1 = plt.plot(time, intensity, label="Middle segment");

intensity = np.abs(E_n[N_seg-1,:])**2;
int3 = plt.plot(time, intensity, label="Right-most segment");



plt.xlabel("time [s]")
plt.ylabel("power [W]")
plt.title("Power $|E_n|^{2}$");
#plt.xlim(0.35e-8,0.45e-8)
#plt.ylim(0,0.1e13)
plt.legend();
plt.show();



# PLOT OF THE CARRIER DENSITY
carrier_density = N[0,:];
plt.plot(time, carrier_density);
plt.xlabel("time [s]");
plt.ylabel("carrier density [$m^{-3}$]");
plt.title("Carrier Density N");
#plt.xlim(0.3e-8,0.35e-8);
#plt.ylim(-1.2e+25, 0.1e+25);
plt.show();



# PLOT OF THE REAL PART OF E-
En0 = np.real(E_n[0,:]);
plt.plot(time, En0);
plt.xlabel("time [s]")
plt.ylabel("Re $E_n$ [$\sqrt{W}$]")
plt.title("Real Part of $E_n$");
#plt.xlim(0.45e-7,0.95e-7)
#plt.ylim(0,0.1e13)
plt.show();



# PLOT OF THE IMAGINARY PART OF E-
En1 = np.imag(E_n[0,:]);
plt.plot(time, En1);
plt.xlabel("time [s]")
plt.ylabel("Im $E_n$ [$\sqrt{W}$]")
plt.title("Imag Part of $E_n$");
#plt.xlim(0.45e-7,0.95e-7)
#plt.ylim(0,0.1e13)
plt.show();
    
    
    