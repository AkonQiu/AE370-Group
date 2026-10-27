import numpy as np

# Vehicle Parameters:
CD = 0.02
CL = 0.080
A = 4.00 # meters^2
Rn = 0.01 # nose radius, meters

#Thermal properties:
cp = 900.0     # Specific heat capacity, J/(kg*K)
rho_s = 1600.0 # material density, kg/m^3
epsilon = 0.85 # surface emissivity
delta = 0.005  # effective heated thickness, m
sigma = 5.670374419e-8 # Stefan-Boltzmann constant

#Atmosphere parameters:
rho30 = 0.0184  # density at 30km, kg/m^3
H = 7000.0      # Atmospheric scale height, m

#Gravity
g0 = 9.81     # gravitational acceleration, m/s^2
RE = 6.371e6  # Earth's radius, m

#thermal capacity per unit area:
CA = rho_s * cp * delta

#functions

def atmo_density(h):
    return rho30*np.exp(-(h-30000.0)/H)

def gravity(h):
    return g0*(RE/(RE+h))**2

def flight_model(t, y):
    """Return derivatives [dh/dt, dV/dt, dgamma/dt, dTs/dt, dm/dt]"""
    h, V, gamma, Ts, m = y
    
    rho = atmo_density(h)
    g = gravity(h)
    
    D = 0.5 * rho * V**2 * CD * A
    L = 0.5 * rho * V**2 * CL * A
    
    q_conv = 1.7415e-4 * np.sqrt(rho/Rn) * V**3
    
    q_rad = epsilon * sigma * Ts**4
    
    dh_dt = V*np.sin(gamma)
    dV_dt = -D / m - g * np.sin(gamma)
    dgamma_dt = L / (m*V) - g * np.cos(gamma) / V
    
    dTs_dt = (q_conv - q_rad) / CA
    
    dm_dt = 0.0
    
    return np.array([dh_dt, dV_dt, dgamma_dt, dTs_dt, dm_dt])
    