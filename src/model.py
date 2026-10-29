import numpy as np
from dataclasses import dataclass


@dataclass(frozen=True)
class Parameters:
    """SI parameters; defaults preserve the group's simplified model."""
    CD: float = 0.02
    CL: float = 0.080
    A: float = 4.0
    Rn: float = 0.01
    cp: float = 900.0
    rho_s: float = 1600.0
    epsilon: float = 0.85
    delta: float = 0.005
    sigma: float = 5.670374419e-8
    rho30: float = 0.0184
    H: float = 7000.0
    g0: float = 9.81
    RE: float = 6.371e6
    heating_coefficient: float = 1.7415e-4

    def __post_init__(self):
        for name, value in vars(self).items():
            if not np.isfinite(value) or value < 0:
                raise ValueError(f'{name} must be finite and nonnegative')
        for name in ('Rn', 'cp', 'rho_s', 'delta', 'H', 'RE'):
            if getattr(self, name) == 0:
                raise ValueError(f'{name} must be positive')
        if self.epsilon > 1:
            raise ValueError('epsilon must lie in [0, 1]')

    @property
    def CA(self):
        return self.rho_s * self.cp * self.delta


DEFAULT_PARAMETERS = Parameters()

# Compatibility names for the original module API. For parameter studies,
# pass a Parameters instance instead of modifying these module attributes.
CD, CL, A, Rn = (getattr(DEFAULT_PARAMETERS, k) for k in ('CD', 'CL', 'A', 'Rn'))
cp, rho_s, epsilon, delta, sigma = (
    getattr(DEFAULT_PARAMETERS, k) for k in ('cp', 'rho_s', 'epsilon', 'delta', 'sigma'))
rho30, H, g0, RE = (getattr(DEFAULT_PARAMETERS, k) for k in ('rho30', 'H', 'g0', 'RE'))
CA = DEFAULT_PARAMETERS.CA

def atmo_density(h, params=DEFAULT_PARAMETERS):
    return params.rho30*np.exp(-(h-30000.0)/params.H)

def gravity(h, params=DEFAULT_PARAMETERS):
    return params.g0*(params.RE/(params.RE+h))**2

def flight_model(t, y, params=DEFAULT_PARAMETERS):
    """Return derivatives [dh/dt, dV/dt, dgamma/dt, dTs/dt, dm/dt]"""
    h, V, gamma, Ts, m = y
    if not np.all(np.isfinite(y)) or V <= 0 or m <= 0 or Ts < 0 or h <= -params.RE:
        raise ValueError('State must be finite with V>0, m>0, Ts>=0, h>-RE')

    rho = atmo_density(h, params)
    g = gravity(h, params)

    D = 0.5 * rho * V**2 * params.CD * params.A
    L = 0.5 * rho * V**2 * params.CL * params.A

    q_conv = convective_flux(h, V, params)

    q_rad = params.epsilon * params.sigma * Ts**4

    dh_dt = V*np.sin(gamma)
    dV_dt = -D / m - g * np.sin(gamma)
    dgamma_dt = L / (m*V) - g * np.cos(gamma) / V

    dTs_dt = (q_conv - q_rad) / params.CA

    dm_dt = 0.0

    return np.array([dh_dt, dV_dt, dgamma_dt, dTs_dt, dm_dt])


def convective_flux(h, V, params=DEFAULT_PARAMETERS):
    """Sutton--Graves-style stagnation heat flux in W/m^2."""
    return params.heating_coefficient * np.sqrt(atmo_density(h, params)/params.Rn) * V**3
