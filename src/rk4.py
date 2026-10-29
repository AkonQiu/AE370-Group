import numpy as np

def rk4_step(f, t, y, dt):
    """
    Advance one timestep using the classical four-stage RK4 method
    
    Parameters:
    f: function that computes dy/dt = f(t,y)
    t: float that represents the current time
    y: array that represents the current state vector
    dt: float the timesteps
    
    Returns:
    y_next: array the state vector at the next timestep
    """
    y = np.asarray(y, dtype=float)
    k1 = f(t, y)
    k2 = f(t + dt/2, y + dt * k1 / 2)
    k3 = f(t + dt/2, y + dt * k2 / 2)
    k4 = f(t + dt, y + dt * k3)
    
    y_next = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    return y_next
