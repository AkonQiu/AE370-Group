import numpy as np
from src.rk4 import rk4_step
from src.model import flight_model

def simulate(y0, dt, t_final, stop_altitude=None):
    """Run the flight simulation using RK4 method"""
    t = 0.0
    y = np.array(y0, dtype=float)
    
    times = [t]
    states = [y.copy()]
    
    while t < t_final:
        step = min(dt, t_final - t)
        y_next = rk4_step(flight_model, t, y, step)
        
        if not np.all(np.isfinite(y_next)):
            raise RuntimeError(
                "Sim produced invalid values, Try a smaller timestep"
            )
        t += step
        y = y_next
        
        times.append(t)
        states.append(y.copy())
        
        if stop_altitude is not None and y[0] <= stop_altitude:
            break
    return np.array(times), np.array(states)