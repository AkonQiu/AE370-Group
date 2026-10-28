import numpy as np
from src.simulation import simulate

y0 = [30000.0, 2112.14, np.deg2rad(-2.0), 300.0, 1400.0]
dt = 0.01
t_final = 400.0

times, states = simulate(y0, dt, t_final, stop_altitude=10000.0)

altitude = states[:, 0]
velocity = states[:, 1]
flight_path_angle = states[:, 2]
surface_temp = states[:, 3]
mass = states[:, 4]

print("Sim completed")
print(f"Time elasped: {times[-1]:.2f} s")
print(f"Final altitude: {altitude[-1]:.2f} m")
print(f"Peak surface temp: {np.max(surface_temp):.2f} K")
print(f"Peak velocity: {np.max(velocity):.2f} m/s")