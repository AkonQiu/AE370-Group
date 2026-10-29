"""Run from the repository root with python -m analysis.baseline."""
import numpy as np
from src.simulation import simulate


def main():
    y0 = [30000.0, 2112.14, np.deg2rad(-2.0), 300.0, 1400.0]
    times, states = simulate(y0, dt=0.1, t_final=400.0, stop_altitude=10000.0)
    print(f'Time elapsed: {times[-1]:.2f} s')
    print(f'Final altitude: {states[-1, 0]:.2f} m')
    print(f'Peak surface temperature: {np.max(states[:, 3]):.2f} K')
    print(f'Peak velocity: {np.max(states[:, 1]):.2f} m/s')
    print('Stop reason:', '10 km altitude reached' if states[-1, 0] <= 10000.0 + 1e-6 else 'time horizon reached')


if __name__ == '__main__':
    main()
