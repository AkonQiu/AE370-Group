"""Generate report-sized figures from the saved analysis CSV files."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def main():
    out = ROOT / 'report' / 'figures'
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'axes.spines.top': False, 'axes.spines.right': False})
    baseline = np.genfromtxt(ROOT / 'results' / 'baseline.csv', delimiter=',', names=True)
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.5), constrained_layout=True)
    for ax, key, scale, ylabel in zip(
        axes.flat,
        ['altitude_m', 'speed_m_s', 'surface_temperature_K', 'convective_flux_W_m2'],
        [1000, 1, 1, 1e6],
        ['Altitude, km', 'Speed, m/s', 'Temperature, K', 'Convective flux, MW/m²']):
        ax.plot(baseline['time_s'], baseline[key]/scale)
        ax.set(xlabel='Time, s', ylabel=ylabel)
        ax.grid(alpha=.3)
    fig.savefig(out / 'baseline.pdf'); plt.close(fig)

    with (ROOT / 'results' / 'parameter_sweeps.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    for study, xlabel in [('velocity', 'Initial speed, m/s'),
                           ('emissivity', 'Emissivity'),
                           ('angle', 'Initial angle, degrees')]:
        group = [r for r in rows if r['study'] == study]
        values = [float(r['value']) for r in group]
        fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.65), constrained_layout=True)
        for ax, key, scale, ylabel in zip(
            axes,
            ['peak_temperature_K', 'peak_heat_flux_W_m2', 'heat_load_J_m2'],
            [1, 1e6, 1e6],
            ['Peak temperature, K', 'Peak flux, MW/m²', 'Heat load, MJ/m²']):
            ax.plot(values, [float(r[key])/scale for r in group], 'o-', markersize=3)
            ax.set(xlabel=xlabel, ylabel=ylabel)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=.3)
        fig.savefig(out / f'{study}_study.pdf'); plt.close(fig)

    c = np.genfromtxt(ROOT/'results'/'flight_convergence.csv', delimiter=',', names=True)
    fig, ax = plt.subplots(figsize=(4.5, 3.1), constrained_layout=True)
    ax.loglog(c['dt_s'], c['temperature_error_K'], 'o-', label='Full-model error')
    ax.loglog(c['dt_s'], c['temperature_error_K'][0]*(c['dt_s']/c['dt_s'][0])**4,
              '--', label='Fourth-order slope')
    ax.set(xlabel='Timestep, s', ylabel='Maximum temperature error, K')
    ax.legend(); ax.grid(which='both', alpha=.3)
    fig.savefig(out/'flight_convergence.pdf'); plt.close(fig)
    print('Report figures written to', out)


if __name__ == '__main__':
    main()
