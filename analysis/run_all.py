"""Run from the repository root: python -m analysis.run_all."""
from dataclasses import asdict, replace
from pathlib import Path
import csv
import json
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.model import Parameters, flight_model, convective_flux, atmo_density
from src.simulation import simulate

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
FIG = OUT / 'figures'
Y0 = np.array([30000., 2112.14, np.deg2rad(-2.), 300., 1400.])
P = Parameters()
DT = 0.1
HORIZON = 400.


def save_table(name, rows):
    with (OUT / name).open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(name):
    plt.tight_layout()
    for extension in ('png', 'pdf'):
        plt.savefig(FIG / f'{name}.{extension}', dpi=180, bbox_inches='tight')
    plt.close()


def metrics(t, y, p):
    q = convective_flux(y[:, 0], y[:, 1], p)
    peak = int(np.argmax(y[:, 3]))
    return dict(peak_temperature_K=float(y[peak, 3]), peak_time_s=float(t[peak]),
                peak_heat_flux_W_m2=float(q.max()),
                heat_load_J_m2=float(np.trapezoid(q, t)),
                final_altitude_m=float(y[-1, 0]), final_speed_m_s=float(y[-1, 1]),
                end_time_s=float(t[-1]), reached_10km=bool(y[-1, 0] <= 10000. + 1e-6))


def exact_studies():
    """Independent analytic limits of the actual five-equation flight RHS."""
    rows = []
    cooling = replace(P, rho30=0., g0=0.)
    yc = Y0.copy(); yc[2] = 0.; yc[3] = 2500.
    # Amplify drag in this analytic limit so error is above roundoff.
    drag = replace(P, CD=2., g0=0., CL=0., epsilon=0., heating_coefficient=0.)
    yd = Y0.copy(); yd[2] = 0.
    for label, y0, p in [('radiative_cooling', yc, cooling), ('constant_altitude_drag', yd, drag)]:
        previous = None
        steps = [2., 1., .5, .25, .125, .0625, .03125] if label == 'radiative_cooling' else [2., 1., .5, .25, .125]
        for dt in steps:
            t, y = simulate(y0, dt, 20., params=p)
            exact = np.tile(y0, (len(t), 1))
            if label == 'radiative_cooling':
                exact[:, 3] = (y0[3]**-3 + 3*p.epsilon*p.sigma*t/p.CA)**(-1/3)
                error = float(np.max(np.abs(y-exact)))
            else:
                a = .5*atmo_density(y0[0], p)*p.CD*p.A/y0[4]
                exact[:, 1] = y0[1]/(1+a*y0[1]*t)
                error = float(np.max(np.abs(y-exact)))
            order = '' if previous is None else float(np.log2(previous/error))
            rows.append(dict(case=label, dt_s=dt, max_error=error, observed_order=order))
            previous = error
    save_table('exact_convergence.csv', rows)
    drag_rows = [r for r in rows if r['case'] == 'constant_altitude_drag']
    assert all(3.7 < r['observed_order'] < 4.3 for r in drag_rows[1:])
    for label in ['radiative_cooling', 'constant_altitude_drag']:
        r = [row for row in rows if row['case'] == label]
        plt.loglog([x['dt_s'] for x in r], [x['max_error'] for x in r], 'o-', label=label.replace('_', ' '))
    plt.xlabel('Time step (s)'); plt.ylabel('Maximum absolute error (K or m/s)')
    plt.legend(); plt.grid(True, which='both', alpha=.3)
    save_figure('exact_convergence')
    return rows


def flight_convergence():
    # A reference made using our own RK4, not an outside time integrator.
    tr, yr = simulate(Y0, .03125, HORIZON)
    rows = []
    previous = None
    for dt in [1., .5, .25, .125]:
        t, y = simulate(Y0, dt, HORIZON)
        reference = yr[np.rint(t/.03125).astype(int)]
        errors = np.max(np.abs(y-reference), axis=0)
        order = '' if previous is None else float(np.log2(previous/errors[3]))
        rows.append(dict(dt_s=dt, altitude_error_m=errors[0], speed_error_m_s=errors[1],
                         angle_error_rad=errors[2], temperature_error_K=errors[3],
                         temperature_order=order, derivative_calls=4*(len(t)-1)))
        previous = errors[3]
    save_table('flight_convergence.csv', rows)
    assert all(3.7 < r['temperature_order'] < 4.5 for r in rows[1:])
    plt.loglog([r['dt_s'] for r in rows], [r['temperature_error_K'] for r in rows], 'o-', label='Full IVP temperature error')
    h = np.array([r['dt_s'] for r in rows])
    plt.loglog(h, rows[0]['temperature_error_K']*(h/h[0])**4, '--', label='Fourth-order slope')
    plt.xlabel('Time step (s)'); plt.ylabel('Maximum temperature error (K)')
    plt.legend(); plt.grid(True, which='both', alpha=.3); save_figure('flight_convergence')
    return rows


def sweeps():
    rows, checks = [], []
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    studies = [('velocity', np.linspace(1700., 2500., 9), 'Initial speed (m/s)'),
               ('angle', np.linspace(-6., 0., 7), 'Initial flight-path angle (degrees)'),
               ('emissivity', np.linspace(.5, 1., 6), 'Surface emissivity')]
    for column, (name, values, xlabel) in enumerate(studies):
        group = []
        for value in values:
            y0 = Y0.copy(); p = P
            if name == 'velocity': y0[1] = value
            elif name == 'angle': y0[2] = np.deg2rad(value)
            else: p = replace(P, epsilon=float(value))
            t, y = simulate(y0, DT, HORIZON, stop_altitude=10000., params=p)
            fine_t, fine_y = simulate(y0, DT/2, HORIZON, stop_altitude=10000., params=p)
            m, fm = metrics(t, y, p), metrics(fine_t, fine_y, p)
            row = dict(study=name, value=float(value), **m)
            rows.append(row); group.append(row)
            checks.append(dict(study=name, value=float(value),
                               peak_temperature_difference_K=abs(m['peak_temperature_K']-fm['peak_temperature_K']),
                               peak_flux_relative_difference=abs(m['peak_heat_flux_W_m2']/fm['peak_heat_flux_W_m2']-1),
                               heat_load_relative_difference=abs(m['heat_load_J_m2']/fm['heat_load_J_m2']-1)))
        axes[0, column].plot(values, [r['peak_temperature_K'] for r in group], 'o-')
        axes[1, column].plot(values, [r['peak_heat_flux_W_m2']/1e6 for r in group], 'o-')
        axes[2, column].plot(values, [r['heat_load_J_m2']/1e6 for r in group], 'o-')
        for row in range(3):
            axes[row, column].set_xlabel(xlabel); axes[row, column].grid(alpha=.3)
    axes[0, 0].set_ylabel('Peak surface temperature (K)')
    axes[1, 0].set_ylabel('Peak convective flux (MW/m²)')
    axes[2, 0].set_ylabel('Integrated convective heat load (MJ/m²)')
    fig.suptitle('One parameter varied at a time; 400 s horizon or first descent to 10 km')
    save_figure('parameter_sweeps')
    save_table('parameter_sweeps.csv', rows)
    save_table('sweep_refinement.csv', checks)
    assert max(r['peak_temperature_difference_K'] for r in checks) < .5
    assert max(r['peak_flux_relative_difference'] for r in checks) < .001
    assert max(r['heat_load_relative_difference'] for r in checks) < .001
    return rows, checks


def baseline():
    t, y = simulate(Y0, DT, HORIZON, stop_altitude=10000.)
    q = convective_flux(y[:, 0], y[:, 1])
    np.savetxt(OUT/'baseline.csv', np.column_stack((t, y, q)), delimiter=',',
               header='time_s,altitude_m,speed_m_s,flight_path_angle_rad,surface_temperature_K,mass_kg,convective_flux_W_m2', comments='')
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, values, label in zip(axes.flat, [y[:, 0]/1000, y[:, 1], y[:, 3], q/1e6],
                               ['Altitude (km)', 'Speed (m/s)', 'Surface temperature (K)', 'Convective flux (MW/m²)']):
        ax.plot(t, values); ax.set_xlabel('Time (s)'); ax.set_ylabel(label); ax.grid(alpha=.3)
    save_figure('baseline')
    # Exact mechanical energy balance: d(V²/2 + g0 RE h/(RE+h))/dt = -DV/m.
    energy = y[:, 1]**2/2 + P.g0*P.RE*y[:, 0]/(P.RE+y[:, 0])
    power = -.5*atmo_density(y[:, 0])*y[:, 1]**3*P.CD*P.A/y[:, 4]
    energy_residual = float(energy[-1]-energy[0]-np.trapezoid(power, t))
    # Local thermal Jacobian eigenvalue. Flight block may contain positive modes;
    # this is a thermal-mode stability diagnostic, not a full nonlinear proof.
    lam = -4*P.epsilon*P.sigma*y[:, 3]**3/P.CA
    thermal_limit = float(2.785293563/np.max(np.abs(lam)))
    z = DT*lam
    amplification = 1+z+z*z/2+z**3/6+z**4/24
    plt.plot(t, np.abs(amplification))
    plt.xlabel('Time (s)'); plt.ylabel('|R(Δt λ thermal)|'); plt.grid(alpha=.3)
    save_figure('thermal_stability')
    return dict(**metrics(t, y, P), mechanical_energy_residual_J_kg=energy_residual,
                local_thermal_dt_limit_s=thermal_limit,
                max_mass_drift_kg=float(np.max(np.abs(y[:, 4]-Y0[4]))))


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    print('Running baseline and independent exact-solution checks...', flush=True)
    b = baseline(); exact = exact_studies(); convergence = flight_convergence()
    print('Running parameter sweeps with step-halving checks...', flush=True)
    rows, checks = sweeps()
    summary = dict(baseline=b, dt_s=DT, horizon_s=HORIZON,
                   initial_state=Y0.tolist(), parameters=asdict(P),
                   python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__,
                   max_sweep_peak_temperature_difference_K=max(r['peak_temperature_difference_K'] for r in checks),
                   full_ivp_convergence=convergence, exact_limit_convergence=exact)
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(b, indent=2)); print('Results written to', OUT)


if __name__ == '__main__':
    main()
