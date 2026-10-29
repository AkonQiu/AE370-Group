# Hypersonic descent: AE 370 group project

A five-state, unpowered, simplified hypersonic descent model inspired by the
X-43A. The implementation preserves the group's original equations and baseline
parameters. It is not a validated reconstruction of an X-43A flight.

## Setup and run

Use Python 3.11 or later. From this repository's root:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m analysis.baseline
.\.venv\Scripts\python.exe -m analysis.run_all
```

If using the existing AE370 virtual environment from this checkout, replace
`.\.venv\Scripts\python.exe` with `..\..\.venv\Scripts\python.exe`.
Use module commands as shown so `src` imports resolve. No external time
integrator is used. Running the full analysis overwrites only generated files
under `results/`. Allow roughly a minute, depending on hardware.

## Code organization

- `src/rk4.py`: classical four-stage RK4 implementation.
- `src/model.py`: immutable `Parameters`, density, gravity, heating, and flight RHS.
- `src/simulation.py`: validated time marching, final partial step, and altitude events.
- `analysis/baseline.py`: short console baseline.
- `analysis/run_all.py`: exact-limit checks, full-IVP convergence, 22 parameter
  cases, step-halving verification, CSV/JSON outputs, and PNG/PDF figures.
- `tests/`: unit and regression tests, including independent analytic checks.
- `docs/verification.md`: equations, test rationale, accuracy criteria, and limitations.
- `results/`: generated data, configuration/version metadata, and figures.

## Model and parameters

The state is `[h, V, gamma, Ts, m]`, in meters, m/s, radians, kelvin, and kg.
The baseline is `[30000, 2112.14, -2*pi/180, 300, 1400]`.
The density is `0.0184*exp(-(h-30000)/7000)` kg/m³.
Gravity is `9.81*(6371000/(6371000+h))**2` m/s².
Drag/lift use `0.5*rho*V**2*C*A`. Heating uses
`1.7415e-4*sqrt(rho/Rn)*V**3` W/m². Radiation is `epsilon*sigma*Ts**4`.
See `docs/verification.md` for all five ODEs.

| Parameter | Default | Units / status |
|---|---:|---|
| CD, CL | 0.02, 0.080 | constant assumed coefficients |
| Reference area | 4.00 | m², assumed |
| Nose radius | 0.01 | m, assumed |
| Specific heat | 900 | J/(kg K), assumed constant |
| Material density | 1600 | kg/m³, assumed |
| Emissivity | 0.85 | assumed |
| Heated thickness | 0.005 | m, assumed |
| Areal heat capacity | 7200 | J/(m² K), derived |
| Mass | 1400 | kg, constant; no ablation |

Change parameters per run without modifying globals:

```python
from dataclasses import replace
from src.model import Parameters
from src.simulation import simulate
p = replace(Parameters(), epsilon=0.7)
t, y = simulate([30000, 2112.14, -0.03490658504, 300, 1400],
                dt=0.1, t_final=400, stop_altitude=10000, params=p)
```

## Studies and generated outputs

Each sweep varies only one parameter and retains baseline values for the others:

- Initial speed: 1700–2500 m/s, 9 cases.
- Initial flight-path angle: -6 to 0 degrees, 7 cases.
- Surface emissivity: 0.5–1.0, 6 cases.

All cases use a common 400 s horizon or the first downward crossing of 10 km.
Current cases all reach the time horizon first. Results describe this observation
window, not a complete descent to landing. Peaks are sampled maxima, with their
sampling uncertainty checked by halving the step.

`results/parameter_sweeps.csv` contains peak temperature, peak time, peak flux,
integrated heat load, terminal state, and stopping information.
`results/sweep_refinement.csv` compares every case at 0.1 and 0.05 s.
`results/baseline.csv` contains the complete baseline history.
`results/summary.json` records the actual parameter values, software versions,
baseline metrics, and convergence tables. Figures are saved as PNG and PDF.

The baseline peak temperature is approximately 2632.82 K at 28.0 s; final altitude
is 21.34 km at 400 s. Higher initial speed increases peak temperature over the
studied range. A steeper initial descent increases peak flux and temperature,
but reduces the integrated heat load over this fixed observation window.
Higher emissivity lowers temperature without changing the trajectory or the
incoming convective flux in this one-way coupled thermal model.

## Sources and scope

The project's supplied premise defines the numerical parameters; unverified
vehicle and material values above are assumptions, not NASA measurements.
NASA's [X-43A overview](https://www.nasa.gov/reference/x-43a/) provides motivation.
The heating-law family is described in Sutton and Graves, NASA TR R-376 (1971),
[A general stagnation-point convective-heating equation for arbitrary gas mixtures](https://ntrs.nasa.gov/citations/19720003329).
The numerical coefficient used here is retained from the supplied premise.
A sharp slender vehicle's entire surface is not represented by a stagnation-point
correlation. There is no thrust, ablation, oxidation, conduction into the body,
Mach-dependent aerodynamics, Earth-curvature flight-path term, or material
failure criterion. Numerical verification does not validate these assumptions.

## AI assistance

OpenAI Codex assisted with reviewing the existing code, parameterizing the
model, implementing input/event handling, writing verification tests and study
scripts, and generating figures and reproducibility documentation. The existing
RK4 and governing equations originated in the group's initial code and supplied
premise. Students should review and understand the implementation and results.
