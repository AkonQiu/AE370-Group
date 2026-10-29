# Numerical verification notes

## Equations preserved from the premise

With `rho=rho(h)`, `g=g(h)`, `D=rho*V²*CD*A/2`, and
`L=rho*V²*CL*A/2`, the implementation advances:

```
h'     = V sin(gamma)
V'     = -D/m - g sin(gamma)
gamma' = L/(m V) - g cos(gamma)/V
Ts'    = (K sqrt(rho/Rn) V³ - epsilon sigma Ts⁴)/CA
m'     = 0
```

Angle inputs use radians; degree sweeps explicitly convert to radians.
`CA=rho_s*cp*delta`. Environment radiation is neglected. This is a flat-flight
geometry model with altitude-dependent gravity, not spherical reentry dynamics.

## RK4

The stage vectors are `k1=f(t,y)`, `k2=f(t+h/2,y+h*k1/2)`,
`k3=f(t+h/2,y+h*k2/2)`, and `k4=f(t+h,y+h*k3)`.
The update is `y_next=y+h*(k1+2*k2+2*k3+k4)/6`.
For smooth solutions it has O(h^5) one-step error and O(h^4) global error;
there are four RHS evaluations per ordinary step.
Its scalar stability polynomial is `R(z)=1+z+z²/2+z³/6+z⁴/24`.
On the negative real axis, stability extends to approximately `z=-2.78529`.
The code uses no SciPy solver to generate any project result.

## Independent correctness checks

1. Radiative cooling: disable atmosphere and gravity and set gamma=0. All flight
   variables are constant, and `Ts(t)=(Ts0^-3+3*epsilon*sigma*t/CA)^(-1/3)`.
   This checks the actual five-state RHS, radiation sign, heat capacity, and
   temperature integration against a closed-form solution. Larger steps exhibit
   apparent orders above four due to small leading error coefficients; smaller
   steps approach order four before floating-point error dominates.
2. Constant-altitude drag: disable gravity, lift, heating, and radiation; set
   gamma=0. Then `V(t)=V0/(1+a*V0*t)` with `a=rho*CD*A/(2*m)`; altitude,
   temperature, and mass are constant. The convergence experiment uses CD=2
   to expose truncation error above roundoff. A separate test retains CD=0.02.
3. Mechanical energy: define `E=V²/2+g0*RE*h/(RE+h)`. The model satisfies
   `E'=-D*V/m`, because the potential derivative is exactly `g(h)` and lift does
   no work. The test verifies this identity at the derivative level; the baseline
   also records the integral balance residual. That residual includes O(h²)
   trapezoidal quadrature error and is not a fourth-order RK4 error estimate.
4. Mass is invariant. Emissivity changes temperature only, so trajectories for
   different emissivities must match exactly.
5. A nonautonomous vector IVP checks stage time arguments; event and invalid-input
   tests cover simulation control flow independently of baseline results.

## Full-model convergence and timestep selection

`flight_convergence.csv` uses steps 1, 0.5, 0.25, and 0.125 s against our own
RK4 reference at 0.03125 s, on matching time nodes, over 400 s. It records maximum
componentwise errors and observed temperature orders `log2(error_h/error_h/2)`.
This verifies convergence of the coupled calculation; the analytic checks above
supply independent evidence that the equations are implemented correctly.

The production step is 0.1 s. Every sweep case is rerun at 0.05 s and must meet:

- Peak-temperature change < 0.5 K.
- Relative peak-flux change < 0.001.
- Relative integrated-convective-load change < 0.001.

The initial execution gave a worst peak-temperature change below 0.001 K,
substantially tighter than the target. These differences include sampling error
in peak location as well as integration error, and are empirical checks, not
rigorous error bounds. The heat-load integral uses trapezoidal quadrature.

The local thermal eigenvalue is `lambda_T=-4*epsilon*sigma*Ts³/CA`.
The baseline records `2.78529/max(abs(lambda_T))` and plots `|R(dt*lambda_T)|`.
This diagnoses only the thermal mode: the flight block can have positive or
complex eigenvalues, and transient Jacobians change. It is not a proof of full
nonlinear stability. Full-IVP convergence and all-case refinement are the
practical accuracy checks.

## Stopping and validity

The final ordinary step is shortened to match the requested horizon. A downward
altitude crossing is located by 50 bisections, each recomputing a shortened RK4
step from the same left endpoint. No linear state interpolation is used.
This assumes an ordinary step brackets a single relevant crossing; sufficiently
large steps could miss an excursion that crosses and returns within one step.
An initial state at or below the stop altitude returns immediately.

The model requires finite inputs, positive speed and mass, nonnegative surface
temperature, and altitude above -RE. A numerically invalid intermediate state
raises an error rather than returning a misleading result. Model validity at
very low speeds or outside the assumed atmosphere is not established by these
checks. The 400 s analysis does not assert a completed descent or thermal safety.
