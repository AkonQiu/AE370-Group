import numpy as np
from src.rk4 import rk4_step
from src.model import flight_model, DEFAULT_PARAMETERS


def simulate(y0, dt, t_final, stop_altitude=None, params=DEFAULT_PARAMETERS):
    """Integrate the five-state IVP; locate downward altitude events by RK4 bisection.

    The final interval is shortened to end at t_final. An initial state at or
    below stop_altitude is returned immediately. All quantities use SI units.
    """
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('dt must be finite and positive')
    if not np.isfinite(t_final) or t_final < 0:
        raise ValueError('t_final must be finite and nonnegative')
    if stop_altitude is not None and not np.isfinite(stop_altitude):
        raise ValueError('stop_altitude must be finite')
    y = np.asarray(y0, dtype=float)
    if y.shape != (5,):
        raise ValueError('y0 must have shape (5,)')
    f = lambda t, state: flight_model(t, state, params)
    f(0.0, y)
    times, states = [0.0], [y.copy()]
    if stop_altitude is not None and y[0] <= stop_altitude:
        return np.array(times), np.array(states)
    for index in range(int(np.ceil(t_final / dt))):
        t = times[-1]
        target = min((index + 1) * dt, t_final)
        step = target - t
        if step <= 0:
            break
        try:
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                y_next = rk4_step(f, t, y, step)
        except (FloatingPointError, ValueError) as exc:
            raise RuntimeError(f'Invalid RK4 state at t={t:g}; reduce dt or check model validity') from exc
        if not np.all(np.isfinite(y_next)):
            raise RuntimeError('Simulation produced nonfinite values; reduce dt')
        event = stop_altitude is not None and y_next[0] <= stop_altitude
        if event:
            low, high = 0.0, step
            for _ in range(50):
                middle = (low + high) / 2
                if rk4_step(f, t, y, middle)[0] > stop_altitude:
                    low = middle
                else:
                    high = middle
            step = (low + high) / 2
            y_next = rk4_step(f, t, y, step)
            target = t + step
        times.append(target)
        states.append(y_next.copy())
        y = y_next
        if event:
            break
    return np.array(times), np.array(states)
