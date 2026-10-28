import numpy as np
from src.simulation import simulate

def test_simulation_advance_time_and_state():
    y0 = [30000.0, 2112.14, np.deg2rad(-2.0), 300.0, 1400.0]
    times, states = simulate(y0, dt=0.01, t_final=0.1)
    
    assert times[-1] == 0.1
    assert states.shape == (len(times), 5)
    assert np.all(np.isfinite(states))
    assert not np.array_equal(states[0], states[-1])