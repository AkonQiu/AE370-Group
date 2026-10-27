import numpy as np
from src.rk4 import rk4_step

def test_rk4_exp_decay():
    def f(t, y):
        return -y

    y = np.array([1.0])
    t = 0
    dt = 0.1

    y_next = rk4_step(f, t, y, dt)
    exact = np.exp(-dt)
    assert np.isclose(y_next[0], exact, atol=1e-6)