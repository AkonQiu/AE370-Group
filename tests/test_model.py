import numpy as np
from src.model import atmo_density, gravity, flight_model

def test_density_at_30_km():
    assert np.isclose(atmo_density(30000), 0.0184)

def test_density_decreases_w_altitude():
    assert atmo_density(40000) < atmo_density(30000)
    
def test_gravity_decreases_w_altitude():
    assert gravity(40000) < gravity(30000)
    
def test_flight_model():
    y = np.array([30000.0, 2122.14, np.deg2rad(-2.0), 300.0, 1400.0])
    dydt = flight_model(0.0, y)
    assert dydt.shape == (5,)
    assert np.all(np.isfinite(dydt))