import numpy as np
import xarray as xr
from datetime import datetime
from .data_ingest import normalize_to_schema
from .data_pipeline import compute_Ri_from_layer

def make_synthetic_profile(z_levels=[1.0,2.0,5.0,10.0], u0=2.0, theta0=270.0):
    z = np.array(z_levels)
    U = u0 + 0.1 * z
    theta = theta0 + 0.2 * z
    ds = xr.Dataset(
        {
            'U': (('z',), U),
            'theta': (('z',), theta)
        },
        coords={'z': z}
    )
    ds.attrs['ObukhovLength'] = 50.0
    return ds

def test_synthetic():
    ds = make_synthetic_profile()
    norm = normalize_to_schema(ds, layer_lower=1.0, layer_upper=2.0)
    metrics = compute_Ri_from_layer(norm)
    assert 'Ri_b' in metrics and 'B' in metrics
    print("Synthetic test metrics:", metrics)

if __name__ == "__main__":
    test_synthetic()
