import os
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from .data_ingest import load_config, fetch_arm_profile, fetch_fluxnet, normalize_to_schema, save_normalized

G = 9.81

def compute_Ri_from_layer(row):
    """Compute Ri_b, Ri_g at z_g and bias B from normalized row."""
    z0 = row['LowerHeight_m']; z1 = row['UpperHeight_m']
    U0 = row['WindSpeed_ms_Z0']; U1 = row['WindSpeed_ms_Z1']
    th0 = row['Temperature_K_Z0']; th1 = row['Temperature_K_Z1']
    L = row.get('ObukhovLength_L', np.nan)
    theta_ref = 0.5*(th0 + th1)
    Dz = z1 - z0
    # bulk Ri_b (guard zeros)
    dU = U1 - U0
    if abs(dU) < 1e-6:
        Ri_b = np.nan
    else:
        Ri_b = (G / theta_ref) * ((th1 - th0) * Dz) / (dU**2)
    # point Ri_g at geometric mean using MOST if L available; fallback to finite diff
    zg = np.sqrt(z0 * z1)
    if not np.isnan(L) and L != 0:
        zeta = zg / L
        # naive phi placeholders; user should replace with calibrated functions
        phi_m = 1.0 + 4.7 * zeta
        phi_h = 1.0 + 7.8 * zeta
        Ri_g_zg = zeta * (phi_h / phi_m**2)
    else:
        # finite difference approx using first layer gradients
        dU_dz = (U1 - U0) / Dz
        dtheta_dz = (th1 - th0) / Dz
        if dU_dz == 0:
            Ri_g_zg = np.nan
        else:
            Ri_g_zg = (G / ((th0+th1)/2.0)) * (dtheta_dz / (dU_dz**2))
    B = Ri_g_zg / Ri_b if (not np.isnan(Ri_b) and Ri_b != 0) else np.nan
    return dict(Ri_b=Ri_b, Ri_g_zg=Ri_g_zg, z_g=zg, B=B)

def run_pipeline_for_arm(site, start, end, layer_lower, layer_upper, out_dir='data/processed'):
    cfg = load_config()
    ds = fetch_arm_profile(site, start, end)
    norm = normalize_to_schema(ds, layer_lower, layer_upper)
    df = pd.DataFrame([norm])
    metrics = df.apply(compute_Ri_from_layer, axis=1, result_type='expand')
    df = pd.concat([df, metrics], axis=1)
    out_path = os.path.join(out_dir, f'consolidated_{site}_{start.strftime("%Y%m%d")}.csv')
    save_normalized(df, out_path)
    return out_path

if __name__ == "__main__":
    # quick CLI demo
    site='sgp'
    start=datetime.utcnow() - timedelta(days=1)
    end=datetime.utcnow()
    out = run_pipeline_for_arm(site, start, end, layer_lower=2.0, layer_upper=10.0)
    print("Wrote:", out)
