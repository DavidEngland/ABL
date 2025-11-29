# ARM NSA Data Extraction Specification

## Target File: `arm_nsa_20200115.csv`

### Date Range Selection (Winter 2019-2020 Campaign)

**Primary Period:**
- **Date:** January 14-16, 2020 (3-day stable episode)
- **Time:** 00:00 - 23:59 UTC each day
- **Focus:** Nighttime hours (21:00 UTC Jan 14 → 12:00 UTC Jan 16)
- **Rationale:** Arctic winter stable boundary layer with minimal solar forcing

**Alternative Periods (if primary unavailable):**
- **January 2020:** Any 3-day period with clear skies, light winds (<6 m/s), continuous data
- **February 2020:** Feb 10-12, 2020 (known stable period)
- **March 2020:** Mar 1-3, 2020 (transition season, still cold)

**Quality Criteria:**
- Surface temperature < -10°C (stable conditions)
- Wind speed 2-5 m/s (avoids both calm/turbulent extremes)
- Cloud cover < 30% (minimize radiative complications)
- Data completeness > 90% for core variables

---

## ARM NSA Site Information

**Location:**
- **Site Code:** NSA-C1 (Central Facility)
- **Coordinates:** 71.323°N, 156.609°W
- **Elevation:** ~8 m ASL
- **Surface Type:** Tundra/sea ice margin (winter); wet tundra (summer)

**Tower Configuration:**
- **Heights:** 2, 5, 10, 20, 40 m (standard ARM setup)
- **Instruments:**
  - Sonic anemometers (3D wind, virtual temperature)
  - Thermocouples/RTDs (temperature)
  - Cup anemometers (wind speed backup)
  - Radiometers (net radiation, LW/SW components)

---

## Required Variables and ARM Data Streams

### Core Variables (Essential)

| Variable | ARM Data Stream | Measurement Name | Units | Notes |
|----------|-----------------|------------------|-------|-------|
| **Time** | `*.cdf` or `*.nc` | `time` or `time_offset` | seconds since epoch or UTC | Convert to `YYYY-MM-DD HH:MM:SS` |
| **Height levels** | Tower metadata | `height` | m | Fixed: [2, 5, 10, 20, 40] |
| **Temperature** | `nsametC1.b1` | `temp_mean` or `atmos_temp` | K or °C | Convert to K if °C |
| **Wind Speed (U, V components)** | `nsaaosmetC1.a1` | `wspd_vec_mean`, `wdir_vec_mean` | m/s, degrees | Calculate U, V from speed/direction |
| **Friction velocity (u*)** | `nsasurfspecalb1mC1.b1` or eddy cov | `u_star` or derive from `u'w'` | m/s | Surface flux product |
| **Sensible heat flux** | `nsasurfspecalb1mC1.b1` | `sensible_heat_flux` | W/m² | For computing L |
| **Obukhov Length (L)** | Derived or direct | `obukhov_length` or compute | m | Compute if not available |

### Supplementary Variables (Highly Recommended)

| Variable | ARM Data Stream | Notes |
|----------|-----------------|-------|
| **Surface temperature** | `nsametC1.b1` | `temp_sfc` or extrapolate from lowest level |
| **Pressure** | `nsametC1.b1` | `atmos_pressure` (for density, virtual temp) |
| **Relative humidity** | `nsametC1.b1` | `rh_mean` (for virtual temperature correction) |
| **Net radiation** | `nsaradC1.b1` | `net_radiation` (context for stable nights) |
| **Roughness length (z₀)** | Metadata or literature | Typical NSA tundra: 0.0001-0.001 m |

### Derived Variables (Compute Post-Download)

| Variable | Formula | Required Inputs |
|----------|---------|-----------------|
| **U component** | `U = -wspd * sin(wdir * π/180)` | wspd, wdir |
| **V component** | `V = -wspd * cos(wdir * π/180)` | wspd, wdir |
| **Virtual temperature** | `Tv = T(1 + 0.61*q)` where q = specific humidity | T, RH, P |
| **Obukhov Length** | `L = -(u*³ · θv) / (κ · g · w'θ')` | u*, sensible heat flux, θv |
| **Bulk Richardson** | `Ri_b = (g/θ̄) · Δθ · Δz / (ΔU)²` | θ, U at two levels |

---

## Data Streams Priority (ARM Data Archive)

### **Option A: Level 1 Met Data (Easiest, 30-min averages)**
```
nsametC1.b1  → Temperature, wind speed at multiple heights
nsaaosmetC1.a1 → Wind components (if vector averaging available)
```

### **Option B: Eddy Covariance (Best quality, 30-min fluxes)**
```
nsasurfspecalb1mC1.b1 → u*, w'θ', momentum/heat flux
```

### **Option C: High-Frequency (Advanced, requires processing)**
```
nsaecor30mC1.b1 → 20 Hz sonic anemometer data (if you want to compute fluxes yourself)
```

**Recommendation:** Start with **Option A + B** (met + flux products) for quickest turnaround.

---

## Sample CSV Structure (Target Format)

**Filename:** `arm_nsa_20200115.csv`

**Columns:**
```csv
timestamp,layer,z0,z1,T0,T1,U0,U1,V0,V1,u_star,theta_star,L,qc
2020-01-15 03:00:00,0,2.0,10.0,250.5,251.2,2.9,3.4,0.8,1.0,0.15,0.025,25.0,0
2020-01-15 03:00:00,1,10.0,30.0,251.2,252.8,3.4,4.0,1.0,1.5,0.15,0.025,25.0,0
2020-01-15 03:30:00,0,2.0,10.0,250.3,251.0,2.8,3.3,0.7,0.9,0.14,0.028,22.0,0
...
```

**Column Definitions:**
- `timestamp`: UTC time (30-min intervals)
- `layer`: Layer index (0=first layer 2-10m, 1=second layer 10-30m, etc.)
- `z0, z1`: Layer bottom/top heights (m)
- `T0, T1`: Temperature at z0, z1 (K)
- `U0, U1`: U-component wind at z0, z1 (m/s)
- `V0, V1`: V-component wind at z0, z1 (m/s)
- `u_star`: Surface friction velocity (m/s) — same for all layers in a timestep
- `theta_star`: Surface temperature scale (K) — derived from sensible heat flux
- `L`: Obukhov length (m) — same for all layers in a timestep
- `qc`: Quality control flag (0=good, 1=suspect, 2=interpolated, 3=missing)

---

## Data Processing Workflow

### Step 1: Download from ARM Archive

**Web Interface:**
1. Go to https://www.arm.gov/data
2. Search: "NSA" + "MET" + "January 2020"
3. Select: `nsametC1.b1` and `nsasurfspecalb1mC1.b1`
4. Date range: 2020-01-14 to 2020-01-16
5. Format: NetCDF (`.nc`)

**Command Line (if you have ARM Live credentials):**
```bash
# Example using ARM Data Discovery API (requires token)
curl -O "https://adc.arm.gov/armlive/data/...nsametC1.b1.20200115.000000.nc"
```

### Step 2: Extract to CSV (Python)

```python
# filepath: /Users/davidengland/Documents/GitHub/ABL/scripts/arm_nsa_to_csv.py
import xarray as xr
import pandas as pd
import numpy as np

# Load ARM NetCDF files
met = xr.open_dataset('nsametC1.b1.20200115.000000.nc')
flux = xr.open_dataset('nsasurfspecalb1mC1.b1.20200115.000000.nc')

# Extract heights (ARM tower standard)
heights = met['height'].values  # e.g., [2, 5, 10, 20, 40]

# Time vector (30-min intervals)
time = pd.to_datetime(met['time'].values)

# Define layer pairs
layers = [
    (2, 10),   # Layer 0: surface to 10m
    (10, 30),  # Layer 1: 10m to mix of 20+40m levels
    (30, 60)   # Layer 2: elevated (if 40+extrapolation)
]

rows = []
for t_idx, t in enumerate(time):
    # Surface flux values (same for all layers at this timestep)
    u_star = flux['u_star'].isel(time=t_idx).values
    H = flux['sensible_heat_flux'].isel(time=t_idx).values  # W/m²
    
    # Compute theta_star and L
    rho_cp = 1.2 * 1005  # kg/m³ * J/(kg·K) ≈ typical cold conditions
    theta_star = -H / (rho_cp * u_star) if u_star > 0 else np.nan
    
    T_ref = met['temp_mean'].isel(time=t_idx, height=0).values + 273.15  # to K
    g = 9.81
    kappa = 0.4
    
    if u_star > 0 and not np.isnan(theta_star):
        L = -(u_star**3 * T_ref) / (kappa * g * theta_star * u_star)
    else:
        L = np.nan
    
    # Quality control (simple)
    qc = 0 if not np.isnan(L) and abs(L) > 1 else 3
    
    # Loop over layer pairs
    for layer_idx, (z0, z1) in enumerate(layers):
        # Find nearest height indices
        idx0 = np.argmin(np.abs(heights - z0))
        idx1 = np.argmin(np.abs(heights - z1))
        
        T0 = met['temp_mean'].isel(time=t_idx, height=idx0).values + 273.15
        T1 = met['temp_mean'].isel(time=t_idx, height=idx1).values + 273.15
        
        wspd0 = met['wspd_vec_mean'].isel(time=t_idx, height=idx0).values
        wspd1 = met['wspd_vec_mean'].isel(time=t_idx, height=idx1).values
        
        wdir0 = met['wdir_vec_mean'].isel(time=t_idx, height=idx0).values
        wdir1 = met['wdir_vec_mean'].isel(time=t_idx, height=idx1).values
        
        # Convert to U, V components
        U0 = -wspd0 * np.sin(np.radians(wdir0))
        V0 = -wspd0 * np.cos(np.radians(wdir0))
        U1 = -wspd1 * np.sin(np.radians(wdir1))
        V1 = -wspd1 * np.cos(np.radians(wdir1))
        
        rows.append({
            'timestamp': t.strftime('%Y-%m-%d %H:%M:%S'),
            'layer': layer_idx,
            'z0': z0,
            'z1': z1,
            'T0': T0,
            'T1': T1,
            'U0': U0,
            'U1': U1,
            'V0': V0,
            'V1': V1,
            'u_star': u_star,
            'theta_star': theta_star,
            'L': L,
            'qc': qc
        })

# Create DataFrame and save
df = pd.DataFrame(rows)
df.to_csv('arm_nsa_20200115.csv', index=False)
print(f"Saved {len(df)} rows to arm_nsa_20200115.csv")
```

### Step 3: Validate Output

```python
# Quick validation checks
df = pd.read_csv('arm_nsa_20200115.csv')

print("Data Summary:")
print(f"Time range: {df['timestamp'].min()} to {df['timestamp'].max()}")
print(f"Unique layers: {df['layer'].unique()}")
print(f"Temperature range: {df['T0'].min():.1f} - {df['T1'].max():.1f} K")
print(f"L range: {df['L'].min():.1f} - {df['L'].max():.1f} m")
print(f"Quality flags: {df['qc'].value_counts().to_dict()}")

# Check for stable conditions
stable = df[df['L'] > 0]
print(f"\nStable periods: {len(stable)} / {len(df)} records ({100*len(stable)/len(df):.1f}%)")
```

---

## Fallback: Synthetic Data Based on ARM NSA Climatology

If real ARM data is unavailable, generate realistic synthetic data:

```python
# filepath: /Users/davidengland/Documents/GitHub/ABL/scripts/generate_synthetic_arm_nsa.py
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

# Arctic winter climatology
np.random.seed(42)
n_timesteps = 48  # 24 hours at 30-min intervals

times = [datetime(2020, 1, 15, 3, 0) + timedelta(minutes=30*i) for i in range(n_timesteps)]

rows = []
for t in times:
    # Stable night: L ~ 20-50 m, u* ~ 0.1-0.2 m/s
    L = np.random.uniform(20, 50)
    u_star = np.random.uniform(0.12, 0.18)
    theta_star = np.random.uniform(0.02, 0.04)
    
    # Surface layer temperatures (strong inversion)
    T_sfc = 248 + np.random.normal(0, 1.5)  # ~-25°C ± noise
    
    for layer_idx, (z0, z1) in enumerate([(2, 10), (10, 30), (30, 60)]):
        # Log-law temperature profile (simplified)
        z_g = np.sqrt(z0 * z1)
        dT = 0.5 * (z1 - z0) * (z_g / L)**0.4  # Empirical stable gradient
        
        T0 = T_sfc + 0.3 * z0**0.5
        T1 = T0 + dT
        
        # Log-wind profile
        kappa = 0.4
        z0_mom = 0.0001  # Tundra roughness
        U0 = (u_star / kappa) * np.log(z0 / z0_mom)
        U1 = (u_star / kappa) * np.log(z1 / z0_mom)
        
        # Simplified V component (cross-wind)
        V0 = U0 * 0.2
        V1 = U1 * 0.2
        
        rows.append({
            'timestamp': t.strftime('%Y-%m-%d %H:%M:%S'),
            'layer': layer_idx,
            'z0': z0, 'z1': z1,
            'T0': T0, 'T1': T1,
            'U0': U0, 'U1': U1,
            'V0': V0, 'V1': V1,
            'u_star': u_star,
            'theta_star': theta_star,
            'L': L,
            'qc': 0
        })

df = pd.DataFrame(rows)
df.to_csv('arm_nsa_20200115_synthetic.csv', index=False)
```

---

## Contact Points

**ARM Data Archive:**
- **Help Desk:** https://www.arm.gov/user-services
- **Email:** armdata@arm.gov
- **Instrument Mentors:** Check ARM data quality reports for NSA-specific contacts

**UAH Collaborators:**
- **McNider:** Historical ARM NSA campaigns (1990s-2000s knowledge)
- **Biazar:** Recent urban/remote sensing; can advise on data fusion

**Recommended ARM Resource:**
- **NSA Data Quality Reports:** https://www.arm.gov/capabilities/observatories/nsa/data-quality
- **User Guide:** ARM NSA Handbook (search "NSA Handbook" on ARM website)

---

## Success Criteria

**Minimal Viable Dataset:**
- ✅ 24-48 hours of continuous stable conditions (L > 0)
- ✅ At least 3 vertical levels (2, 10, 30+ m)
- ✅ Temperature, wind speed, u*, sensible heat flux
- ✅ <10% missing data (QC flag 0 or 1)

**Ideal Dataset:**
- ✅ 72+ hours (multi-night episode)
- ✅ 5 vertical levels (2, 5, 10, 20, 40 m)
- ✅ All variables including humidity, radiation
- ✅ High-frequency (10 Hz) sonic data for flux validation
- ✅ Independent lidar/radiometer profiles for cross-check

---

**Document Status:** Data acquisition guide v1.0  
**Last Updated:** 2025-01-26  
**Maintainer:** david.england@uah.edu
