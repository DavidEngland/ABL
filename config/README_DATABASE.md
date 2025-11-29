# ABL Database Usage Guide

## Quick Start

**PostgreSQL:**
```bash
# Create database
createdb abl_database

# Load schema
psql -d abl_database -f schema_complete.sql

# Load sample data
psql -d abl_database -f sample_data.sql

# Query ML features
psql -d abl_database -c "SELECT * FROM ML_TrainingSet LIMIT 10;"
```

**SQLite (simpler, no server needed):**
```bash
# Create database
sqlite3 abl.db < schema_complete.sql
sqlite3 abl.db < sample_data.sql

# Query
sqlite3 abl.db "SELECT * FROM Bias_Summary;"
```

---

## Import Real Data (CSV → SQL)

### Python Helper Script

```python
# filepath: /Users/davidengland/Documents/GitHub/ABL/scripts/csv_to_sql.py
import pandas as pd
import psycopg2
from datetime import datetime

# Load ARM NSA CSV (example format)
df = pd.read_csv('arm_nsa_20200115.csv')

# Connect to database
conn = psycopg2.connect("dbname=abl_database user=postgres")
cur = conn.cursor()

# Insert layers
for _, row in df.iterrows():
    obs_id = f"ARM_NSA_{row['timestamp'].replace(' ','_')}_L{row['layer']}"
    cur.execute("""
        INSERT INTO Layers (
            ObservationID, SiteID, TimestampUTC, LayerIndex,
            LowerHeight_Z0, UpperHeight_Z1,
            Temperature_K_Z0, Temperature_K_Z1,
            WindSpeed_ms_Z0, WindSpeed_ms_Z1,
            U_Star_FrictionVel, ObukhovLength_L, QualityFlag
        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (ObservationID) DO NOTHING
    """, (
        obs_id, 'ARM_NSA_C1', row['timestamp'], row['layer'],
        row['z0'], row['z1'], row['T0'], row['T1'],
        row['U0'], row['U1'], row['u_star'], row['L'], row['qc']
    ))

conn.commit()
conn.close()
```

### Direct SQL COPY (fastest for large files)

```sql
-- ARM NSA tower data (pre-formatted CSV)
\COPY Layers(ObservationID, SiteID, TimestampUTC, LayerIndex, 
             LowerHeight_Z0, UpperHeight_Z1, Temperature_K_Z0, Temperature_K_Z1,
             WindSpeed_ms_Z0, WindSpeed_ms_Z1, U_Star_FrictionVel, ObukhovLength_L)
FROM 'arm_nsa_winter_2020.csv' CSV HEADER;

-- SHEBA (add via Python for timestamp parsing)
```

---

## ML Training Workflow

### Export Training Set

```sql
-- PostgreSQL
\COPY (SELECT * FROM ML_TrainingSet WHERE CampaignName='SHEBA') 
TO 'sheba_training.csv' CSV HEADER;

-- SQLite
.mode csv
.output sheba_training.csv
SELECT * FROM ML_TrainingSet WHERE CampaignName='SHEBA';
.output stdout
```

### Python Integration

```python
import pandas as pd
from sqlalchemy import create_engine

# Connect
engine = create_engine('postgresql://user:pass@localhost/abl_database')

# Load ML features
df = pd.read_sql_query("""
    SELECT * FROM ML_TrainingSet 
    WHERE BiasRatio_B > 1.1 AND QualityFlag = 0
""", engine)

# Train model
from sklearn.ensemble import RandomForestRegressor
X = df[['Ri_Bulk_Rib', 'Zeta_NondimHeight', 'LayerThickness_Dz']]
y = df['BiasRatio_B']
model = RandomForestRegressor().fit(X, y)
```

---

## Validation Queries

### Check Data Quality

```sql
-- Sites with most observations
SELECT SiteID, CampaignName, COUNT(*) AS N_Obs
FROM Layers
GROUP BY SiteID, CampaignName
ORDER BY N_Obs DESC;

-- Stable cases by site
SELECT SiteID, COUNT(*) AS N_Stable
FROM Layers
WHERE ObukhovLength_L > 0 AND QualityFlag = 0
GROUP BY SiteID;

-- High-bias cases needing correction
SELECT SiteID, COUNT(*) AS N_HighBias
FROM Layers l
JOIN Diagnostics d ON l.ObservationID = d.ObservationID
WHERE d.BiasRatio_B > 1.5
GROUP BY SiteID;
```

### Time Series Extraction

```sql
-- ARM NSA stable episode (Jan 14-16, 2020)
SELECT TimestampUTC, LayerIndex, 
       Temperature_K_Z0, Ri_Bulk_Rib, BiasRatio_B, GridDampingFactor_G
FROM Time_Series_Validation
WHERE SiteID = 'ARM_NSA_C1'
  AND TimestampUTC BETWEEN '2020-01-14' AND '2020-01-16'
ORDER BY TimestampUTC, LayerIndex;
```

---

## Maintenance

### Update Statistics (PostgreSQL)

```sql
ANALYZE MetaData;
ANALYZE Layers;
ANALYZE Diagnostics;
```

### Backup

```bash
# PostgreSQL
pg_dump abl_database > abl_backup_$(date +%Y%m%d).sql

# SQLite
cp abl.db abl_backup_$(date +%Y%m%d).db
```

### Add New Site

```sql
INSERT INTO MetaData (SiteID, CampaignName, ...) VALUES (...);
-- Then add Layers and Diagnostics as above
```

---

## Real-World Data Sources

| Dataset | Format | Location | Notes |
|---------|--------|----------|-------|
| ARM NSA | NetCDF | https://www.arm.gov/data | Request via ARM Data Discovery |
| SHEBA | ASCII | NCAR EOL | Tower + flux data 1997-1998 |
| GABLS LES | NetCDF | KNMI server | Benchmark cases 1-4 |
| Dallas Tower | CSV | Biazar (UAH) | Urban validation, contact APB |

**Conversion scripts** in `scripts/` handle NetCDF → CSV → SQL pipeline.

---

## Contact

- **Schema questions:** david.england@uah.edu
- **Data access:** McNider/Biazar for ARM/SHEBA protocols
- **Issues:** https://github.com/DavidEngland/ABL/issues
