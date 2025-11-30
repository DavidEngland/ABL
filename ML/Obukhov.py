"""
Obukhov.py — synthetic Obukhov-length experiment and tiny ML demo

Purpose
- Generate simple synthetic profiles where the Obukhov length L varies with height.
- Demonstrate via small regression tasks why the constant-L (MOST) assumption can break.
- Provide a reusable API for generating data, training quick regressors, and saving results.

Usage (CLI)
- python Obukhov.py --n 500 --out sample_obukhov.csv
- Import functions in notebooks:
    from Obukhov import generate_obukhov_profiles, train_and_evaluate

Notes
- Intended as pedagogical demo; not a production ML pipeline.
- Requires: numpy, pandas, scikit-learn
"""

from __future__ import annotations
import argparse
import logging
from typing import Dict, Tuple, Optional

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import train_test_split

# Simple logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def generate_obukhov_profiles(
    n_obs: int = 500,
    seed: Optional[int] = 42,
    layer_scale_factors: Tuple[float, float, float] = (1.0, 1.3, 1.8),
    ri_b_noise: float = 0.005,
    ri_b_clip: Tuple[float, float] = (0.05, 0.5),
) -> pd.DataFrame:
    """
    Generate synthetic soundings with layer-varying Obukhov length L(z).

    Returns a DataFrame with columns:
      - Ri_b : bulk Richardson number proxy (feature)
      - L_z1_Target, L_z2_Target, L_z3_Target : true L at layers

    Parameters
    - n_obs: number of soundings
    - seed: RNG seed for reproducibility
    - layer_scale_factors: multiplicative factors for L at z1,z2,z3 relative to surface L
    - ri_b_noise: std dev of additive noise on Ri_b
    - ri_b_clip: min/max clip for Ri_b
    """
    rng = np.random.default_rng(seed)
    # surface Obukhov length (positive = stable)
    L_surface = rng.uniform(5.0, 150.0, size=n_obs)

    rows = []
    for i in range(n_obs):
        L_s = L_surface[i]
        L_z1 = L_s * layer_scale_factors[0]
        L_z2 = L_s * layer_scale_factors[1]
        L_z3 = L_s * layer_scale_factors[2]

        # simple proxy: Ri_b correlated with 1/L_surface (small noise)
        ri_b = 0.5 * (1.0 / L_s) + rng.normal(0.0, ri_b_noise)
        ri_b = float(np.clip(ri_b, ri_b_clip[0], ri_b_clip[1]))

        rows.append(
            {
                "Ri_b": ri_b,
                "L_z1_Target": float(L_z1),
                "L_z2_Target": float(L_z2),
                "L_z3_Target": float(L_z3),
            }
        )

    df = pd.DataFrame(rows)
    return df


def train_and_evaluate(
    df: pd.DataFrame,
    feature_col: str = "Ri_b",
    test_size: float = 0.30,
    random_state: int = 42,
) -> Dict[str, Dict[str, float]]:
    """
    Train a linear regressor for each layer target and report metrics.

    Returns:
      metrics: dict keyed by target column -> {r2, rmse, rmse_rel}
    """
    X = df[[feature_col]].values
    metrics = {}
    # We'll use same train/test split for comparability across targets
    X_train, X_test, _, _ = train_test_split(X, X, test_size=test_size, random_state=random_state)
    # But we need corresponding y splits per target; do deterministic splitting using indices
    n = len(df)
    train_idx, test_idx = train_test_split(np.arange(n), test_size=test_size, random_state=random_state)

    for target in [c for c in df.columns if c != feature_col]:
        y = df[target].values
        y_train = y[train_idx]
        y_test = y[test_idx]
        X_tr = X[train_idx]
        X_te = X[test_idx]

        model = LinearRegression()
        model.fit(X_tr, y_train)
        y_pred = model.predict(X_te)

        r2 = float(r2_score(y_test, y_pred))
        rmse = float(np.sqrt(mean_squared_error(y_test, y_pred)))
        mean_target = float(np.mean(y_test)) if len(y_test) > 0 else 0.0
        rmse_rel = float(rmse / mean_target) if mean_target != 0 else float("nan")

        metrics[target] = {"r2": r2, "rmse": rmse, "rmse_rel": rmse_rel}
        logger.info(
            f"Trained model -> {target}: R2={r2:.4f}, RMSE={rmse:.3f} (rel {rmse_rel:.1%})"
        )

    return metrics


def save_df_csv(df: pd.DataFrame, path: str) -> None:
    """Save DataFrame to CSV (overwrite)."""
    df.to_csv(path, index=False)
    logger.info(f"Saved synthetic dataset -> {path}")


def quick_demo(
    n_obs: int = 500,
    seed: int = 42,
    out_csv: Optional[str] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Convenience wrapper: generate data, run train/eval, optionally save CSV.
    Returns metrics dict.
    """
    logger.info("Generating synthetic Obukhov profiles...")
    df = generate_obukhov_profiles(n_obs=n_obs, seed=seed)

    if out_csv:
        save_df_csv(df, out_csv)

    logger.info("Training simple linear regressions and reporting metrics...")
    metrics = train_and_evaluate(df)
    return metrics


def main():
    parser = argparse.ArgumentParser(description="Obukhov synthetic data + tiny ML demo")
    parser.add_argument("--n", type=int, default=500, help="Number of synthetic soundings")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed")
    parser.add_argument("--out", type=str, default=None, help="Optional CSV output path")
    args = parser.parse_args()

    metrics = quick_demo(n_obs=args.n, seed=args.seed, out_csv=args.out)

    logger.info("\nSummary metrics:")
    for target, m in metrics.items():
        logger.info(f" - {target}: R2={m['r2']:.4f}, RMSE={m['rmse']:.3f}, RMSE_rel={m['rmse_rel']:.2%}")


if __name__ == "__main__":
    main()