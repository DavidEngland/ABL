"""
Train a Random Forest model to predict the Richardson number bias ratio B.

Usage:
    python train_bias_predictor.py --n_samples 10000 --output models/bias_predictor.pkl
"""

import numpy as np
import pandas as pd
from scipy.integrate import quad
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
import pickle
import argparse
import matplotlib.pyplot as plt

def generate_training_data(n_samples=10000, seed=42):
    """Generate synthetic MOST profiles with known bias ratio B."""
    np.random.seed(seed)
    data = []
    
    # Parameter ranges (stable conditions)
    delta_z_range = np.logspace(0.5, 2.5, 50)  # 3-300 m
    L_range = np.logspace(0.5, 2.5, 50)        # 3-300 m (Obukhov length)
    z0_range = np.logspace(-3, 0, 20)          # 0.001-1 m (roughness)
    
    # MOST parameters (linear stable: Businger-Dyer SBL)
    a_m, a_h = 4.7, 7.8
    Delta = a_h - 2*a_m  # Neutral curvature invariant
    
    print(f"Generating {n_samples} synthetic profiles...")
    print(f"MOST parameters: a_m={a_m}, a_h={a_h}, Δ={Delta:.2f}")
    
    for i in range(n_samples):
        if i % 1000 == 0:
            print(f"  Progress: {i}/{n_samples}")
        
        z0 = np.random.choice(z0_range)
        z_bottom = np.random.uniform(z0 + 1, 50)
        delta_z = np.random.choice(delta_z_range)
        z_top = z_bottom + delta_z
        z_geom = np.sqrt(z_bottom * z_top)
        
        L = np.random.choice(L_range)
        zeta_geom = z_geom / L
        
        # Skip if outside valid MOST range
        if zeta_geom > 2.0 or zeta_geom < 0.01:
            continue
        
        # Analytic Ri_g at geometric mean height
        phi_m = 1 + a_m * zeta_geom
        phi_h = 1 + a_h * zeta_geom
        Ri_g_geom = zeta_geom * phi_h / phi_m**2
        
        # Analytic Ri_b (layer-averaged)
        def integrand(z):
            zeta = z / L
            pm = 1 + a_m * zeta
            ph = 1 + a_h * zeta
            return zeta * ph / pm**2
        
        try:
            Ri_b_true, _ = quad(integrand, z_bottom, z_top, limit=100)
            Ri_b_true /= delta_z
        except:
            continue  # Integration failed
        
        if Ri_b_true < 1e-6:
            continue
        
        # Bias ratio (TARGET)
        B_true = Ri_g_geom / Ri_b_true
        
        # Synthetic bulk quantities (what a coarse model would see)
        theta_bottom = 280.0
        theta_top = theta_bottom + 5.0 * (zeta_geom**0.5)
        dtheta = theta_top - theta_bottom
        
        U_bottom = 5.0
        u_star = 0.3
        kappa = 0.4
        U_top = U_bottom + (u_star / kappa) * np.log(z_top / z_bottom)
        dU = U_top - U_bottom
        
        g = 9.81
        theta_ref = 0.5 * (theta_bottom + theta_top)
        Ri_b_raw = (g / theta_ref) * (dtheta * delta_z) / (dU**2 + 1e-10)
        
        dtheta_dz = dtheta / delta_z
        du_dz = dU / delta_z
        
        data.append({
            'delta_z': delta_z,
            'z_bottom': z_bottom,
            'z0': z0,
            'Ri_b_raw': Ri_b_raw,
            'zeta_geom': zeta_geom,
            'dtheta_dz': dtheta_dz,
            'du_dz': du_dz,
            'L': L,
            'Ri_g_geom': Ri_g_geom,
            'Ri_b_true': Ri_b_true,
            'B_true': B_true
        })
    
    df = pd.DataFrame(data)
    print(f"Generated {len(df)} valid samples")
    print(f"B_true range: [{df['B_true'].min():.3f}, {df['B_true'].max():.3f}]")
    return df

def train_model(df, test_size=0.2, random_state=42):
    """Train Random Forest regressor."""
    feature_cols = ['delta_z', 'z_bottom', 'z0', 'Ri_b_raw', 
                    'zeta_geom', 'dtheta_dz', 'du_dz']
    target_col = 'B_true'
    
    X = df[feature_cols]
    y = df[target_col]
    
    # Split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state
    )
    
    print(f"\nTraining Random Forest on {len(X_train)} samples...")
    model = RandomForestRegressor(
        n_estimators=100,
        max_depth=12,
        min_samples_leaf=5,
        random_state=random_state,
        n_jobs=-1  # Use all CPU cores
    )
    model.fit(X_train, y_train)
    
    # Validation
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)
    
    print("\n=== Training Set ===")
    print(f"MAE:  {mean_absolute_error(y_train, y_pred_train):.4f}")
    print(f"RMSE: {np.sqrt(mean_squared_error(y_train, y_pred_train)):.4f}")
    print(f"R²:   {r2_score(y_train, y_pred_train):.4f}")
    
    print("\n=== Test Set (Generalization) ===")
    print(f"MAE:  {mean_absolute_error(y_test, y_pred_test):.4f}")
    print(f"RMSE: {np.sqrt(mean_squared_error(y_test, y_pred_test)):.4f}")
    print(f"R²:   {r2_score(y_test, y_pred_test):.4f}")
    
    # Feature importance
    importances = pd.DataFrame({
        'feature': feature_cols,
        'importance': model.feature_importances_
    }).sort_values('importance', ascending=False)
    print("\n=== Feature Importance ===")
    print(importances.to_string(index=False))
    
    return model, X_test, y_test, y_pred_test

def plot_validation(y_test, y_pred_test, output_path='validation_plot.png'):
    """Create validation scatter plot."""
    plt.figure(figsize=(8, 6))
    plt.scatter(y_test, y_pred_test, alpha=0.3, s=10)
    plt.plot([y_test.min(), y_test.max()], 
             [y_test.min(), y_test.max()], 
             'r--', lw=2, label='Perfect prediction')
    plt.xlabel('True B (from analytic MOST)', fontsize=12)
    plt.ylabel('Predicted B (ML model)', fontsize=12)
    plt.title('Bias Ratio Prediction Validation', fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"\nValidation plot saved: {output_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Train ML model for Richardson number bias prediction'
    )
    parser.add_argument('--n_samples', type=int, default=10000,
                       help='Number of synthetic profiles to generate')
    parser.add_argument('--output', type=str, default='bias_predictor.pkl',
                       help='Output path for trained model')
    parser.add_argument('--plot', type=str, default='validation_plot.png',
                       help='Output path for validation plot')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    # Generate data
    df = generate_training_data(n_samples=args.n_samples, seed=args.seed)
    
    # Train
    model, X_test, y_test, y_pred_test = train_model(df, random_state=args.seed)
    
    # Save model
    with open(args.output, 'wb') as f:
        pickle.dump(model, f)
    print(f"\nModel saved: {args.output}")
    
    # Validation plot
    plot_validation(y_test, y_pred_test, output_path=args.plot)
    
    print("\n=== Quick Usage Example ===")
    print("from ml_inference import predict_bias_ratio")
    print("B_pred = predict_bias_ratio(model, delta_z=100, z_bottom=10, ...)")
