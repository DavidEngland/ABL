"""
Runtime ML inference functions for atmospheric models.

Lightweight, vectorized operations suitable for operational NWP integration.
"""

import numpy as np
import pickle

def load_model(model_path):
    """Load pre-trained ML model."""
    with open(model_path, 'rb') as f:
        return pickle.load(f)

def predict_bias_ratio(model, delta_z, z_bottom, z0, Ri_b_raw, 
                       zeta_geom, dtheta_dz, du_dz):
    """
    Predict Richardson number bias ratio B using trained ML model.
    
    Parameters
    ----------
    model : sklearn model
        Pre-trained Random Forest or similar
    delta_z : float or array
        Layer thickness (m)
    z_bottom : float or array
        Bottom height of layer (m)
    z0 : float or array
        Surface roughness length (m)
    Ri_b_raw : float or array
        Bulk Richardson number (raw model calculation)
    zeta_geom : float or array
        Stability parameter ζ at geometric mean height
    dtheta_dz : float or array
        Temperature gradient (K/m)
    du_dz : float or array
        Wind shear (1/s)
    
    Returns
    -------
    B_pred : float or array
        Predicted bias ratio (clipped to physical bounds)
    """
    # Ensure arrays
    is_scalar = np.isscalar(delta_z)
    
    features = np.column_stack([
        delta_z, z_bottom, z0, Ri_b_raw, 
        zeta_geom, dtheta_dz, du_dz
    ])
    
    # Predict
    B_pred = model.predict(features)
    
    # Physical bounds: B must be in [1.0, 3.0]
    B_pred = np.clip(B_pred, 1.0, 3.0)
    
    return float(B_pred) if is_scalar else B_pred

def should_apply_correction_ml(clf, B_pred, zeta, Ri_b, delta_z, 
                               threshold=0.7):
    """
    Classify whether a layer needs curvature correction.
    
    Parameters
    ----------
    clf : sklearn classifier
        Pre-trained logistic regression or similar
    B_pred : float or array
        Predicted bias ratio
    zeta : float or array
        Stability parameter
    Ri_b : float or array
        Bulk Richardson number
    delta_z : float or array
        Layer thickness
    threshold : float
        Probability threshold for positive class (default 0.7)
    
    Returns
    -------
    flags : bool or array of bool
        True if correction should be applied
    """
    features = np.column_stack([B_pred, zeta, Ri_b, delta_z])
    probs = clf.predict_proba(features)[:, 1]  # Probability of class 1
    return probs > threshold

def compute_mcnider_correction(B, zeta, delta_z, alpha=1.0, gamma=1.5,
                                delta_z_ref=10.0, zeta_ref=0.5, q=2.0):
    """
    Compute McNider curvature correction factor f_c.
    
    Parameters
    ----------
    B : float or array
        Bias ratio (from ML or analytic)
    zeta : float or array
        Stability parameter
    delta_z : float or array
        Layer thickness
    alpha : float
        Correction strength parameter
    gamma : float
        Dynamic Ri_c modifier
    delta_z_ref : float
        Reference grid spacing (m)
    zeta_ref : float
        Reference stability
    q : float
        Stability exponent (default 2.0 for neutral preservation)
    
    Returns
    -------
    f_c : float or array
        Correction factor (0 < f_c <= 1)
    """
    # Exponent (neutral-preserving if q >= 2)
    exponent = -alpha * (delta_z / delta_z_ref) * (zeta / zeta_ref)**q
    
    # McNider factor
    f_c = np.exp(exponent)
    
    # Optional: dynamic Ri_c contribution
    # (simplified; full version would use inversion strength, TKE memory)
    if gamma != 0:
        Ric_mod = 1.0 + gamma * (B - 1.0)
        f_c *= Ric_mod
    
    # Bounds
    f_c = np.clip(f_c, 0.1, 1.0)
    
    return f_c

# Vectorized operational workflow
def apply_ml_correction_full(model_B, clf_layer, 
                              delta_z, z_bottom, z0, Ri_b_raw,
                              zeta_geom, dtheta_dz, du_dz,
                              K_original, alpha=1.0, gamma=1.5):
    """
    Full ML-augmented correction pipeline (vectorized for vertical column).
    
    Parameters
    ----------
    model_B : sklearn model
        Trained bias predictor
    clf_layer : sklearn classifier
        Trained layer classifier
    delta_z : array (nz,)
        Layer thicknesses
    z_bottom : array (nz,)
        Bottom heights
    ... (other arrays of same shape)
    K_original : array (nz,)
        Original eddy diffusivities
    
    Returns
    -------
    K_corrected : array (nz,)
        Corrected diffusivities
    diagnostics : dict
        Diagnostic outputs (B_pred, flags, f_c values)
    """
    nz = len(delta_z)
    
    # Step 1: Predict bias ratio (fast vectorized)
    B_pred = predict_bias_ratio(
        model_B, delta_z, z_bottom, z0, Ri_b_raw,
        zeta_geom, dtheta_dz, du_dz
    )
    
    # Step 2: Classify layers
    flags = should_apply_correction_ml(
        clf_layer, B_pred, zeta_geom, Ri_b_raw, delta_z
    )
    
    # Step 3: Compute correction factors (only where flagged)
    f_c = np.ones(nz)  # Default: no correction
    f_c[flags] = compute_mcnider_correction(
        B_pred[flags], zeta_geom[flags], delta_z[flags],
        alpha=alpha, gamma=gamma
    )
    
    # Step 4: Apply correction
    K_corrected = K_original * f_c
    
    diagnostics = {
        'B_pred': B_pred,
        'flags': flags,
        'f_c': f_c,
        'n_corrected': flags.sum()
    }
    
    return K_corrected, diagnostics
