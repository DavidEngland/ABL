#!/usr/bin/env python3
"""Quick diagnostic tool for Richardson number bias assessment.

Usage:
    python diagnose_ri_bias.py profile.csv
    python diagnose_ri_bias.py --help

Input CSV format (with header):
    z,theta,U,V
    2.0,265.0,3.0,0.0
    10.0,265.5,5.0,0.2
    ...

Output:
    - Bias ratio B = Ri_g(z_g) / Ri_b for each layer
    - Red/green flags for correction need
    - Recommended correction strategy
    - Parameter suggestions
"""

import argparse
import sys
import numpy as np
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from rct.core.ri_estimators import ri_gradient, ri_bulk, bias_ratio


def load_profile(filepath):
    """Load vertical profile from CSV (no pandas dependency)."""
    try:
        # Read CSV manually
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Parse header
        header = lines[0].strip().split(',')
        required = ['z', 'theta', 'U', 'V']
        
        # Check required columns
        for col in required:
            if col not in header:
                print(f"ERROR: Missing required column '{col}'")
                print(f"Found columns: {header}")
                sys.exit(1)
        
        # Get column indices
        indices = {col: header.index(col) for col in required}
        
        # Parse data
        data = {col: [] for col in required}
        for line in lines[1:]:
            if line.strip():  # Skip empty lines
                values = line.strip().split(',')
                for col in required:
                    data[col].append(float(values[indices[col]]))
        
        z = np.array(data['z'])
        theta = np.array(data['theta'])
        U = np.array(data['U'])
        V = np.array(data['V'])
        
        return z, theta, U, V
    
    except FileNotFoundError:
        print(f"ERROR: File not found: {filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR loading file: {e}")
        sys.exit(1)


def diagnose_layer(z1, z2, theta1, theta2, U1, U2, V1, V2):
    """Diagnose single layer for bias."""
    
    # Compute bulk Richardson number
    Ri_b = ri_bulk(theta1, theta2, U1, U2, V1, V2, z1, z2)
    
    # Compute gradient Richardson at geometric mean height
    z_g = np.sqrt(z1 * z2)
    
    # Need z0, theta0, etc. for gradient Ri (requires triplet)
    # For 2-point layer, approximate with centered differences
    dtheta_dz = (theta2 - theta1) / (z2 - z1)
    dU_dz = (U2 - U1) / (z2 - z1)
    dV_dz = (V2 - V1) / (z2 - z1)
    
    theta_mean = 0.5 * (theta1 + theta2)
    shear_sq = dU_dz**2 + dV_dz**2
    
    if shear_sq < 1e-10:
        Ri_g = 0.0
    else:
        g = 9.81
        Ri_g = (g / theta_mean) * dtheta_dz / shear_sq
    
    # Bias ratio
    B = bias_ratio(Ri_g, Ri_b)
    
    # Grid spacing
    dz = z2 - z1
    
    return {
        'z1': z1,
        'z2': z2,
        'z_g': z_g,
        'dz': dz,
        'Ri_b': Ri_b,
        'Ri_g': Ri_g,
        'B': B,
    }


def assess_correction_need(results):
    """Determine if corrections are needed based on diagnostics."""
    
    B_values = [r['B'] for r in results]
    dz_values = [r['dz'] for r in results]
    
    # Criteria for needing corrections:
    # 1. Large bias: B > 1.5 in any layer
    # 2. Moderate bias on coarse grid: B > 1.3 and dz > 20m
    # 3. Systematic bias: mean(B) > 1.2
    
    max_B = np.max(B_values)
    mean_B = np.mean(B_values)
    max_dz = np.max(dz_values)
    
    flags = []
    
    if max_B > 1.5:
        flags.append(('RED', f"Large bias detected: B_max = {max_B:.2f}"))
    
    if mean_B > 1.2:
        flags.append(('YELLOW', f"Systematic bias: B_mean = {mean_B:.2f}"))
    
    coarse_bias = [(r['B'], r['dz']) for r in results if r['B'] > 1.3 and r['dz'] > 20]
    if coarse_bias:
        flags.append(('RED', f"{len(coarse_bias)} coarse layers with B > 1.3"))
    
    if max_dz > 50:
        flags.append(('YELLOW', f"Coarse grid: Δz_max = {max_dz:.0f} m"))
    
    if not flags:
        flags.append(('GREEN', "No significant bias detected"))
    
    return flags, max_B, mean_B


def recommend_strategy(max_B, mean_B, max_dz):
    """Recommend correction strategy based on diagnostics."""
    
    if max_B < 1.2 and mean_B < 1.1:
        return "NONE", "Bias is negligible. No corrections needed."
    
    elif max_B < 1.5 and max_dz < 30:
        return "SIMPLE", (
            "Use simple_multiplicative_correction() with D=0.3-0.5.\n"
            "   Minimal code changes, sufficient for moderate bias."
        )
    
    elif max_B < 2.0 and max_dz < 50:
        return "MODERATE", (
            "Use stability_function_correction() with D=0.5.\n"
            "   Extends tail of f(Ri) for coarse grids. More accurate than simple."
        )
    
    else:
        return "ADVANCED", (
            "Use dynamic_critical_richardson() or full CorrectionODE.\n"
            "   Large bias and/or very coarse grid requires physics-based approach.\n"
            "   Consider D=0.6-0.8 and grid-aware Ri_c adjustment."
        )


def suggest_parameters(max_B, mean_B, max_dz):
    """Suggest parameter values based on diagnostics."""
    
    # D parameter (damping amplitude)
    if max_B < 1.3:
        D = 0.3
    elif max_B < 1.7:
        D = 0.5
    else:
        D = 0.7
    
    # alpha, beta for ODE (if using advanced strategy)
    if mean_B > 1.4:
        alpha = 0.4  # Stronger curvature weight
        beta = 0.6
    else:
        alpha = 0.3  # Default
        beta = 0.7
    
    # Reference grid spacing (use median of profile)
    # For production, this should be the "target" resolution
    dz_ref = 10.0  # Default
    
    return {
        'D': D,
        'alpha': alpha,
        'beta': beta,
        'dz_ref': dz_ref,
        'p': 0.5,  # Grid scaling exponent (empirical)
        'q': 1.0,  # Ri scaling exponent (empirical)
    }


def print_report(results, flags, strategy, params):
    """Print formatted diagnostic report."""
    
    print("\n" + "="*70)
    print("RICHARDSON NUMBER BIAS DIAGNOSTIC REPORT")
    print("="*70)
    
    # Layer-by-layer results
    print("\nLayer-by-Layer Analysis:")
    print("-"*70)
    print(f"{'Layer':<8} {'Height (m)':<12} {'Δz (m)':<8} {'Ri_b':<10} {'Ri_g':<10} {'Bias B':<10}")
    print("-"*70)
    
    for i, r in enumerate(results):
        flag = "🚩" if r['B'] > 1.5 else "⚠️ " if r['B'] > 1.3 else "  "
        print(f"{flag} {i:<5} {r['z1']:5.1f}-{r['z2']:5.1f}   {r['dz']:6.1f}   "
              f"{r['Ri_b']:8.4f}  {r['Ri_g']:8.4f}  {r['B']:8.4f}")
    
    # Summary statistics
    B_values = [r['B'] for r in results]
    print("\n" + "-"*70)
    print(f"Summary Statistics:")
    print(f"  Bias ratio:  mean = {np.mean(B_values):.3f}, "
          f"max = {np.max(B_values):.3f}, min = {np.min(B_values):.3f}")
    print(f"  Grid:        Δz_max = {np.max([r['dz'] for r in results]):.1f} m, "
          f"Δz_mean = {np.mean([r['dz'] for r in results]):.1f} m")
    
    # Assessment flags
    print("\n" + "="*70)
    print("ASSESSMENT:")
    print("="*70)
    for color, message in flags:
        symbol = "🔴" if color == "RED" else "🟡" if color == "YELLOW" else "🟢"
        print(f"{symbol} {message}")
    
    # Recommendation
    print("\n" + "="*70)
    print("RECOMMENDED STRATEGY:")
    print("="*70)
    strategy_name, description = strategy
    print(f"\n{strategy_name}:")
    print(f"{description}")
    
    # Parameters
    if strategy_name != "NONE":
        print("\n" + "="*70)
        print("SUGGESTED PARAMETERS:")
        print("="*70)
        for key, value in params.items():
            print(f"  {key:10} = {value}")
        
        print("\n  Start with these values, then tune based on validation.")
        print("  See MODEL_DEVELOPER_GUIDE.md § Parameter Tuning for workflow.")
    
    print("\n" + "="*70)


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose Richardson number bias in vertical profiles.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python diagnose_ri_bias.py my_profile.csv
  python diagnose_ri_bias.py wrf_output.csv --verbose

Input file format (CSV with header):
  z,theta,U,V
  2.0,265.0,3.0,0.0
  10.0,265.5,5.0,0.2
  30.0,267.0,7.0,0.4
  ...
        """
    )
    
    parser.add_argument('input_file', help='Path to CSV file with profile data')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Print detailed diagnostics')
    parser.add_argument('-o', '--output', help='Save report to file (optional)')
    
    args = parser.parse_args()
    
    # Load profile
    print(f"Loading profile from: {args.input_file}")
    z, theta, U, V = load_profile(args.input_file)
    print(f"  Found {len(z)} vertical levels")
    print(f"  Height range: {z[0]:.1f} - {z[-1]:.1f} m")
    
    # Diagnose each layer
    results = []
    for i in range(len(z) - 1):
        result = diagnose_layer(
            z[i], z[i+1],
            theta[i], theta[i+1],
            U[i], U[i+1],
            V[i], V[i+1]
        )
        results.append(result)
    
    # Assess correction need
    flags, max_B, mean_B = assess_correction_need(results)
    
    # Recommend strategy
    max_dz = np.max([r['dz'] for r in results])
    strategy = recommend_strategy(max_B, mean_B, max_dz)
    
    # Suggest parameters
    params = suggest_parameters(max_B, mean_B, max_dz)
    
    # Print report
    print_report(results, flags, strategy, params)
    
    # Save to file if requested
    if args.output:
        import sys
        original_stdout = sys.stdout
        with open(args.output, 'w') as f:
            sys.stdout = f
            print_report(results, flags, strategy, params)
        sys.stdout = original_stdout
        print(f"\nReport saved to: {args.output}")


if __name__ == "__main__":
    main()
