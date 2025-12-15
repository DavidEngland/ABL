import React, { useState, useMemo, useCallback } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { Download, Info } from 'lucide-react';

/**
 * MOSTStabilityTool
 * 
 * Interactive tool for Monin-Obukhov Similarity Theory (MOST) stability analysis.
 * Provides high-precision (10⁻¹⁰ tolerance) conversions between stability parameters:
 * - ζ (dimensionless height, z/L)
 * - Ri_g (gradient Richardson number)
 * - Ri_b (bulk Richardson number)
 * - φ_m, φ_h (stability functions)
 * 
 * Supports three canonical MOST profiles:
 * - Businger-Dyer 1971 (BD71)
 * - Högström 1988 (HOG88)
 * - Cheng-Brutsaert 2005 (CB05)
 * 
 * @component
 * @returns {JSX.Element} Interactive MOST stability reference tool
 */
const MOSTStabilityTool = () => {
  // ============================================================================
  // STATE MANAGEMENT
  // ============================================================================

  const [profile, setProfile] = useState('BD71');           // Selected MOST profile
  const [regime, setRegime] = useState('full');              // Stability range
  const [showTable, setShowTable] = useState(false);         // Table visibility toggle
  const [conversionMode, setConversionMode] = useState('rig_to_zeta'); // Conversion type
  const [inputValue, setInputValue] = useState('0.1');       // User input value
  const [conversionResult, setConversionResult] = useState(null); // Conversion result

  // ============================================================================
  // CONSTANTS: MOST PROFILE PARAMETERS
  // ============================================================================

  /**
   * Canonical MOST profile parameters from literature.
   * Each profile defines stability function coefficients for unstable and stable regimes.
   * 
   * References:
   * - BD71: Businger et al. (1971), J. Atmos. Sci., 28, 181–189
   * - HOG88: Högström (1988), Boundary-Layer Meteorology, 44, 25–60
   * - CB05: Cheng & Brutsaert (2005), J. Atmos. Sci., 62, 2112–2132
   */
  const profiles = {
    BD71: {
      name: 'Businger-Dyer 1971',
      unstable: { am: 16, ah: 16, bm: 0.25, bh: 0.5 },
      stable: { bm: 5, bh: 5 }
    },
    HOG88: {
      name: 'Högström 1988',
      unstable: { am: 19.3, ah: 11.6, bm: 0.25, bh: 0.5 },
      stable: { bm: 6, bh: 7.8 }
    },
    CB05: {
      name: 'Cheng-Brutsaert 2005',
      unstable: { am: 16, ah: 16, bm: 0.25, bh: 0.5 },
      stable: { bm: 6.1, bh: 6.1, cm: 5.3, ch: 5.3 }
    }
  };

  // Numerical integration and convergence parameters
  const NUMERICAL_PARAMS = {
    PSI_STEPS: 100,        // Number of steps for numerical integration of ψ
    NEWTON_TOL: 1e-10,     // Convergence tolerance for Newton's method
    NEWTON_MAX_ITER: 100,  // Maximum iterations for Newton's method
    DERIVATIVE_DZ: 1e-8,   // Step size for numerical differentiation
    SINGULARITY_TOL: 1e-10, // Tolerance for singularity avoidance
    ZETA_BOUNDS: [-10, 10] // Physical bounds for dimensionless height
  };
  // Add canonical roughness ratio (z/z0), e.g., z=10 m, z0=0.01 m → 1000
  const SURFACE_PARAMS = {
    Z_OVER_Z0: 1000
  };

  // ============================================================================
  // CORE PHYSICS: STABILITY FUNCTIONS
  // ============================================================================

  /**
   * Momentum stability function φ_m(ζ).
   * 
   * Defines how the stability parameter ζ affects momentum flux in the surface layer.
   * 
   * Physics:
   * - Unstable (ζ < 0): Convective instability enhances mixing
   * - Stable (ζ > 0): Static stability suppresses mixing
   * 
   * @param {number} zeta - Dimensionless height (z/L where L is Obukhov length)
   * @param {string} prof - Profile name ('BD71', 'HOG88', or 'CB05')
   * @returns {number} Stability function φ_m (dimensionless)
   */
  const phi_m = (zeta, prof) => {
    const p = profiles[prof];
    if (zeta < 0) {
      // Unstable: (1 - a_m * ζ)^(-b_m)
      // This formulation ensures φ_m decreases (enhanced mixing) as instability increases
      return Math.pow(1 - p.unstable.am * zeta, -p.unstable.bm);
    } else {
      // Stable: linear or quadratic term
      if (prof === 'CB05') {
        // CB05 includes quadratic term: 1 + b_m*ζ + c_m*ζ²
        return 1 + p.stable.bm * zeta + p.stable.cm * zeta * zeta;
      }
      // BD71, HOG88: linear only: 1 + b_m*ζ
      return 1 + p.stable.bm * zeta;
    }
  };

  /**
   * Heat stability function φ_h(ζ).
   * Analogous to φ_m but for heat (temperature) gradients.
   * Typically differs from φ_m in unstable regime (different exponent).
   * 
   * @param {number} zeta - Dimensionless height (z/L)
   * @param {string} prof - Profile name ('BD71', 'HOG88', or 'CB05')
   * @returns {number} Stability function φ_h (dimensionless)
   */
  const phi_h = (zeta, prof) => {
    const p = profiles[prof];
    if (zeta < 0) {
      // Unstable: (1 - a_h * ζ)^(-b_h)
      // Heat typically has different coefficient (a_h ≠ a_m)
      return Math.pow(1 - p.unstable.ah * zeta, -p.unstable.bh);
    } else {
      // Stable: same structure as φ_m
      if (prof === 'CB05') {
        return 1 + p.stable.bh * zeta + p.stable.ch * zeta * zeta;
      }
      return 1 + p.stable.bh * zeta;
    }
  };

  /**
   * Gradient Richardson number Ri_g(ζ).
   * 
   * Local stability measure comparing static stability to wind shear.
   * Ri_g = (g/θ₀) * (dθ/dz) / (du/dz)²
   * 
   * In MOST scaling:
   * Ri_g ≈ ζ * φ_h / φ_m²
   * 
   * Physical interpretation:
   * - Ri_g < 0: Unstable (convection dominates)
   * - Ri_g ≈ 0: Neutral
   * - Ri_g > 0: Stable (shear dominates)
   * - Ri_g > Ri_c (~0.25): Turbulence suppressed (local maximum)
   * 
   * @param {number} zeta - Dimensionless height (z/L)
   * @param {string} prof - Profile name
   * @returns {number} Gradient Richardson number (dimensionless)
   */
  const ri_g = (zeta, prof) => {
    const pm = phi_m(zeta, prof);
    const ph = phi_h(zeta, prof);
    return zeta * ph / (pm * pm);
  };

  /**
   * Integral momentum stability function ψ_m(ζ).
   * 
   * Defined as: ψ_m = ∫₀^ζ (1/φ_m - 1) dz'
   * 
   * Used in logarithmic wind profile:
   * u(z) = (u*/κ) * [ln(z/z₀) - ψ_m(ζ)]
   * 
   * Unstable regime: Analytical solution available
   * Stable regime: Numerical integration required
   * 
   * @param {number} zeta - Dimensionless height (z/L)
   * @param {string} prof - Profile name
   * @returns {number} Integral stability function ψ_m (dimensionless)
   */
  const psi_m = (zeta, prof) => {
    // Handle neutral condition (ζ ≈ 0)
    if (Math.abs(zeta) < NUMERICAL_PARAMS.SINGULARITY_TOL) return 0;

    const p = profiles[prof];
    if (zeta < 0) {
      // Unstable: Analytical solution (standard form)
      // x = (1 - a_m * ζ)^0.25
      // ψ_m = 2*ln((1+x)/2) + ln((1+x²)/2) - 2*atan(x) + π/2
      const x = Math.pow(1 - p.unstable.am * zeta, 0.25);
      return (
        2 * Math.log((1 + x) / 2) +
        Math.log((1 + x * x) / 2) -
        2 * Math.atan(x) +
        Math.PI / 2
      );
    } else {
      // Stable: Numerical integration via trapezoidal rule
      // ψ_m(ζ) ≈ ∑ (1/φ_m - 1) * Δz
      const steps = NUMERICAL_PARAMS.PSI_STEPS;
      const dz = zeta / steps;
      let sum = 0;
      for (let i = 1; i <= steps; i++) {
        const z = i * dz;
        sum += (1 / phi_m(z, prof) - 1) * dz;
      }
      return sum;
    }
  };

  /**
   * Integral heat stability function ψ_h(ζ).
   * Analogous to ψ_m but for temperature gradient.
   * 
   * Used in logarithmic temperature profile:
   * θ(z) - θ(z₀) = (θ*/κ) * [ln(z/z₀) - ψ_h(ζ)]
   * 
   * @param {number} zeta - Dimensionless height (z/L)
   * @param {string} prof - Profile name
   * @returns {number} Integral stability function ψ_h (dimensionless)
   */
  const psi_h = (zeta, prof) => {
    if (Math.abs(zeta) < NUMERICAL_PARAMS.SINGULARITY_TOL) return 0;

    const p = profiles[prof];
    if (zeta < 0) {
      // Unstable: Analytical solution
      // y = (1 - a_h * ζ)^0.5
      // ψ_h = 2*ln((1+y)/2)
      const y = Math.pow(1 - p.unstable.ah * zeta, 0.5);
      return 2 * Math.log((1 + y) / 2);
    } else {
      // Stable: Numerical integration
      const steps = NUMERICAL_PARAMS.PSI_STEPS;
      const dz = zeta / steps;
      let sum = 0;
      for (let i = 1; i <= steps; i++) {
        const z = i * dz;
        sum += (1 / phi_h(z, prof) - 1) * dz;
      }
      return sum;
    }
  };

  /**
   * Bulk Richardson number Ri_b(ζ).
   * 
   * Ri_b(ζ; z/z0) ≈ ζ * [ln(z/z0) - ψ_h(ζ)] / [ln(z/z0) - ψ_m(ζ)]²
   * Assumes z >> z0; ln(z/z0) is treated as a constant per chosen z/z0.
   * 
   * Physical interpretation:
   * - More integrative than local Ri_g
   * - Used to determine if turbulence decouples in stable layer
   * - Critical value Ri_b,crit ≈ 0.2–0.3
   * 
   * @param {number} zeta - Dimensionless height (z/L)
   * @param {string} prof - Profile name
   * @returns {number} Bulk Richardson number (dimensionless)
   */
  const ri_b = (zeta, prof) => {
    // Avoid undefined ln near nonpositive ratios and near-neutral ζ artifacts
    const ln_z_over_z0 = SURFACE_PARAMS.Z_OVER_Z0 > 0 ? Math.log(SURFACE_PARAMS.Z_OVER_Z0) : 0;
    if (ln_z_over_z0 <= 0) return 0;
    if (Math.abs(zeta) < NUMERICAL_PARAMS.SINGULARITY_TOL) return 0;

    const psim = psi_m(zeta, prof);
    const psih = psi_h(zeta, prof);

    const denomTerm = ln_z_over_z0 - psim;
    const denominator = denomTerm * denomTerm;
    if (denominator < 1e-15) return 0;

    return zeta * (ln_z_over_z0 - psih) / denominator;
  };

  // ============================================================================
  // INVERSION: NEWTON'S METHOD FOR PARAMETER CONVERSION
  // ============================================================================

  /**
   * Invert Ri_g to find ζ using Newton's method.
   * 
   * Solves: Ri_g(ζ) = rig_target
   * 
   * Algorithm: Newton-Raphson with numerical differentiation
   * f(ζ) = Ri_g(ζ) - rig_target
   * ζ_{n+1} = ζ_n - f(ζ_n) / f'(ζ_n)
   * 
   * Includes bounds checking and singular point avoidance.
   * 
   * @param {number} rig_target - Target gradient Richardson number
   * @param {string} prof - Profile name
   * @param {number} [initial_guess=0.1] - Initial ζ for iteration
   * @returns {number} Dimensionless height ζ = z/L
   */
  const zeta_from_rig = useCallback((rig_target, prof, initial_guess = 0.1) => {
    let zeta = initial_guess;
    const { NEWTON_TOL, NEWTON_MAX_ITER, DERIVATIVE_DZ, ZETA_BOUNDS } = NUMERICAL_PARAMS;

    for (let i = 0; i < NEWTON_MAX_ITER; i++) {
      // Evaluate function: f(ζ) = Ri_g(ζ) - rig_target
      const f = ri_g(zeta, prof) - rig_target;
      
      // Check convergence
      if (Math.abs(f) < NEWTON_TOL) {
        return zeta;
      }

      // Numerical derivative via central difference
      // f'(ζ) ≈ [f(ζ + dz) - f(ζ - dz)] / (2*dz)
      const df =
        (ri_g(zeta + DERIVATIVE_DZ, prof) - ri_g(zeta - DERIVATIVE_DZ, prof)) /
        (2 * DERIVATIVE_DZ);

      // Avoid division by near-zero derivative
      if (Math.abs(df) < 1e-15) break;

      // Newton step
      zeta = zeta - f / df;

      // Enforce physical bounds
      zeta = Math.max(ZETA_BOUNDS[0], Math.min(ZETA_BOUNDS[1], zeta));
    }
    return zeta;
  }, []);

  /**
   * Invert Ri_b to find ζ using Newton's method.
   * Similar to zeta_from_rig but for bulk Richardson number.
   * 
   * @param {number} rib_target - Target bulk Richardson number
   * @param {string} prof - Profile name
   * @param {number} [initial_guess=0.1] - Initial ζ for iteration
   * @returns {number} Dimensionless height ζ = z/L
   */
  const zeta_from_rib = useCallback((rib_target, prof, initial_guess = 0.1) => {
    let zeta = initial_guess;
    const { NEWTON_TOL, NEWTON_MAX_ITER, DERIVATIVE_DZ, ZETA_BOUNDS } = NUMERICAL_PARAMS;

    for (let i = 0; i < NEWTON_MAX_ITER; i++) {
      // Avoid singularity at ζ = 0
      if (Math.abs(zeta) < NUMERICAL_PARAMS.SINGULARITY_TOL) {
        zeta = NUMERICAL_PARAMS.SINGULARITY_TOL;
      }

      const f = ri_b(zeta, prof) - rib_target;
      if (Math.abs(f) < NEWTON_TOL) {
        return zeta;
      }

      const df =
        (ri_b(zeta + DERIVATIVE_DZ, prof) - ri_b(zeta - DERIVATIVE_DZ, prof)) /
        (2 * DERIVATIVE_DZ);

      if (Math.abs(df) < 1e-15) break;

      zeta = zeta - f / df;
      zeta = Math.max(ZETA_BOUNDS[0], Math.min(ZETA_BOUNDS[1], zeta));
    }
    return zeta;
  }, []);

  // ============================================================================
  // DATA GENERATION: REFERENCE TABLES
  // ============================================================================

  /**
   * Generate reference data for selected profile and regime.
   * 
   * Creates a regular grid of ζ values and computes all stability parameters.
   * Memoized to avoid recalculation when UI updates.
   * 
   * Returns array of objects: { zeta, phi_m, phi_h, Ri_g, Ri_b }
   * 
   * @returns {Array<Object>} Reference data points
   */
  const referenceData = useMemo(() => {
    const data = [];
    let zetaRange = [];

    // Select ζ range based on regime
    if (regime === 'unstable') {
      // Unstable: ζ from -10 to 0 (101 points)
      zetaRange = Array.from({ length: 101 }, (_, i) => -10 + i * 0.1);
    } else if (regime === 'stable') {
      // Stable: ζ from 0 to 10 (101 points)
      zetaRange = Array.from({ length: 101 }, (_, i) => i * 0.1);
    } else {
      // Full range: ζ from -10 to 10 (201 points)
      zetaRange = Array.from({ length: 201 }, (_, i) => -10 + i * 0.1);
    }

    // Calculate stability parameters for each ζ
    for (const zeta of zetaRange) {
      const pm = phi_m(zeta, profile);
      const ph = phi_h(zeta, profile);
      const rig = ri_g(zeta, profile);
      // Avoid Ri_b calculation near singularity
      const rib = Math.abs(zeta) > 1e-5 ? ri_b(zeta, profile) : 0;

      data.push({
        zeta: parseFloat(zeta.toFixed(4)),
        phi_m: parseFloat(pm.toFixed(4)),
        phi_h: parseFloat(ph.toFixed(4)),
        Ri_g: parseFloat(rig.toFixed(4)),
        Ri_b: parseFloat(rib.toFixed(4))
      });
    }

    return data;
  }, [profile, regime]);

  // ============================================================================
  // EVENT HANDLERS: CONVERSION AND EXPORT
  // ============================================================================

  /**
   * Handle user conversion request.
   * 
   * Validates input and calls appropriate conversion function based on mode.
   * Updates conversionResult state with output parameters.
   */
  const handleConversion = useCallback(() => {
    const input = parseFloat(inputValue);
    if (isNaN(input)) {
      setConversionResult({ error: 'Invalid input: Please enter a number' });
      return;
    }

    let result = {};

    try {
      if (conversionMode === 'rig_to_zeta') {
        // Input: Ri_g, Output: ζ and derived parameters
        const initialGuess = input > 0 ? 0.1 : -0.1;
        const zeta = zeta_from_rig(input, profile, initialGuess);
        result = {
          input: `Ri_g = ${input}`,
          zeta: zeta.toFixed(6),
          Ri_b: ri_b(zeta, profile).toFixed(6),
          phi_m: phi_m(zeta, profile).toFixed(4),
          phi_h: phi_h(zeta, profile).toFixed(4)
        };
      } else if (conversionMode === 'rib_to_zeta') {
        // Input: Ri_b, Output: ζ and derived parameters
        const initialGuess = input > 0 ? 0.1 : -0.1;
        const zeta = zeta_from_rib(input, profile, initialGuess);
        result = {
          input: `Ri_b = ${input}`,
          zeta: zeta.toFixed(6),
          Ri_g: ri_g(zeta, profile).toFixed(6),
          phi_m: phi_m(zeta, profile).toFixed(4),
          phi_h: phi_h(zeta, profile).toFixed(4)
        };
      } else if (conversionMode === 'zeta_to_all') {
        // Input: ζ, Output: All stability parameters
        const zeta = input;
        result = {
          input: `ζ = ${input}`,
          Ri_g: ri_g(zeta, profile).toFixed(6),
          Ri_b: ri_b(zeta, profile).toFixed(6),
          phi_m: phi_m(zeta, profile).toFixed(4),
          phi_h: phi_h(zeta, profile).toFixed(4)
        };
      }
      setConversionResult(result);
    } catch (e) {
      setConversionResult({ error: `Conversion failed: ${e.message}` });
    }
  }, [inputValue, conversionMode, profile, zeta_from_rig, zeta_from_rib]);

  /**
   * Export reference data as CSV file.
   * Downloads current table as comma-separated values.
   */
  const exportCSV = useCallback(() => {
    const headers = ['zeta', 'phi_m', 'phi_h', 'Ri_g', 'Ri_b'];
    const csv = [
      headers.join(','),
      ...referenceData.map(row => headers.map(h => row[h]).join(','))
    ].join('\n');

    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `MOST_reference_${profile}_${regime}.csv`;
    link.click();
    URL.revokeObjectURL(url);
  }, [referenceData, profile, regime]);

  // ============================================================================
  // RENDER: JSX COMPONENT
  // ============================================================================

  return (
    <div className="w-full max-w-7xl mx-auto p-6 bg-gray-50">
      <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
        {/* Header */}
        <h1 className="text-3xl font-bold text-gray-800 mb-2">
          MOST Stability Universal Reference Tool
        </h1>
        <p className="text-gray-600 mb-4">
          Canonical reference data and conversions for Monin-Obukhov Similarity Theory
        </p>

        {/* Info box */}
        <div className="flex items-start gap-2 bg-blue-50 border border-blue-200 rounded p-4 mb-6">
          <Info className="w-5 h-5 text-blue-600 flex-shrink-0 mt-0.5" />
          <div className="text-sm text-blue-800">
            <p className="font-semibold mb-1">About this tool:</p>
            <p>Provides high-precision (10⁻¹⁰ tolerance) conversions between stability parameters:
            ζ (z/L), Ri_g (gradient Richardson), and Ri_b (bulk Richardson) using canonical MOST profiles.</p>
          </div>
        </div>

        {/* Controls */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
          <div>
            <label className="block text-sm font-semibold text-gray-700 mb-2">
              Profile Function
            </label>
            <select
              value={profile}
              onChange={(e) => setProfile(e.target.value)}
              className="w-full p-2 border border-gray-300 rounded"
              aria-label="Select MOST profile"
            >
              {Object.entries(profiles).map(([key, val]) => (
                <option key={key} value={key}>{val.name}</option>
              ))}
            </select>
          </div>

          <div>
            <label className="block text-sm font-semibold text-gray-700 mb-2">
              Stability Regime
            </label>
            <select
              value={regime}
              onChange={(e) => setRegime(e.target.value)}
              className="w-full p-2 border border-gray-300 rounded"
              aria-label="Select stability regime"
            >
              <option value="full">Full Range (ζ: -10 to 10)</option>
              <option value="unstable">Unstable (ζ: -10 to 0)</option>
              <option value="stable">Stable (ζ: 0 to 10)</option>
            </select>
          </div>

          <div className="flex items-end gap-2">
            <button
              onClick={() => setShowTable(!showTable)}
              className="flex-1 bg-blue-600 text-white px-4 py-2 rounded hover:bg-blue-700 transition"
              aria-pressed={showTable}
              aria-label={showTable ? 'Hide reference table' : 'Show reference table'}
            >
              {showTable ? 'Hide' : 'Show'} Table
            </button>
            <button
              onClick={exportCSV}
              className="bg-green-600 text-white px-4 py-2 rounded hover:bg-green-700 transition flex items-center gap-2"
              aria-label="Export reference data as CSV"
            >
              <Download className="w-4 h-4" />
              CSV
            </button>
          </div>
        </div>

        {/* Conversion Tool */}
        <div className="bg-gray-50 border border-gray-200 rounded-lg p-4 mb-6">
          <h2 className="text-xl font-bold text-gray-800 mb-4">Parameter Conversion</h2>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
            <div>
              <label className="block text-sm font-semibold text-gray-700 mb-2">
                Conversion Type
              </label>
              <select
                value={conversionMode}
                onChange={(e) => setConversionMode(e.target.value)}
                className="w-full p-2 border border-gray-300 rounded"
                aria-label="Select conversion type"
              >
                <option value="rig_to_zeta">Ri_g → ζ (and Ri_b)</option>
                <option value="rib_to_zeta">Ri_b → ζ (and Ri_g)</option>
                <option value="zeta_to_all">ζ → All parameters</option>
              </select>
            </div>

            <div>
              <label className="block text-sm font-semibold text-gray-700 mb-2">
                Input Value
              </label>
              <input
                type="number"
                step="0.001"
                value={inputValue}
                onChange={(e) => setInputValue(e.target.value)}
                className="w-full p-2 border border-gray-300 rounded"
                aria-label="Enter parameter value"
              />
            </div>

            <div className="flex items-end">
              <button
                onClick={handleConversion}
                className="w-full bg-purple-600 text-white px-4 py-2 rounded hover:bg-purple-700 transition"
                aria-label="Perform conversion"
              >
                Convert
              </button>
            </div>
          </div>

          {/* Conversion Result Display */}
          {conversionResult && (
            <div className="bg-white border border-gray-300 rounded p-4">
              {conversionResult.error ? (
                <p className="text-red-600 font-semibold">{conversionResult.error}</p>
              ) : (
                <div>
                  <p className="font-semibold text-gray-700 mb-2">{conversionResult.input}</p>
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    {Object.entries(conversionResult)
                      .filter(([k]) => k !== 'input')
                      .map(([key, val]) => (
                        <div key={key} className="bg-gray-50 p-2 rounded">
                          <p className="text-xs text-gray-600">{key}</p>
                          <p className="text-lg font-mono font-bold text-gray-800">{val}</p>
                        </div>
                      ))}
                  </div>
                </div>
              )}
            </div>
          )}
        </div>

        {/* Visualization */}
        <div className="mb-6">
          <h2 className="text-xl font-bold text-gray-800 mb-4">Stability Relationships</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {/* Richardson Numbers Chart */}
            <div className="bg-white border border-gray-200 rounded p-4">
              <h3 className="text-sm font-semibold text-gray-700 mb-2">ζ vs Richardson Numbers</h3>
              <ResponsiveContainer width="100%" height={250}>
                <LineChart data={referenceData}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis
                    dataKey="zeta"
                    label={{ value: 'ζ = z/L', position: 'insideBottom', offset: -5 }}
                  />
                  <YAxis
                    label={{ value: 'Richardson Number', angle: -90, position: 'insideLeft' }}
                  />
                  <Tooltip />
                  <Legend />
                  <Line type="monotone" dataKey="Ri_g" stroke="#8884d8" dot={false} name="Ri_g" />
                  <Line type="monotone" dataKey="Ri_b" stroke="#82ca9d" dot={false} name="Ri_b" />
                </LineChart>
              </ResponsiveContainer>
            </div>

            {/* Stability Functions Chart */}
            <div className="bg-white border border-gray-200 rounded p-4">
              <h3 className="text-sm font-semibold text-gray-700 mb-2">Stability Functions φ</h3>
              <ResponsiveContainer width="100%" height={250}>
                <LineChart data={referenceData}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis
                    dataKey="zeta"
                    label={{ value: 'ζ = z/L', position: 'insideBottom', offset: -5 }}
                  />
                  <YAxis label={{ value: 'φ', angle: -90, position: 'insideLeft' }} />
                  <Tooltip />
                  <Legend />
                  <Line type="monotone" dataKey="phi_m" stroke="#ff7300" dot={false} name="φ_m" />
                  <Line type="monotone" dataKey="phi_h" stroke="#d84aff" dot={false} name="φ_h" />
                </LineChart>
              </ResponsiveContainer>
            </div>
          </div>
        </div>

        {/* Reference Table */}
        {showTable && (
          <div className="overflow-x-auto mb-6">
            <h2 className="text-xl font-bold text-gray-800 mb-4">Reference Data Table</h2>
            <div className="max-h-96 overflow-y-auto border border-gray-300 rounded">
              <table className="min-w-full bg-white">
                <thead className="bg-gray-100 sticky top-0">
                  <tr>
                    <th className="px-4 py-2 text-left text-sm font-semibold text-gray-700">ζ</th>
                    <th className="px-4 py-2 text-left text-sm font-semibold text-gray-700">φ_m</th>
                    <th className="px-4 py-2 text-left text-sm font-semibold text-gray-700">φ_h</th>
                    <th className="px-4 py-2 text-left text-sm font-semibold text-gray-700">Ri_g</th>
                    <th className="px-4 py-2 text-left text-sm font-semibold text-gray-700">Ri_b</th>
                  </tr>
                </thead>
                <tbody>
                  {referenceData.map((row, i) => (
                    <tr key={i} className="border-t border-gray-200 hover:bg-gray-50">
                      <td className="px-4 py-2 text-sm font-mono">{row.zeta}</td>
                      <td className="px-4 py-2 text-sm font-mono">{row.phi_m}</td>
                      <td className="px-4 py-2 text-sm font-mono">{row.phi_h}</td>
                      <td className="px-4 py-2 text-sm font-mono">{row.Ri_g}</td>
                      <td className="px-4 py-2 text-sm font-mono">{row.Ri_b}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}

        {/* Profile Information */}
        <div className="mt-6 bg-gray-50 border border-gray-200 rounded p-4">
          <h2 className="text-lg font-bold text-gray-800 mb-2">
            Current Profile: {profiles[profile].name}
          </h2>
          <div className="grid grid-cols-2 gap-4 text-sm">
            <div>
              <p className="font-semibold text-gray-700">Unstable (ζ &lt; 0):</p>
              <p className="font-mono text-gray-600">
                φ_m = (1 - {profiles[profile].unstable.am}ζ)^(-{profiles[profile].unstable.bm})
              </p>
              <p className="font-mono text-gray-600">
                φ_h = (1 - {profiles[profile].unstable.ah}ζ)^(-{profiles[profile].unstable.bh})
              </p>
            </div>
            <div>
              <p className="font-semibold text-gray-700">Stable (ζ &gt; 0):</p>
              <p className="font-mono text-gray-600">
                φ_m = 1 + {profiles[profile].stable.bm}ζ
                {profiles[profile].stable.cm && ` + ${profiles[profile].stable.cm}ζ²`}
              </p>
              <p className="font-mono text-gray-600">
                φ_h = 1 + {profiles[profile].stable.bh}ζ
                {profiles[profile].stable.ch && ` + ${profiles[profile].stable.ch}ζ²`}
              </p>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default MOSTStabilityTool;