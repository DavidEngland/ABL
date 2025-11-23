# Beamer slide deck with SVG graphics and speaker notes

Below are ready-to-use SVG graphics plus a LaTeX Beamer deck that embeds them and includes presenter notes. The graphics focus on SBL curvature (concave-down) and contrast with BD USL (unstable) behavior, plus grid-averaging bias and validation scaffolding.

---

## SVG graphics

### SBL vs BD USL curvature fingerprints (Ri_g vs ζ)

```svg
<?xml version="1.0" encoding="UTF-8"?>
<svg width="720" height="420" viewBox="0 0 720 420" xmlns="http://www.w3.org/2000/svg">
  <style>
    .axis { stroke: #333; stroke-width: 2; }
    .grid { stroke: #ddd; stroke-width: 1; }
    .label { fill: #333; font-size: 14px; font-family: Arial, sans-serif; }
    .title { fill: #111; font-size: 18px; font-family: Arial, sans-serif; font-weight: bold; }
    .curve-sbl { stroke: #1f77b4; stroke-width: 3; fill: none; }
    .curve-usl { stroke: #d62728; stroke-width: 3; fill: none; stroke-dasharray: 8 4; }
    .neutral { stroke: #2ca02c; stroke-width: 2; fill: none; stroke-dasharray: 4 6; }
    .legend { font-size: 13px; }
  </style>

  <!-- Title -->
  <text class="title" x="20" y="30">Curvature fingerprints: SBL vs BD USL (Ri_g vs ζ)</text>

  <!-- Axes -->
  <line class="axis" x1="60" y1="360" x2="680" y2="360"/>
  <line class="axis" x1="60" y1="360" x2="60" y2="60"/>
  <text class="label" x="340" y="390">ζ (z/L)</text>
  <text class="label" x="20" y="60">Ri_g</text>

  <!-- Grid lines -->
  <g>
    <line class="grid" x1="60" y1="300" x2="680" y2="300"/>
    <line class="grid" x1="60" y1="240" x2="680" y2="240"/>
    <line class="grid" x1="60" y1="180" x2="680" y2="180"/>
    <line class="grid" x1="60" y1="120" x2="680" y2="120"/>
    <line class="grid" x1="140" y1="360" x2="140" y2="60"/>
    <line class="grid" x1="220" y1="360" x2="220" y2="60"/>
    <line class="grid" x1="300" y1="360" x2="300" y2="60"/>
    <line class="grid" x1="380" y1="360" x2="380" y2="60"/>
    <line class="grid" x1="460" y1="360" x2="460" y2="60"/>
    <line class="grid" x1="540" y1="360" x2="540" y2="60"/>
    <line class="grid" x1="620" y1="360" x2="620" y2="60"/>
  </g>

  <!-- Neutral reference line -->
  <path class="neutral" d="M60,340 C180,320 300,300 420,280 C520,260 600,240 680,220"/>
  <text class="legend" x="500" y="245" fill="#2ca02c">Neutral reference</text>

  <!-- SBL concave-down curve (Δ < 0) -->
  <path class="curve-sbl" d="M60,340 C200,300 320,270 420,250 C520,235 600,230 680,228"/>
  <text class="legend" x="420" y="255" fill="#1f77b4">SBL (Δ &lt; 0, concave-down)</text>

  <!-- BD USL curve (unstable; Δ ≳ 0) -->
  <path class="curve-usl" d="M60,340 C160,335 260,330 360,324 C480,316 580,300 680,280"/>
  <text class="legend" x="450" y="285" fill="#d62728">BD USL (unstable; Δ ≳ 0)</text>

  <!-- Annotation -->
  <text class="label" x="80" y="80">SBL curvature is analyzed here; BD USL is shown only as contrast.</text>
</svg>
```

---

### Grid-averaging bias: Ri_b vs true Ri_g

```svg
<?xml version="1.0" encoding="UTF-8"?>
<svg width="720" height="420" viewBox="0 0 720 420" xmlns="http://www.w3.org/2000/svg">
  <style>
    .axis { stroke: #333; stroke-width: 2; }
    .label { fill: #333; font-size: 14px; font-family: Arial, sans-serif; }
    .title { fill: #111; font-size: 18px; font-family: Arial, sans-serif; font-weight: bold; }
    .curve-true { stroke: #1f77b4; stroke-width: 3; fill: none; }
    .avg-layer { stroke: #ff7f0e; stroke-width: 3; fill: none; stroke-dasharray: 6 4; }
    .point { fill: #1f77b4; }
    .point-avg { fill: #ff7f0e; }
    .brace { stroke: #999; stroke-width: 2; fill: none; }
  </style>

  <!-- Title -->
  <text class="title" x="20" y="30">Layer averaging bias: Ri_b underestimates Ri_g in SBL</text>

  <!-- Axes -->
  <line class="axis" x1="100" y1="360" x2="660" y2="360"/>
  <line class="axis" x1="100" y1="360" x2="100" y2="60"/>
  <text class="label" x="360" y="390">Height z</text>
  <text class="label" x="60" y="60">Ri</text>

  <!-- True SBL profile (concave-down) -->
  <path class="curve-true" d="M100,340 C240,280 360,240 480,220 C560,210 620,208 660,208"/>

  <!-- Layer bounds -->
  <line x1="100" y1="300" x2="660" y2="300" stroke="#bbb" stroke-dasharray="4 4"/>
  <line x1="100" y1="220" x2="660" y2="220" stroke="#bbb" stroke-dasharray="4 4"/>
  <text class="label" x="660" y="305" text-anchor="end">z₀</text>
  <text class="label" x="660" y="225" text-anchor="end">z₁</text>

  <!-- Geometric mean height marker -->
  <line x1="100" y1="260" x2="660" y2="260" stroke="#999" stroke-dasharray="2 6"/>
  <text class="label" x="660" y="265" text-anchor="end">z_g = √(z₀ z₁)</text>

  <!-- Average (bulk) Ri_b line between layer bounds -->
  <line class="avg-layer" x1="100" y1="260" x2="660" y2="260"/>
  <circle class="point-avg" cx="600" cy="260" r="4"/>
  <text class="label" x="610" y="255">Ri_b</text>

  <!-- True Ri_g at z_g -->
  <circle class="point" cx="600" cy="235" r="4"/>
  <text class="label" x="610" y="235">Ri_g(z_g)</text>

  <!-- Brace showing bias -->
  <path class="brace" d="M600,235 C590,240 590,255 600,260"/>
  <text class="label" x="560" y="250" text-anchor="end">Bias: Ri_b &lt; Ri_g</text>
</svg>
```

---

### Δ fingerprint illustration (concave-up vs concave-down)

```svg
<?xml version="1.0" encoding="UTF-8"?>
<svg width="720" height="360" viewBox="0 0 720 360" xmlns="http://www.w3.org/2000/svg">
  <style>
    .title { fill: #111; font-size: 18px; font-family: Arial, sans-serif; font-weight: bold; }
    .axis { stroke: #333; stroke-width: 2; }
    .label { fill: #333; font-size: 14px; font-family: Arial, sans-serif; }
    .up { stroke: #2ca02c; stroke-width: 3; fill: none; }
    .down { stroke: #1f77b4; stroke-width: 3; fill: none; }
    .zero { stroke: #999; stroke-width: 3; fill: none; stroke-dasharray: 6 4; }
  </style>

  <text class="title" x="20" y="30">Δ fingerprint of stability functions</text>
  <line class="axis" x1="80" y1="300" x2="640" y2="300"/>
  <line class="axis" x1="80" y1="300" x2="80" y2="60"/>
  <text class="label" x="340" y="325">ζ</text>
  <text class="label" x="40" y="60">Curvature</text>

  <path class="zero" d="M80,280 C200,275 320,270 440,265 C560,260 600,255 640,250"/>
  <text class="label" x="520" y="255">Δ = 0</text>

  <path class="up" d="M80,290 C200,300 320,305 440,310 C560,312 600,314 640,315"/>
  <text class="label" x="520" y="315">Δ &gt; 0 (concave-up)</text>

  <path class="down" d="M80,290 C200,270 320,250 440,235 C560,225 600,220 640,218"/>
  <text class="label" x="500" y="230">Δ &lt; 0 (concave-down, SBL)</text>
</svg>
```

---

### Validation layout scaffold (observations vs model)

```svg
<?xml version="1.0" encoding="UTF-8"?>
<svg width="720" height="420" viewBox="0 0 720 420" xmlns="http://www.w3.org/2000/svg">
  <style>
    .title { fill: #111; font-size: 18px; font-family: Arial, sans-serif; font-weight: bold; }
    .axis { stroke: #333; stroke-width: 2; }
    .label { fill: #333; font-size: 14px; font-family: Arial, sans-serif; }
    .obs { stroke: #000; stroke-width: 2; fill: none; }
    .unc { stroke: #d62728; stroke-width: 3; fill: none; }
    .cor { stroke: #1f77b4; stroke-width: 3; fill: none; }
    .rmse { fill: #666; font-size: 13px; }
  </style>

  <text class="title" x="20" y="30">Tower validation: obs vs uncorrected vs corrected</text>
  <line class="axis" x1="80" y1="360" x2="640" y2="360"/>
  <line class="axis" x1="80" y1="360" x2="80" y2="80"/>
  <text class="label" x="340" y="385">Time</text>
  <text class="label" x="40" y="80">Metric</text>

  <!-- Placeholder polylines -->
  <polyline class="obs" points="80,300 140,290 200,280 260,275 320,270 380,268 440,270 500,275 560,280 620,290"/>
  <polyline class="unc" points="80,320 140,315 200,310 260,308 320,306 380,305 440,306 500,310 560,315 620,318"/>
  <polyline class="cor" points="80,305 140,295 200,285 260,280 320,278 380,277 440,279 500,282 560,286 620,292"/>

  <!-- Legend -->
  <rect x="420" y="90" width="220" height="90" fill="#f8f8f8" stroke="#ccc"/>
  <line class="obs" x1="430" y1="110" x2="470" y2="110"/>
  <text class="label" x="480" y="114">Observations</text>
  <line class="unc" x1="430" y1="135" x2="470" y2="135"/>
  <text class="label" x="480" y="139">Model (uncorrected)</text>
  <line class="cor" x1="430" y1="160" x2="470" y2="160"/>
  <text class="label" x="480" y="164">Model (corrected)</text>

  <!-- RMSE placeholders -->
  <text class="rmse" x="100" y="100">RMSE uncorrected: 0.45</text>
  <text class="rmse" x="100" y="120">RMSE corrected: 0.26</text>
</svg>
```

---

## Beamer .tex file with framed slides and notes

Copy into a file like slides.tex, place the SVGs in the same directory, and compile with lualatex or xelatex for SVG support (or convert SVGs to PDF/PNG if using pdflatex).

```tex
\documentclass[aspectratio=169]{beamer}
\usetheme{metropolis}
\usepackage{graphicx}
\usepackage{svg}
\usepackage{amsmath}
\title{MOST, Stable Boundary Layer, and Curvature-aware Corrections}
\author{Your Name}
\date{\today}

\begin{document}
\maketitle

\begin{frame}{Motivation}
  \begin{itemize}
    \item \textbf{Problem:} Coarse vertical grids misclassify stability, yielding warm biases and pollutant errors.
    \item \textbf{Goal:} Resolution-aware correction that preserves neutral physics and reduces bias.
    \item \textbf{Scope:} Curvature analysis for the \textbf{stable boundary layer (SBL)}, not BD USL parameters.
  \end{itemize}
  \note{
    Set the stage: stability misclassification near the surface has outsized impacts.
    Emphasize that our curvature discussion is SBL-specific. BD USL applies to unstable daytime turbulence and is used only for contrast.
  }
\end{frame}

\begin{frame}{MOST essentials}
  \begin{itemize}
    \item \textbf{Key variable:} \(\zeta = z/L\) (dimensionless height).
    \item \textbf{Functions:} \(\phi_m(\zeta), \phi_h(\zeta)\) govern shear/heat transfer.
    \item \textbf{Related measures:} \(Ri_g\), stability via shear and buoyancy.
  \end{itemize}
  \note{
    Keep intuition: near-surface turbulence collapses to \(\zeta\).
    We'll link \(\phi_m, \phi_h\) to curvature via a compact fingerprint \(\Delta\).
  }
\end{frame}

\begin{frame}{Curvature fingerprints: SBL vs BD USL}
  \begin{center}
    \includesvg[width=0.9\linewidth]{curvature_sbl_usl.svg}
  \end{center}
  \begin{itemize}
    \item \textbf{SBL:} Concave-down Ri\(_g\) vs \(\zeta\) (negative curvature; \(\Delta < 0\)).
    \item \textbf{BD USL:} Unstable branch often nearer neutral/concave-up; different parameter regime.
    \item \textbf{Clarification:} Our correction targets SBL curvature only.
  \end{itemize}
  \note{
    Point at the blue curve: SBL concave-down behavior.
    Red dashed is BD USL—showing contrast but not used for parameterization here.
  }
\end{frame}

\begin{frame}{Neutral curvature invariant}
  \begin{itemize}
    \item \textbf{Fingerprint:} \(\Delta = \alpha_h \beta_h - 2\,\alpha_m \beta_m\).
    \item \textbf{Interpretation:} Sign of \(\Delta\) sets concavity of stability curves.
    \item \textbf{SBL fits:} Typically yield \(\Delta < 0\) (concave-down).
  \end{itemize}
  \note{
    You can motivate \(\Delta\) as a compact way to preserve neutral behavior while correcting curvature effects.
    No need to deep-derive; focus on physical meaning.
  }
\end{frame}

\begin{frame}{Layer averaging bias in SBL}
  \begin{center}
    \includesvg[width=0.9\linewidth]{grid_bias.svg}
  \end{center}
  \begin{itemize}
    \item \textbf{Issue:} Averaging across a concave-down profile yields \(Ri_b < Ri_g(z_g)\).
    \item \textbf{Geometric mean height:} \(z_g = \sqrt{z_0 z_1}\) is the natural evaluation point.
    \item \textbf{Impact:} Systematic underestimation of stability on coarse grids.
  \end{itemize}
  \note{
    Invoke Jensen’s inequality informally: concave-down -> average below the curve.
    Emphasize evaluating at \(z_g\) rather than linear midpoints.
  }
\end{frame}

\begin{frame}{Correction strategy}
  \begin{itemize}
    \item \textbf{Preserve:} Neutral behavior and the curvature fingerprint \(2\Delta\).
    \item \textbf{Adjust:} Bulk Ri\(_b\) toward the true Ri\(_g(z_g)\) with grid-aware damping.
    \item \textbf{Converge:} Correction vanishes as \(\Delta z \to 0\).
  \end{itemize}
  \note{
    Core principle: don’t break neutrality; correct only the curvature-induced bias.
    Stress resolution-awareness and graceful convergence.
  }
\end{frame}

\begin{frame}{Validation scaffold}
  \begin{center}
    \includesvg[width=0.9\linewidth]{validation_layout.svg}
  \end{center}
  \begin{itemize}
    \item \textbf{Comparison:} Observations vs model (uncorrected) vs corrected.
    \item \textbf{Metrics:} RMSE reduction indicative of improved stability depiction.
    \item \textbf{Cases:} Apply across diverse SBL nights to test robustness.
  \end{itemize}
  \note{
    Replace placeholders with your tower/time-series data.
    Narrate RMSE changes and qualitative fit improvements (phase, amplitude).
  }
\end{frame}

\begin{frame}{Limitations and open questions}
  \begin{itemize}
    \item \textbf{Transferability:} Parameter portability across terrains and urban canopies.
    \item \textbf{Inflections:} Handling mixed regimes and nocturnal jets.
    \item \textbf{L(z):} Height-varying Obukhov length effects.
  \end{itemize}
  \note{
    Keep it candid: where the approach might struggle.
    Invite collaboration on canopy and complex terrain extensions.
  }
\end{frame}

\begin{frame}{Takeaways}
  \begin{itemize}
    \item \textbf{Curvature matters:} It sets departures from neutrality.
    \item \textbf{Evaluate smartly:} Use \(z_g\) and preserve \(2\Delta\).
    \item \textbf{Correct gently:} Grid-aware adjustments reduce biases without breaking physics.
  \end{itemize}
  \note{
    Summarize and reiterate SBL-only applicability. BD USL parameters are not the source of the curvature used in this correction.
  }
\end{frame}

\end{document}
```

---

## Speaker notes script cues

- **Motivation:** Biases cluster in the lowest 100 m; stability misclassification propagates into temperature and pollutant forecasts.
- **MOST essentials:** Anchor on ζ = z/L; φ_m, φ_h dictate transport; Ri_g interprets the balance.
- **Curvature fingerprints:** SBL shows concave-down Ri_g(ζ) due to Δ < 0; BD USL is a different regime and not used here.
- **Neutral invariant:** Δ captures the curvature character while maintaining neutral behavior.
- **Layer averaging bias:** Jensen’s inequality intuition—averages fall below concave-down curves; evaluate at z_g.
- **Correction strategy:** Preserve neutrality and 2Δ; grid-aware correction fades with resolution.
- **Validation:** Show RMSE improvements and shape fidelity across cases.
- **Limitations:** Terrain, urban canopy, nocturnal jets, L(z) variability.
- **Takeaways:** Respect curvature, honor physics, SBL-specific.

---

## Practical tips

- **Compilation:** Prefer lualatex/xelatex for native SVG inclusion via the svg package; otherwise, convert SVGs to PDF using inkscape and include with includegraphics.
- **Colors/fonts:** Keep blue for SBL, red for BD USL contrast, green for neutral references for consistent visual semantics.
- **Data swap:** Replace placeholder polylines in SVGs with your data-driven paths for case studies.

If you want, I can tailor the SVG shapes to match your exact parameter values or generate additional figures (e.g., φ_m/φ_h overlays, Ri_g second derivative sketches).