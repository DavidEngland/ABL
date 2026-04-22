# SCMSkeleton vs. Parameterization Scaffold

This note evaluates `julia/SCMSkeleton.jl` against `param/SCAFFOLDING.md` and
the slab-column extension in `param/core/slab_column.md`.

## Current Alignment

The Julia SCM already matched the scaffold well in a few structural ways:

1. clear separation of grid, state, forcing, closure, timestepper, and diagnostics
2. local Richardson-number closure hooks, including a curvature-aware variant
3. flexible vertical grid support for shallow stable boundary-layer problems
4. explicit surface-flux coupling point that can be upgraded without rewriting the column solver

## Grid Design Note

The original Julia SCM used a single exponential stretch function that places
extra resolution near the surface. That remains useful for LLJ and stable ABL
experiments where the dominant gradients are expected in the lowest few hundred
meters.

Biazar's criticism is valid for the more general case: a single global stretch
is not the same thing as a physics-targeted nonuniform grid. If the important
gradients are expected near known heights such as the surface inversion, PBL top,
LLJ nose, or tropopause, the better abstraction is to refine around those target
heights directly.

`SCMSkeleton.jl` now supports both approaches:

1. `grid_strategy=:stretched` for the original exponential near-surface refinement
2. `grid_strategy=:targeted` for Gaussian-weighted refinement around estimated feature heights

The targeted grid is closer to the workflow Biazar described for WRF: estimate
important heights first, then allocate levels preferentially there.

## Main Mismatches Identified

Before this update, the largest mismatches relative to the scaffold were:

1. the lower boundary was a single bulk slab rather than a layered substrate column
2. no explicit representation of snow over land or snow over ice
3. no under-ice layer option for sea-ice cases
4. surface exchange used neutral bulk coefficients only, not regime-aware MOST logic
5. no explicit displacement-height handling
6. no prognostic snow-thickness evolution

## What Was Implemented In This Patch

The Julia skeleton is now closer to the slab-column design note.

### Added substrate structure

`SCMSkeleton.jl` now includes a `SubstrateColumnState` with:

1. substrate mode (`:land`, `:ice_ocean`, `:water`)
2. diagnosed surface state (`:bare_land`, `:snow_land`, `:ice_ocean`, etc.)
3. ordered active layers with per-layer material properties, thickness, temperature, and liquid fraction
4. snow depth constrained to remain nonnegative
5. optional thin under-ice water state with temperature and salinity placeholders

### Added material support

The default material catalog now supports:

1. `:soil`
2. `:snow`
3. `:ice`
4. `:seaice`
5. `:water`
6. `:molten_liquid`

### Updated surface coupling

The lower-boundary update now:

1. computes conductive fluxes between substrate layers
2. advances the substrate column from the net surface energy input `R_n - H - LE`
3. uses the active top material for albedo, emissivity, roughness, and moisture availability
4. keeps the old `slab` and `slab_params` views synchronized for compatibility

### Updated driver interface

`run_model` and `initialize_model` now accept:

1. `surface_material`
2. `snow_depth`
3. `has_under_ice_layer`
4. `tower`
5. `grid_strategy`
6. `stretch`
7. `target_heights`, `target_widths`, `target_strengths`

This makes snow-land and snow-on-sea-ice experiments possible without editing the model internals.
It also makes it possible to compare LLJ-style near-surface stretching against
height-targeted refinement without changing the solver core.

## Remaining Scientific Gaps

The Julia model still falls short of the full scaffold in a few important ways:

1. air-surface exchange remains neutral-bulk rather than regime-aware MOST
2. there is no iterative `zeta` or Obukhov-length solve
3. displacement height `d` is not yet carried through profiles and fluxes
4. the substrate phase-change treatment is still approximate rather than full enthalpy inversion
5. snow depth is initialized and constrained, but not yet forced by snowfall, melt, sublimation, or compaction tendencies
6. the under-ice salinity variable is present but not yet dynamically forced by freezing or mixing physics

## Recommended Next Julia Steps

The next highest-value upgrades are:

1. replace `bulk_transfer_coefficients` with a regime-aware surface-layer module using the scaffold's exact log-ratio and stability logic
2. add prognostic snow-thickness forcing terms to the `Forcing` struct and lower-boundary update
3. convert the mixed-phase substrate update from apparent liquid-fraction logic to enthalpy-based phase change
4. compute surface Richardson and transfer coefficients using tower heights plus optional displacement height
5. expose substrate diagnostics in the simulation history for snow depth, ground flux, and under-ice state
6. add adaptive target-height updates so the targeted grid can follow a diagnosed inversion, LLJ nose, or PBL top over time

## Practical Summary

The Julia SCM now has a usable lower-boundary column for snow, soil, ice, and
sea-ice prototyping. The next major physics step is not more substrate detail;
it is replacing the neutral surface exchange with the grid-aware MOST framework
defined in `param/SCAFFOLDING.md`.
