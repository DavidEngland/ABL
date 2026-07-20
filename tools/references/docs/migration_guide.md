# Code Folder Reorganization (Phase 1)

This phase introduces a standalone references toolkit without breaking existing scripts.

## What Was Added

- `tools/references/refs.py`: new CLI entry point
- `tools/references/refkit/*`: modular implementation for extraction, DOI verification, improvement, and paper sync
- `tools/references/README.md`: usage guide

## What Was Not Changed Yet

- Existing scripts in `code/` remain in place.
- No destructive moves or path removals were performed.

## Planned Next Reorganization Steps

1. Keep `code/doi.py`, `code/refs_to_bib.py`, and `code/pdf_extract.py` as legacy paths during transition.
2. Add compatibility wrappers in `code/` that call `tools/references/refs.py` commands.
3. Move legacy Fortran and non-reference scripts to dedicated folders in later phases.
4. Publish an old-path to new-path map and deprecation timeline.

## Safety Policy

- Do not delete old scripts until the new CLI has passed real-paper verification on your `.tex` and `.bib` inputs.
- Preserve backward compatibility during at least one transition cycle.
