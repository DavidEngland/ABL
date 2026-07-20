# Code Folder Organization Plan

This note tracks migration from mixed `code/` scripts to clearer folders.

## Current Direction

- New reference workflow lives in `tools/references/`.
- Existing `code/` reference scripts are preserved as legacy entry points during transition.

## Reference Workflow (New)

Use:

```bash
python tools/references/refs.py extract ...
python tools/references/refs.py verify-doi ...
python tools/references/refs.py improve ...
python tools/references/refs.py sync-paper ...
```

## Next Planned Cleanup

1. Add compatibility wrappers for old reference scripts.
2. Separate physics prototypes and legacy Fortran from reference tooling.
3. Add final migration map before any removals.
