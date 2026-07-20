# References Toolkit (Standalone)

Standalone utilities for reference extraction, DOI verification, bibliography improvement, and paper citation sync.

## Quick Start

```bash
# From repository root
python tools/references/refs.py --help
```

## Commands

### 1) Extract references to BibTeX

```bash
python tools/references/refs.py extract \
  code/ refs/ manuscripts/ \
  --recursive \
  --formats pdf,tex,bib,md,html \
  --output-bib tools/references/output/extracted.bib
```

### 2) Verify DOI metadata against BibTeX fields

```bash
python tools/references/refs.py verify-doi \
  --bib ABL_refs.bib complete_bibliography.bib \
  --email David.England@UAH.Edu \
  --report-md tools/references/output/doi_verification_report.md \
  --report-json tools/references/output/doi_verification_report.json
```

### 3) Improve bibliography (balanced mode)

Balanced mode behavior:
- Fill missing fields from Crossref DOI metadata.
- Do not overwrite populated conflicting local fields.
- Report conflicts for manual review.

```bash
python tools/references/refs.py improve \
  --bib ABL_refs.bib \
  --email David.England@UAH.Edu \
  --output-bib tools/references/output/ABL_refs_improved.bib \
  --report-md tools/references/output/improve_report.md \
  --report-json tools/references/output/improve_report.json
```

### 4) Sync TeX citations with BibTeX keys

```bash
python tools/references/refs.py sync-paper \
  --tex "manuscripts/*.tex" \
  --bib ABL_refs.bib \
  --report-md tools/references/output/sync_paper_report.md \
  --report-json tools/references/output/sync_paper_report.json
```

## Output Artifacts

- `tools/references/output/*.bib` improved/extracted bibliography files
- `tools/references/output/*_report.md` human-readable reports
- `tools/references/output/*_report.json` machine-readable reports
- `tools/references/output/crossref_cache.json` DOI metadata cache

## Notes and Limits

- Crossref verification requires internet access.
- PDF extraction is lightweight and best-effort.
- Install `pymupdf` (preferred) or `PyPDF2` to enable PDF extraction.
- BibTeX parsing is robust for common entries and may need manual cleanup for highly nonstandard files.

## Dependency Additions

Optional packages:

```bash
pip install pymupdf PyPDF2
```

Crossref access uses Python standard library HTTP (no extra package required).
