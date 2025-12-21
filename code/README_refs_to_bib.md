# References to BibTeX Converter

Extract references from multiple file formats and convert to BibTeX bibliography.

## Installation

```bash
# Recommended: Install PyMuPDF for best PDF extraction
pip install pymupdf

# Alternative: PyPDF2 (fallback)
pip install PyPDF2
```

## Basic Usage

### Extract from PDFs in refs/ directory

```bash
python refs_to_bib.py refs/*.pdf --output references.bib
```

### Extract from all supported formats recursively

```bash
python refs_to_bib.py refs/ --recursive --formats pdf,tex,md,html --output all_refs.bib
```

### Merge with existing bibliography

```bash
python refs_to_bib.py new_papers/*.pdf --output updated.bib --merge existing.bib
```

### Deduplicate entries

```bash
python refs_to_bib.py refs/ --recursive --deduplicate --output clean_refs.bib -v
```

## Supported Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| PDF | `.pdf` | Requires `pymupdf` or `PyPDF2` |
| LaTeX | `.tex` | Extracts `\bibitem` and `\bibliography{}` |
| BibTeX | `.bib` | Parses and cleans existing entries |
| Markdown | `.md` | Extracts from References section |
| HTML | `.html`, `.htm` | Strips tags, finds References |

## Advanced Options

```bash
# Only process specific formats
python refs_to_bib.py papers/ -f pdf,tex -o subset.bib

# Verbose output with statistics
python refs_to_bib.py refs/*.pdf -o refs.bib -v

# Recursive directory scan
python refs_to_bib.py . -r --formats pdf,md -o all.bib
```

## Output Format

Generated BibTeX entries follow standard format:

```bibtex
@article{businger1971,
  author = {Businger, J. A. and Wyngaard, J. C. and Izumi, Y. and Bradley, E. F.},
  title = {Flux-Profile Relationships in the Atmospheric Surface Layer},
  year = {1971},
  journal = {Journal of the Atmospheric Sciences},
  volume = {28},
  pages = {181--189},
  doi = {10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2},
}
```

## Citation Key Generation

Default format: `firstauthor_year` (e.g., `businger1971`)

Duplicate keys are auto-numbered: `businger1971_2`, `businger1971_3`, etc.

## Deduplication

Uses DOI as primary key (if available), otherwise normalized title.

```bash
# Remove duplicates across multiple sources
python refs_to_bib.py paper.tex paper.pdf refs/*.bib --deduplicate -o unique.bib
```

## Limitations

- **PDF extraction**: Quality depends on PDF structure. Scanned PDFs may fail.
- **Author parsing**: Heuristic-based; may miss complex name formats.
- **Title detection**: Works best with quoted or italicized titles.
- **Journal names**: Pattern-based; may not catch all variations.

## Workflow for ABL Project

```bash
# Step 1: Extract from all PDFs in refs/
python refs_to_bib.py refs/*.pdf --output refs_extracted.bib -v

# Step 2: Merge with existing BibTeX
python refs_to_bib.py refs_extracted.bib \
    drafts/*.tex \
    --merge existing_refs.bib \
    --deduplicate \
    --output ABL_bibliography.bib

# Step 3: Verify output
grep -c "^@" ABL_bibliography.bib  # Count entries
```

## Integration with LaTeX

```tex
% In your .tex file:
\bibliography{ABL_bibliography}
\bibliographystyle{plain}
```

## Troubleshooting

**No references extracted from PDF:**
- Check if References section is clearly marked
- Try `--verbose` to see extraction process
- PDF may be scanned (OCR required)

**Duplicate entries:**
- Use `--deduplicate` flag
- Manually review output for false positives

**Missing DOIs:**
- Tool extracts DOIs when present in text
- Manual addition may be needed for older papers

## Example Output Statistics

```
Found 45 files to process
Processing: businger1971.pdf...
  ✓ success: 234 chars
Processing: hogstrom1988.pdf...
  ✓ success: 189 chars
...
Deduplicated: 127 → 89 references

✓ Wrote 89 references to ABL_bibliography.bib

Entry type breakdown:
  article: 67
  inproceedings: 15
  book: 7
```
