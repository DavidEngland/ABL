# McNider & England (1995) Reference Finder

Find all references to your stability functions paper across multiple file formats.

## Installation

```bash
# Required for PDF support
pip install pymupdf  # Recommended (faster, better text extraction)
# OR
pip install PyPDF2   # Fallback option

# Required for DOCX support
pip install python-docx

# All in one
pip install pymupdf python-docx
```

## Basic Usage

### Search current directory and all subdirectories

```bash
python find_mcnider_refs.py
```

### Search specific directory

```bash
python find_mcnider_refs.py /Users/davidengland/Documents/GitHub/ABL
```

### Search only PDFs and Markdown files

```bash
python find_mcnider_refs.py --formats pdf,md
```

### Save report to file

```bash
python find_mcnider_refs.py --output mcnider_refs_report.txt
```

### Search without recursing into subdirectories

```bash
python find_mcnider_refs.py --no-recursive
```

### Verbose output (see all files being searched)

```bash
python find_mcnider_refs.py --verbose
```

## Search Patterns

The tool searches for multiple citation patterns:

1. **Title matches:**
   - "Stability functions based upon shear functions"
   
2. **Author-year citations:**
   - McNider & England (1995)
   - McNider and England 1995
   - England & McNider (1995)
   
3. **BibTeX keys:**
   - mcnider1995
   - england1995
   - mcnider_1995
   
4. **Journal citations:**
   - Boundary-Layer Meteorology, Vol 72-74
   
5. **Combined patterns:**
   - Any mention of "McNider" and "England" with "1995" and keywords like "stability" or "shear"

## Supported File Formats

| Format | Extension | Library Required |
|--------|-----------|------------------|
| PDF | `.pdf` | pymupdf or PyPDF2 |
| Markdown | `.md` | (built-in) |
| LaTeX | `.tex`, `.latex` | (built-in) |
| BibTeX | `.bib` | (built-in) |
| Text | `.txt` | (built-in) |
| Word | `.docx` | python-docx |
| reStructuredText | `.rst` | (built-in) |

## Output Format

```
======================================================================
McNider & England (1995) Reference Search Report
======================================================================

Total files with references: 12
Total references found: 23

Files by type:
  .bib: 1 file(s)
  .md: 7 file(s)
  .pdf: 2 file(s)
  .tex: 2 file(s)

======================================================================
Detailed Results
======================================================================

/Users/davidengland/Documents/GitHub/ABL/drafts/Curvature_Correction_Draft.md
----------------------------------------------------------------------
  Line 487: ...McNider & England (1995) sensitivity study showed...
  Line 523: ...based on the stability functions derived in...

/Users/davidengland/Documents/GitHub/ABL/ABL_refs.bib
----------------------------------------------------------------------
  Line 1042: ...title = {Stability functions based upon shear functions}...
```

## Example Workflows

### Find all references in your ABL project

```bash
cd /Users/davidengland/Documents/GitHub/ABL
python code/find_mcnider_refs.py . --output reports/mcnider_citations.txt
```

### Check only draft manuscripts

```bash
python code/find_mcnider_refs.py drafts/ --formats md,tex --verbose
```

### Search publications folder for PDFs only

```bash
python code/find_mcnider_refs.py papers/ --formats pdf --output pdf_refs.txt
```

### Quick check of BibTeX files

```bash
python code/find_mcnider_refs.py . --formats bib --no-recursive
```

## Output Interpretation

**Line numbers:** 
- For text files (MD, TEX, TXT): actual line numbers in the file
- For PDF files: approximate page breaks (may not be exact)
- For DOCX files: paragraph numbers

**Context:**
- Shows ~40 characters before and after the match
- Helps identify false positives vs. actual citations

## Common Issues

**"No PDF library found":**
```bash
pip install pymupdf
```

**"python-docx not found":**
```bash
pip install python-docx
```

**PDF text extraction poor quality:**
- Scanned PDFs may not extract well
- Try pymupdf instead of PyPDF2 (better quality)
- Consider using OCR for scanned documents

**Too many false positives:**
- Review the context column to distinguish actual citations from coincidental matches
- The tool intentionally casts a wide net to avoid missing references

## Advanced Usage

### Create custom search patterns

Edit the `PATTERNS` list in `McNiderRefFinder` class to add project-specific patterns:

```python
# Add custom patterns
PATTERNS = [
    r'stability\s+shear\s+function',  # Custom phrase
    r'England.*McNider.*shear',        # Custom author order
]
```

### Exclude specific directories

```bash
# Use find with exclude patterns
find . -type f -name "*.pdf" ! -path "*/archive/*" | xargs python find_mcnider_refs.py
```

### Integration with version control

```bash
# Find references in changed files only
git diff --name-only | grep -E '\.(md|tex|bib)$' | xargs python find_mcnider_refs.py
```

## Related Tools

- `refs_to_bib.py`: Extract references from documents and create BibTeX
- `pdf_extract.py`: Extract abstracts and metadata from PDFs

## Citation Format Reminder

Your paper:
```bibtex
@article{england1995,
  author = {England, D. E. and McNider, R. T.},
  title = {Stability functions based upon shear functions},
  journal = {Boundary-Layer Meteorology},
  year = {1995},
  volume = {74},
  pages = {113--130},
  doi = {10.1007/BF00715712}
}
```
