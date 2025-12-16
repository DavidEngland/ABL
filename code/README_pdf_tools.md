# PDF Abstract Extraction Tools

## Quick Start

### Install dependencies (in your virtual environment)

```bash
# Activate your venv first
source /Users/davidengland/Documents/GitHub/ABL/.venv/bin/activate

# Install PyMuPDF (recommended - faster, better text extraction)
pip install pymupdf

# OR install PyPDF2 (fallback)
pip install PyPDF2
```

### Usage Examples

```bash
# Extract abstract from single PDF
python code/pdf_extract.py refs/paper.pdf

# Process multiple PDFs and save to JSON
python code/pdf_extract.py refs/*.pdf --output paper_summaries.json

# Verbose mode (show progress)
python code/pdf_extract.py refs/arXiv-*/*.pdf -v --output arxiv_abstracts.json

# Process specific folder
python code/pdf_extract.py "refs/arXiv-2212.06889v3/manuscript.pdf" -v
```

### Output Format (JSON)

```json
[
  {
    "file": "manuscript.pdf",
    "title": "Paper Title Here",
    "year": "2022",
    "doi": "10.1234/example",
    "venue": "Journal Name",
    "abstract": "Full abstract text extracted from the PDF...",
    "status": "success"
  }
]
```

### Integration with your workflow

Create a batch script to process all your references:

```bash
#!/bin/bash
# Extract abstracts from all PDFs in refs/

python code/pdf_extract.py \
  refs/arXiv-*/manuscript*.pdf \
  refs/*.pdf \
  --output data/paper_abstracts.json \
  --verbose
```

### Troubleshooting

**"No PDF library found" error:**
```bash
pip install pymupdf
```

**Abstract not detected:**
- Script scans first 3 pages by default
- Increase with `--max-pages 5`
- Some PDFs have non-standard formatting (scanned images won't work without OCR)

**For scanned PDFs (OCR needed):**
```bash
pip install pytesseract
# (requires tesseract binary: brew install tesseract)
```

## Alternative: Command-line PDF tools

### Using pdftotext (fast, simple)

```bash
# Install poppler tools
brew install poppler

# Extract text from PDF
pdftotext -l 2 refs/paper.pdf - | head -50

# Extract first page to text file
pdftotext -f 1 -l 1 refs/paper.pdf paper_abstract.txt
```

### Using grep to find abstracts

```bash
# Quick abstract extraction
pdftotext -l 2 refs/paper.pdf - | \
  grep -A 20 -i "abstract" | \
  head -30
```

## Batch Processing Script Example

Create `scripts/extract_all_abstracts.sh`:

```bash
#!/bin/bash
# Extract abstracts from all PDFs in refs/

REFS_DIR="refs"
OUTPUT="data/abstracts_$(date +%Y%m%d).json"

find "$REFS_DIR" -name "*.pdf" -type f | \
while read pdf; do
  echo "Processing: $pdf"
  python code/pdf_extract.py "$pdf" --output temp_abs.json
done

echo "Done! Check $OUTPUT"
```

## Python Integration

Use in your own scripts:

```python
from code.pdf_extract import extract_text, find_abstract, extract_metadata
from pathlib import Path

pdf_path = Path('refs/paper.pdf')
text = extract_text(pdf_path, max_pages=3)
abstract = find_abstract(text)
meta = extract_metadata(text, pdf_path)

print(f"Title: {meta.get('title')}")
print(f"Abstract: {abstract}")
```
