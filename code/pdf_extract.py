#!/usr/bin/env python3
"""Extract abstracts and metadata from PDF papers.

Usage:
    python pdf_extract.py paper.pdf
    python pdf_extract.py refs/*.pdf --output summaries.json
"""

import re
import json
import argparse
from pathlib import Path
from typing import Optional, Dict, List

try:
    import fitz  # PyMuPDF
    HAS_FITZ = True
except ImportError:
    HAS_FITZ = False

try:
    import PyPDF2
    HAS_PYPDF2 = True
except ImportError:
    HAS_PYPDF2 = False


def extract_text_fitz(pdf_path: Path, max_pages: int = 3) -> str:
    """Extract text using PyMuPDF (faster, better formatting)."""
    doc = fitz.open(pdf_path)
    text = ""
    for page_num in range(min(max_pages, len(doc))):
        text += doc[page_num].get_text()
    doc.close()
    return text


def extract_text_pypdf2(pdf_path: Path, max_pages: int = 3) -> str:
    """Extract text using PyPDF2 (fallback)."""
    text = ""
    with open(pdf_path, 'rb') as f:
        reader = PyPDF2.PdfReader(f)
        for page_num in range(min(max_pages, len(reader.pages))):
            text += reader.pages[page_num].extract_text()
    return text


def extract_text(pdf_path: Path, max_pages: int = 3) -> str:
    """Extract text from PDF using available library."""
    if HAS_FITZ:
        return extract_text_fitz(pdf_path, max_pages)
    elif HAS_PYPDF2:
        return extract_text_pypdf2(pdf_path, max_pages)
    else:
        raise ImportError("Install PyMuPDF (pip install pymupdf) or PyPDF2 (pip install PyPDF2)")


def find_abstract(text: str) -> Optional[str]:
    """Extract abstract using common patterns."""
    # Normalize whitespace
    text = re.sub(r'\s+', ' ', text)
    
    # Pattern 1: "Abstract" followed by text until "Introduction" or next section
    patterns = [
        r'Abstract\s*[:\-]?\s*(.+?)(?:Introduction|Keywords|1\.|§\s*1|INTRODUCTION)',
        r'ABSTRACT\s*[:\-]?\s*(.+?)(?:INTRODUCTION|Keywords|1\.|§\s*1)',
        r'Summary\s*[:\-]?\s*(.+?)(?:Introduction|Keywords|1\.|§\s*1)',
    ]
    
    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE | re.DOTALL)
        if match:
            abstract = match.group(1).strip()
            # Clean up common artifacts
            abstract = re.sub(r'\s+', ' ', abstract)
            abstract = re.sub(r'^\W+', '', abstract)
            # Limit length (abstracts typically < 2000 chars)
            if 50 < len(abstract) < 2000:
                return abstract
    
    # Fallback: first substantial paragraph after removing header junk
    lines = text.split('\n')
    substantial = []
    for line in lines[5:50]:  # skip first few lines (headers)
        line = line.strip()
        if len(line) > 100 and not line.isupper():
            substantial.append(line)
            if len(' '.join(substantial)) > 500:
                break
    
    if substantial:
        candidate = ' '.join(substantial)
        if 100 < len(candidate) < 2000:
            return candidate
    
    return None


def extract_metadata(text: str, pdf_path: Path) -> Dict:
    """Extract paper metadata (title, authors, year, DOI)."""
    metadata = {'file': pdf_path.name}
    
    # Extract title (usually largest text on first page, or first line)
    lines = [l.strip() for l in text.split('\n') if l.strip()]
    if lines:
        # Heuristic: title is often first substantial line
        for line in lines[:20]:
            if 20 < len(line) < 200 and not line.isupper():
                metadata['title'] = line
                break
    
    # Extract year (4-digit number, typically 19xx or 20xx)
    year_match = re.search(r'\b(19|20)\d{2}\b', text[:2000])
    if year_match:
        metadata['year'] = year_match.group(0)
    
    # Extract DOI
    doi_match = re.search(r'10\.\d{4,}/[^\s]+', text[:3000])
    if doi_match:
        metadata['doi'] = doi_match.group(0)
    
    # Extract journal/venue (common patterns)
    journal_patterns = [
        r'(?:published in|appeared in)\s+([A-Z][^.]{10,80})',
        r'Journal of\s+([^.]{5,50})',
        r'([A-Z][a-z]+\s+[A-Z][a-z]+)\s+(?:Journal|Review|Letters)',
    ]
    for pattern in journal_patterns:
        match = re.search(pattern, text[:2000], re.IGNORECASE)
        if match:
            metadata['venue'] = match.group(1).strip()
            break
    
    return metadata


def process_pdf(pdf_path: Path, verbose: bool = False) -> Dict:
    """Extract abstract and metadata from a single PDF."""
    if verbose:
        print(f"Processing: {pdf_path.name}...")
    
    try:
        # Extract text from first few pages
        text = extract_text(pdf_path, max_pages=3)
        
        # Find abstract
        abstract = find_abstract(text)
        
        # Extract metadata
        metadata = extract_metadata(text, pdf_path)
        
        result = {
            **metadata,
            'abstract': abstract if abstract else '[Abstract not found]',
            'status': 'success' if abstract else 'partial'
        }
        
        if verbose:
            print(f"  ✓ {result['status']}: {len(abstract) if abstract else 0} chars")
        
        return result
        
    except Exception as e:
        if verbose:
            print(f"  ✗ Error: {e}")
        return {
            'file': pdf_path.name,
            'status': 'error',
            'error': str(e)
        }


def main():
    parser = argparse.ArgumentParser(description='Extract abstracts from PDF papers')
    parser.add_argument('pdfs', nargs='+', help='PDF file(s) to process')
    parser.add_argument('-o', '--output', help='Output JSON file (default: print to stdout)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--max-pages', type=int, default=3, help='Max pages to scan (default: 3)')
    
    args = parser.parse_args()
    
    # Check dependencies
    if not (HAS_FITZ or HAS_PYPDF2):
        print("Error: No PDF library found. Install one of:")
        print("  pip install pymupdf  (recommended, faster)")
        print("  pip install PyPDF2   (fallback)")
        return 1
    
    # Process PDFs
    results = []
    for pdf_pattern in args.pdfs:
        pdf_paths = list(Path('.').glob(pdf_pattern)) if '*' in pdf_pattern else [Path(pdf_pattern)]
        
        for pdf_path in pdf_paths:
            if not pdf_path.exists():
                print(f"Warning: {pdf_path} not found, skipping...")
                continue
            
            result = process_pdf(pdf_path, verbose=args.verbose)
            results.append(result)
    
    # Output
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Wrote {len(results)} summaries to {args.output}")
    else:
        # Pretty print to stdout
        for r in results:
            print(f"\n{'='*70}")
            print(f"File: {r['file']}")
            if 'title' in r:
                print(f"Title: {r.get('title', 'N/A')}")
            if 'year' in r:
                print(f"Year: {r.get('year', 'N/A')}")
            if 'doi' in r:
                print(f"DOI: {r.get('doi', 'N/A')}")
            print(f"\nAbstract:\n{r.get('abstract', '[Not found]')}")
    
    return 0


if __name__ == '__main__':
    exit(main())
