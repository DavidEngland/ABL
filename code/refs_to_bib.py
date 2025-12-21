#!/usr/bin/env python3
"""Extract references from multiple file formats and convert to BibTeX.

Usage:
    python refs_to_bib.py refs/*.pdf --output references.bib
    python refs_to_bib.py paper.tex --output extracted.bib --merge existing.bib
    python refs_to_bib.py refs/ --recursive --formats pdf,tex,md
"""

import re
import json
import argparse
import hashlib
from pathlib import Path
from typing import List, Dict, Optional, Set
from dataclasses import dataclass, field
from collections import defaultdict

# Try to import PDF libraries
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


@dataclass
class Reference:
    """Structured reference entry."""
    authors: List[str] = field(default_factory=list)
    title: str = ""
    year: Optional[str] = None
    journal: Optional[str] = None
    volume: Optional[str] = None
    pages: Optional[str] = None
    doi: Optional[str] = None
    url: Optional[str] = None
    booktitle: Optional[str] = None
    publisher: Optional[str] = None
    raw_text: str = ""
    entry_type: str = "article"  # article, book, inproceedings, etc.
    
    def to_bibtex(self, cite_key: str = None) -> str:
        """Convert to BibTeX format."""
        if not cite_key:
            cite_key = self.generate_cite_key()
        
        lines = [f"@{self.entry_type}{{{cite_key},"]
        
        # Author formatting
        if self.authors:
            author_str = " and ".join(self.authors)
            lines.append(f"  author = {{{author_str}}},")
        
        # Required fields
        if self.title:
            lines.append(f"  title = {{{self.title}}},")
        if self.year:
            lines.append(f"  year = {{{self.year}}},")
        
        # Optional fields
        if self.journal:
            lines.append(f"  journal = {{{self.journal}}},")
        if self.booktitle:
            lines.append(f"  booktitle = {{{self.booktitle}}},")
        if self.volume:
            lines.append(f"  volume = {{{self.volume}}},")
        if self.pages:
            lines.append(f"  pages = {{{self.pages}}},")
        if self.publisher:
            lines.append(f"  publisher = {{{self.publisher}}},")
        if self.doi:
            lines.append(f"  doi = {{{self.doi}}},")
        if self.url:
            lines.append(f"  url = {{{self.url}}},")
        
        lines.append("}")
        return "\n".join(lines)
    
    def generate_cite_key(self) -> str:
        """Generate citation key (author_year format)."""
        if not self.authors or not self.year:
            # Fallback to hash of title
            return "ref_" + hashlib.md5(self.title.encode()).hexdigest()[:8]
        
        # Extract last name of first author
        first_author = self.authors[0].split()[-1].lower()
        first_author = re.sub(r'[^a-z]', '', first_author)
        
        return f"{first_author}{self.year}"


class ReferenceExtractor:
    """Extract references from various file formats."""
    
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.references: List[Reference] = []
    
    def extract_from_pdf(self, pdf_path: Path) -> List[Reference]:
        """Extract references from PDF."""
        if self.verbose:
            print(f"Extracting from PDF: {pdf_path.name}...")
        
        try:
            # Extract full text
            if HAS_FITZ:
                doc = fitz.open(pdf_path)
                text = "\n".join(page.get_text() for page in doc)
                doc.close()
            elif HAS_PYPDF2:
                text = ""
                with open(pdf_path, 'rb') as f:
                    reader = PyPDF2.PdfReader(f)
                    text = "\n".join(page.extract_text() for page in reader.pages)
            else:
                raise ImportError("No PDF library available")
            
            # Find references section
            refs_text = self._extract_references_section(text)
            
            # Parse individual references
            return self._parse_reference_list(refs_text)
            
        except Exception as e:
            if self.verbose:
                print(f"  Error: {e}")
            return []
    
    def extract_from_tex(self, tex_path: Path) -> List[Reference]:
        """Extract references from LaTeX file."""
        if self.verbose:
            print(f"Extracting from TeX: {tex_path.name}...")
        
        text = tex_path.read_text(encoding='utf-8', errors='ignore')
        refs = []
        
        # Extract \bibitem entries
        bibitem_pattern = r'\\bibitem(?:\[.*?\])?\{([^}]+)\}\s*(.+?)(?=\\bibitem|\\end\{thebibliography\}|$)'
        for match in re.finditer(bibitem_pattern, text, re.DOTALL):
            cite_key, ref_text = match.groups()
            ref = self._parse_single_reference(ref_text.strip())
            ref.raw_text = ref_text.strip()
            refs.append(ref)
        
        # Extract existing \bibliography or .bib file references
        bib_pattern = r'\\bibliography\{([^}]+)\}'
        for match in re.finditer(bib_pattern, text):
            bib_file = match.group(1)
            if not bib_file.endswith('.bib'):
                bib_file += '.bib'
            bib_path = tex_path.parent / bib_file
            if bib_path.exists():
                refs.extend(self.extract_from_bib(bib_path))
        
        return refs
    
    def extract_from_bib(self, bib_path: Path) -> List[Reference]:
        """Extract references from existing BibTeX file."""
        if self.verbose:
            print(f"Reading BibTeX: {bib_path.name}...")
        
        text = bib_path.read_text(encoding='utf-8', errors='ignore')
        refs = []
        
        # Parse BibTeX entries
        entry_pattern = r'@(\w+)\{([^,]+),\s*(.+?)\n\}'
        for match in re.finditer(entry_pattern, text, re.DOTALL):
            entry_type, cite_key, fields = match.groups()
            ref = Reference(entry_type=entry_type.lower())
            
            # Parse fields
            field_pattern = r'(\w+)\s*=\s*\{([^}]+)\}|(\w+)\s*=\s*"([^"]+)"'
            for field_match in re.finditer(field_pattern, fields):
                if field_match.group(1):
                    field_name, field_value = field_match.group(1), field_match.group(2)
                else:
                    field_name, field_value = field_match.group(3), field_match.group(4)
                
                field_name = field_name.lower()
                
                if field_name == 'author':
                    ref.authors = [a.strip() for a in field_value.split(' and ')]
                elif field_name == 'title':
                    ref.title = field_value
                elif field_name == 'year':
                    ref.year = field_value
                elif field_name == 'journal':
                    ref.journal = field_value
                elif field_name == 'volume':
                    ref.volume = field_value
                elif field_name == 'pages':
                    ref.pages = field_value
                elif field_name == 'doi':
                    ref.doi = field_value
                elif field_name == 'url':
                    ref.url = field_value
                elif field_name == 'booktitle':
                    ref.booktitle = field_value
                elif field_name == 'publisher':
                    ref.publisher = field_value
            
            refs.append(ref)
        
        return refs
    
    def extract_from_md(self, md_path: Path) -> List[Reference]:
        """Extract references from Markdown file."""
        if self.verbose:
            print(f"Extracting from Markdown: {md_path.name}...")
        
        text = md_path.read_text(encoding='utf-8', errors='ignore')
        refs = []
        
        # Find reference sections (common patterns)
        ref_section_patterns = [
            r'##?\s*References?\s*\n(.+?)(?=##|\Z)',
            r'##?\s*Bibliography\s*\n(.+?)(?=##|\Z)',
            r'##?\s*Citations?\s*\n(.+?)(?=##|\Z)',
        ]
        
        refs_text = ""
        for pattern in ref_section_patterns:
            match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
            if match:
                refs_text = match.group(1)
                break
        
        if refs_text:
            # Parse numbered or bulleted list of references
            ref_items = re.split(r'\n(?:\d+\.|[\-\*])\s+', refs_text)
            for item in ref_items:
                if len(item.strip()) > 20:
                    ref = self._parse_single_reference(item.strip())
                    refs.append(ref)
        
        # Also extract DOIs from anywhere in document
        doi_pattern = r'(?:https?://)?(?:doi\.org/|DOI:\s*)(10\.\d{4,}/[^\s\)]+)'
        for match in re.finditer(doi_pattern, text, re.IGNORECASE):
            doi = match.group(1)
            # Create minimal reference entry
            ref = Reference(doi=doi, raw_text=f"DOI: {doi}")
            refs.append(ref)
        
        return refs
    
    def extract_from_html(self, html_path: Path) -> List[Reference]:
        """Extract references from HTML file."""
        if self.verbose:
            print(f"Extracting from HTML: {html_path.name}...")
        
        text = html_path.read_text(encoding='utf-8', errors='ignore')
        
        # Remove HTML tags
        text = re.sub(r'<script[^>]*>.*?</script>', '', text, flags=re.DOTALL)
        text = re.sub(r'<style[^>]*>.*?</style>', '', text, flags=re.DOTALL)
        text = re.sub(r'<[^>]+>', ' ', text)
        
        # Find references section
        refs_text = self._extract_references_section(text)
        
        return self._parse_reference_list(refs_text)
    
    def _extract_references_section(self, text: str) -> str:
        """Extract the References/Bibliography section from text."""
        patterns = [
            r'References?\s*\n(.+?)(?=\n\s*(?:Appendix|Acknowledgments?|Figure|Table|\Z))',
            r'Bibliography\s*\n(.+?)(?=\n\s*(?:Appendix|Acknowledgments?|Figure|Table|\Z))',
            r'REFERENCES?\s*\n(.+?)(?=\n\s*(?:APPENDIX|ACKNOWLEDGMENTS?|FIGURE|TABLE|\Z))',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
            if match:
                return match.group(1)
        
        # Fallback: last 20% of document
        return text[int(len(text)*0.8):]
    
    def _parse_reference_list(self, refs_text: str) -> List[Reference]:
        """Parse a block of reference text into individual references."""
        refs = []
        
        # Split by numbered references (common pattern: [1], 1., (1), etc.)
        # Also handle cases with author names at start
        ref_patterns = [
            r'\[\d+\]\s*(.+?)(?=\[\d+\]|\Z)',  # [1] format
            r'^\d+\.\s*(.+?)(?=^\d+\.|\Z)',     # 1. format
            r'^\(\d+\)\s*(.+?)(?=^\(\d+\)|\Z)', # (1) format
            r'^[A-Z][^.]+\.\s*\(\d{4}\)\.(.+?)(?=^[A-Z][^.]+\.\s*\(\d{4}\)\.|\Z)',  # Author (year) format
        ]
        
        for pattern in ref_patterns:
            matches = list(re.finditer(pattern, refs_text, re.MULTILINE | re.DOTALL))
            if matches:
                for match in matches:
                    ref_text = match.group(1).strip()
                    if len(ref_text) > 20:  # Filter out spurious matches
                        ref = self._parse_single_reference(ref_text)
                        ref.raw_text = ref_text
                        refs.append(ref)
                break
        
        return refs
    
    def _parse_single_reference(self, ref_text: str) -> Reference:
        """Parse a single reference string into structured fields."""
        ref = Reference(raw_text=ref_text)
        
        # Extract DOI
        doi_match = re.search(r'(?:https?://)?(?:doi\.org/|DOI:\s*)(10\.\d{4,}/[^\s,\.]+)', ref_text, re.IGNORECASE)
        if doi_match:
            ref.doi = doi_match.group(1)
        
        # Extract year (4-digit number in parens or standalone)
        year_match = re.search(r'\(?(19|20)\d{2}\)?', ref_text)
        if year_match:
            ref.year = re.sub(r'[()]', '', year_match.group(0))
        
        # Extract title (text in quotes or italics)
        title_patterns = [
            r'"([^"]{20,})"',
            r'«([^»]{20,})»',
            r'\*([^\*]{20,})\*',
        ]
        for pattern in title_patterns:
            title_match = re.search(pattern, ref_text)
            if title_match:
                ref.title = title_match.group(1).strip()
                break
        
        # If no quoted title, try to extract longest sentence-like segment
        if not ref.title:
            # Remove year and author clues first
            clean_text = re.sub(r'\(?\d{4}\)?', '', ref_text)
            sentences = re.split(r'[\.]\s+', clean_text)
            for sent in sentences:
                if 30 < len(sent) < 300:
                    ref.title = sent.strip()
                    break
        
        # Extract journal name (common patterns)
        journal_patterns = [
            r'\b([A-Z][a-z]+ (?:of|and) [A-Z][a-z]+(?:\s+[A-Z][a-z]+)?)\b',  # "Journal of X"
            r'\b([A-Z][a-z]+\s+[A-Z][a-z]+(?:\s+[A-Z][a-z]+)?)\s+\d+',  # "Name Volume"
        ]
        for pattern in journal_patterns:
            journal_match = re.search(pattern, ref_text)
            if journal_match:
                candidate = journal_match.group(1)
                if 'Journal' in candidate or 'Review' in candidate or 'Letters' in candidate:
                    ref.journal = candidate
                    break
        
        # Extract volume and pages
        vol_pages_match = re.search(r'(\d+)\s*[\(:]\s*(\d+(?:[-–]\d+)?)', ref_text)
        if vol_pages_match:
            ref.volume = vol_pages_match.group(1)
            ref.pages = vol_pages_match.group(2).replace('–', '--')
        
        # Extract authors (heuristic: Names before year or title)
        if ref.year:
            author_text = ref_text.split(ref.year)[0]
        elif ref.title:
            author_text = ref_text.split(ref.title)[0]
        else:
            author_text = ref_text[:200]
        
        # Parse author names (simplified)
        author_candidates = re.findall(r'\b([A-Z][a-z]+(?:\s+[A-Z]\.)*)\b', author_text)
        if author_candidates:
            ref.authors = author_candidates[:5]  # Limit to first 5 authors
        
        # Determine entry type
        if 'Proceedings' in ref_text or 'Conference' in ref_text:
            ref.entry_type = 'inproceedings'
        elif any(word in ref_text for word in ['Book', 'Press', 'Publishers']):
            ref.entry_type = 'book'
        
        return ref


def merge_bibtex_files(files: List[Path]) -> str:
    """Merge multiple BibTeX files, removing duplicates."""
    extractor = ReferenceExtractor()
    all_refs = []
    seen_keys = set()
    
    for bib_file in files:
        refs = extractor.extract_from_bib(bib_file)
        for ref in refs:
            key = ref.generate_cite_key()
            if key not in seen_keys:
                all_refs.append(ref)
                seen_keys.add(key)
    
    return "\n\n".join(ref.to_bibtex() for ref in all_refs)


def main():
    parser = argparse.ArgumentParser(description='Extract references and convert to BibTeX')
    parser.add_argument('inputs', nargs='+', help='Input files or directories')
    parser.add_argument('-o', '--output', required=True, help='Output .bib file')
    parser.add_argument('-f', '--formats', default='pdf,tex,md,html,bib', 
                        help='File formats to process (comma-separated)')
    parser.add_argument('-r', '--recursive', action='store_true', 
                        help='Recursively search directories')
    parser.add_argument('--merge', help='Merge with existing .bib file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--deduplicate', action='store_true', 
                        help='Remove duplicate entries (by DOI/title)')
    
    args = parser.parse_args()
    
    # Check PDF dependencies
    if 'pdf' in args.formats.lower() and not (HAS_FITZ or HAS_PYPDF2):
        print("Warning: No PDF library found. Install pymupdf or PyPDF2 to process PDFs")
        args.formats = args.formats.replace('pdf', '')
    
    # Collect input files
    formats = [f.strip() for f in args.formats.split(',')]
    input_files = []
    
    for input_path in args.inputs:
        path = Path(input_path)
        if path.is_dir():
            for fmt in formats:
                pattern = f"**/*.{fmt}" if args.recursive else f"*.{fmt}"
                input_files.extend(path.glob(pattern))
        elif path.is_file():
            input_files.append(path)
    
    if not input_files:
        print("No input files found!")
        return 1
    
    if args.verbose:
        print(f"Found {len(input_files)} files to process")
    
    # Extract references
    extractor = ReferenceExtractor(verbose=args.verbose)
    all_refs: List[Reference] = []
    
    for file_path in input_files:
        suffix = file_path.suffix.lower()
        
        if suffix == '.pdf':
            refs = extractor.extract_from_pdf(file_path)
        elif suffix == '.tex':
            refs = extractor.extract_from_tex(file_path)
        elif suffix == '.bib':
            refs = extractor.extract_from_bib(file_path)
        elif suffix == '.md':
            refs = extractor.extract_from_md(file_path)
        elif suffix in ['.html', '.htm']:
            refs = extractor.extract_from_html(file_path)
        else:
            continue
        
        all_refs.extend(refs)
    
    # Deduplicate if requested
    if args.deduplicate:
        unique_refs = []
        seen = set()
        
        for ref in all_refs:
            # Use DOI as primary key, fallback to normalized title
            if ref.doi:
                key = ref.doi.lower()
            else:
                key = re.sub(r'\s+', '', ref.title.lower())[:100]
            
            if key and key not in seen:
                unique_refs.append(ref)
                seen.add(key)
        
        if args.verbose:
            print(f"Deduplicated: {len(all_refs)} → {len(unique_refs)} references")
        
        all_refs = unique_refs
    
    # Merge with existing file if requested
    if args.merge and Path(args.merge).exists():
        existing_refs = extractor.extract_from_bib(Path(args.merge))
        # Add existing refs that aren't duplicates
        existing_keys = {ref.generate_cite_key() for ref in all_refs}
        for ref in existing_refs:
            if ref.generate_cite_key() not in existing_keys:
                all_refs.append(ref)
        
        if args.verbose:
            print(f"Merged with {len(existing_refs)} existing references")
    
    # Generate BibTeX output
    bibtex_entries = []
    seen_cite_keys = set()
    
    for ref in all_refs:
        cite_key = ref.generate_cite_key()
        
        # Handle duplicate cite keys
        if cite_key in seen_cite_keys:
            counter = 2
            while f"{cite_key}_{counter}" in seen_cite_keys:
                counter += 1
            cite_key = f"{cite_key}_{counter}"
        
        seen_cite_keys.add(cite_key)
        bibtex_entries.append(ref.to_bibtex(cite_key))
    
    # Write output
    output_path = Path(args.output)
    output_path.write_text("\n\n".join(bibtex_entries), encoding='utf-8')
    
    print(f"✓ Wrote {len(bibtex_entries)} references to {output_path}")
    
    # Statistics
    if args.verbose:
        stats = defaultdict(int)
        for ref in all_refs:
            stats[ref.entry_type] += 1
        
        print("\nEntry type breakdown:")
        for entry_type, count in sorted(stats.items()):
            print(f"  {entry_type}: {count}")
    
    return 0


if __name__ == '__main__':
    exit(main())
