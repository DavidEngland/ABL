#!/usr/bin/env python3
import re
import argparse
import hashlib
import os
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass, field

# PDF Library handling
try:
    import fitz  # PyMuPDF (Highly recommended for better results)
    HAS_FITZ = True
except ImportError:
    HAS_FITZ = False

@dataclass
class Reference:
    authors: List[str] = field(default_factory=list)
    title: str = ""
    year: Optional[str] = None
    journal: Optional[str] = None
    volume: Optional[str] = None
    pages: Optional[str] = None
    doi: Optional[str] = None
    url: Optional[str] = None
    abstract: Optional[str] = None
    source_file: Optional[str] = None  # New field for relative path
    entry_type: str = "article"
    
    def to_bibtex(self, cite_key: str = None) -> str:
        if not cite_key:
            cite_key = self.generate_cite_key()
        
        lines = [f"@{self.entry_type}{{{cite_key},"]
        
        # BibTeX standard: ALL authors separated by " and "
        if self.authors:
            # Clean authors of extra whitespace
            clean_authors = [a.strip() for a in self.authors if len(a.strip()) > 1]
            lines.append(f"  author = {{{' and '.join(clean_authors)}}},")
        
        if self.title: lines.append(f"  title = {{{self.title}}},")
        if self.year: lines.append(f"  year = {{{self.year}}},")
        if self.journal: lines.append(f"  journal = {{{self.journal}}},")
        if self.abstract: lines.append(f"  abstract = {{{self.abstract}}},")
        if self.source_file: lines.append(f"  file = {{{self.source_file}}},")
        if self.doi: lines.append(f"  doi = {{{self.doi}}},")
        if self.url: lines.append(f"  url = {{{self.url}}},")
        
        lines.append("}")
        return "\n".join(lines)

    def generate_cite_key(self) -> str:
        if not self.authors or not self.year:
            return "ref_" + hashlib.md5(self.title.encode()).hexdigest()[:8]
        # Get last name of first author
        last_name = self.authors[0].split()[-1].lower()
        last_name = re.sub(r'[^a-z]', '', last_name)
        return f"{last_name}{self.year}"

class ReferenceExtractor:
    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def _get_relative_path(self, absolute_path: Path) -> str:
        """Calculates path relative to current working directory."""
        try:
            return str(absolute_path.relative_to(Path.cwd()))
        except ValueError:
            return str(absolute_path)

    def extract_from_pdf(self, pdf_path: Path) -> List[Reference]:
        if not HAS_FITZ:
            if self.verbose: print(f"Skipping {pdf_path}: PyMuPDF not installed.")
            return []

        try:
            doc = fitz.open(pdf_path)
            # We usually find title/abstract on the first two pages
            header_text = "\n".join(page.get_text() for page in doc[:2])
            full_text = "\n".join(page.get_text() for page in doc)
            
            # 1. Extract Abstract
            abstract = self._extract_abstract(header_text)
            
            # 2. Extract References Section
            refs_text = self._extract_references_section(full_text)
            parsed_refs = self._parse_reference_list(refs_text)
            
            # Enrich with source metadata
            for ref in parsed_refs:
                ref.source_file = self._get_relative_path(pdf_path)
                if not ref.abstract: # If the reference is the paper itself
                    ref.abstract = abstract

            doc.close()
            return parsed_refs
        except Exception as e:
            if self.verbose: print(f"Error processing {pdf_path}: {e}")
            return []

    def _extract_abstract(self, text: str) -> Optional[str]:
        """Looks for text between 'Abstract' and the next major section."""
        pattern = r'(?i)abstract[:\s\n]+(.*?)(?=\n(?:1\.?\s+|Introduction|Keywords|I\.\s+))'
        match = re.search(pattern, text, re.DOTALL)
        if match:
            # Clean up newlines and extra spaces
            return re.sub(r'\s+', ' ', match.group(1)).strip()
        return None

    def _extract_references_section(self, text: str) -> str:
        patterns = [
            r'(?i)\nReferences\s*\n(.*)',
            r'(?i)\nBibliography\s*\n(.*)',
            r'(?i)\nLITERATURE CITED\s*\n(.*)'
        ]
        for p in patterns:
            match = re.search(p, text, re.DOTALL)
            if match: return match.group(1)
        return text[int(len(text)*0.8):] # Fallback

    def _parse_reference_list(self, text: str) -> List[Reference]:
        # Improved splitting for [1], 1., or Author (Year)
        entries = re.split(r'\n(?:\[\d+\]|\d+\.)\s+', text)
        refs = []
        for entry in entries:
            if len(entry.strip()) < 15: continue
            refs.append(self._parse_single_reference(entry.strip()))
        return refs

    def _parse_single_reference(self, ref_text: str) -> Reference:
        ref = Reference()
        
        # Extract Year
        year_match = re.search(r'\b(19|20)\d{2}\b', ref_text)
        if year_match:
            ref.year = year_match.group(0)

        # Better Author Extraction: looks for "Name, F." or "First Last" before the year/title
        author_part = ref_text.split(ref.year)[0] if ref.year else ref_text[:100]
        # Regex for common academic author formats
        authors = re.findall(r'([A-Z][a-z]+,?\s+[A-Z]\.(?:\s+[A-Z]\.)?|[A-Z][a-z]+\s+[A-Z][a-z]+)', author_part)
        ref.authors = authors if authors else ["Unknown Author"]

        # Title heuristic: Text between quotes or the first long string after authors
        title_match = re.search(r'["“](.*?)["”]', ref_text)
        if title_match:
            ref.title = title_match.group(1)
        else:
            # Take the longest segment that doesn't look like authors or journal
            segments = re.split(r'\.\s+', ref_text)
            for seg in segments:
                if len(seg) > 20 and not any(a in seg for a in ref.authors):
                    ref.title = seg.strip()
                    break
        
        return ref

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', nargs='+')
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    extractor = ReferenceExtractor(verbose=args.verbose)
    all_refs = []

    for input_str in args.inputs:
        path = Path(input_str)
        files = path.glob('**/*.pdf') if path.is_dir() else [path]

        for f in files:
            if f.suffix.lower() == '.pdf':
                all_refs.extend(extractor.extract_from_pdf(f))

    with open(args.output, 'w', encoding='utf-8') as f:
        output = "\n\n".join(r.to_bibtex() for r in all_refs)
        f.write(output)
    
    print(f"Done! Saved to {args.output}")

if __name__ == "__main__":
    main()