#!/usr/bin/env python3
import re
import argparse
import hashlib
import time
from pathlib import Path
from typing import List, Optional
from dataclasses import dataclass, field

# New requirement: pip install requests
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    import fitz  # PyMuPDF
    HAS_FITZ = True
except ImportError:
    HAS_FITZ = False

@dataclass
class Reference:
    authors: List[str] = field(default_factory=list)
    title: str = ""
    year: Optional[str] = None
    journal: Optional[str] = None
    doi: Optional[str] = None
    abstract: Optional[str] = None
    source_file: Optional[str] = None
    entry_type: str = "article"
    raw_bibtex: Optional[str] = None  # To store official BibTeX from API

    def to_bibtex(self, cite_key: str = None) -> str:
        # If we successfully fetched the official BibTeX from DOI, return that
        if self.raw_bibtex:
            # Preserve original entry type and only replace the cite key.
            patched = re.sub(
                r'^(@\w+\{)[^,]+,',
                r'\1' + f"{cite_key or self.generate_cite_key()}" + ',',
                self.raw_bibtex
            )
            # Inject extra fields (file/abstract) if available.
            extra = ""
            if self.source_file:
                extra += f"  file = {{{self.source_file}}},\n"
            if self.abstract:
                extra += f"  abstract = {{{self.abstract}}},\n"
            if extra:
                patched = re.sub(r'\n\}$', f"\n{extra}" + "}", patched)
            return patched

        # Fallback to manual construction
        if not cite_key:
            cite_key = self.generate_cite_key()

        lines = [f"@{self.entry_type}{{{cite_key},"]
        if self.authors:
            lines.append(f"  author = {{{' and '.join(self.authors)}}},")
        if self.title: lines.append(f"  title = {{{self.title}}},")
        if self.year: lines.append(f"  year = {{{self.year}}},")
        if self.journal: lines.append(f"  journal = {{{self.journal}}},")
        if self.abstract: lines.append(f"  abstract = {{{self.abstract}}},")
        if self.source_file: lines.append(f"  file = {{{self.source_file}}},")
        if self.doi: lines.append(f"  doi = {{{self.doi}}},")
        lines.append("}")
        return "\n".join(lines)

    def generate_cite_key(self) -> str:
        if not self.authors or not self.year:
            return "ref_" + hashlib.md5(self.title.encode()).hexdigest()[:8]
        last_name = self.authors[0].split()[-1].lower()
        last_name = re.sub(r'[^a-z]', '', last_name)
        return f"{last_name}{self.year}"

class ReferenceExtractor:
    def __init__(self, verbose: bool = False, email: Optional[str] = None):
        self.verbose = verbose
        ua = "ABL-DOI/1.0"
        if email:
            ua += f" ({email})"
        self.headers = {'Accept': 'application/x-bibtex', 'User-Agent': ua}
        if self.verbose and not HAS_REQUESTS:
            print("Note: requests not installed; DOI verification disabled. Install with: pip install requests")

    def find_doi(self, text: str) -> Optional[str]:
        """Find and clean a DOI from arbitrary text."""
        m = re.search(r'(10\.\d{4,9}/[-._;()/:A-Z0-9]+)', text, re.I)
        if not m:
            return None
        doi = m.group(1)
        # Trim common trailing punctuation/brackets
        doi = doi.rstrip(').,;] ').lstrip('(')
        return doi

    def fetch_metadata_by_doi(self, doi: str) -> Optional[str]:
        """Queries the DOI resolver for an official BibTeX entry."""
        if not doi:
            return None
        if not HAS_REQUESTS:
            if self.verbose:
                print("    Skipping DOI lookup (requests not installed). Run: pip install requests")
            return None
        clean_doi = self.find_doi(doi)
        if not clean_doi:
            return None
        url = f"https://doi.org/{clean_doi}"
        try:
            if self.verbose: print(f"    Searching internet for DOI: {clean_doi}...")
            response = requests.get(url, headers=self.headers, timeout=10)
            if response.status_code == 200:
                return response.text
            if self.verbose: print(f"    DOI lookup HTTP {response.status_code}")
        except Exception as e:
            if self.verbose: print(f"    Web search failed: {e}")
        return None

    def extract_from_pdf(self, pdf_path: Path) -> List[Reference]:
        if not HAS_FITZ: return []
        try:
            pdf_path = pdf_path.resolve()
            doc = fitz.open(pdf_path)
            # Prefer document metadata title when available
            meta = doc.metadata or {}
            header_text = "\n".join(page.get_text() for page in doc[:2])
            full_text = "\n".join(page.get_text() for page in doc)

            abstract = self._extract_abstract(header_text)

            # Look for DOI in the first few pages (often in footer/header)
            found_doi = self.find_doi(header_text)

            try:
                rel = pdf_path.relative_to(Path.cwd())
            except ValueError:
                rel = pdf_path

            paper = Reference(
                title=(meta.get('title') or pdf_path.stem),
                doi=found_doi,
                abstract=abstract,
                source_file=str(rel)
            )

            if paper.doi:
                official_bib = self.fetch_metadata_by_doi(paper.doi)
                if official_bib:
                    paper.raw_bibtex = official_bib
                    if self.verbose: print(f"    ✓ Verified via DOI")

            # 2. Extract bibliography (References at end)
            refs_text = self._extract_references_section(full_text)
            parsed_refs = self._parse_reference_list(refs_text)

            all_found = [paper] + parsed_refs
            doc.close()
            return all_found
        except Exception as e:
            if self.verbose: print(f"Error: {e}")
            return []

    def _extract_abstract(self, text: str) -> Optional[str]:
        match = re.search(r'(?i)abstract[:\s\n]+(.*?)(?=\n(?:1\.?\s+|Introduction|Keywords))', text, re.DOTALL)
        return re.sub(r'\s+', ' ', match.group(1)).strip() if match else None

    def _extract_references_section(self, text: str) -> str:
        match = re.search(r'(?i)\nReferences\s*\n(.*)', text, re.DOTALL)
        return match.group(1) if match else text[int(len(text)*0.8):]

    def _parse_reference_list(self, text: str) -> List[Reference]:
        entries = re.split(r'\n(?:\[\d+\]|\d+\.)\s+', text)
        return [self._parse_single_reference(e.strip()) for e in entries if len(e.strip()) > 15]

    def _parse_single_reference(self, ref_text: str) -> Reference:
        ref = Reference()
        year_match = re.search(r'\b(19|20)\d{2}\b', ref_text)
        if year_match: ref.year = year_match.group(0)

        # Look for DOI within the reference text itself
        doi_match = re.search(r'10\.\d{4,}/[^\s]+', ref_text)
        if doi_match:
            ref.doi = doi_match.group(0).rstrip('.')
            # Attempt to verify this reference online too
            official = self.fetch_metadata_by_doi(ref.doi)
            if official:
                ref.raw_bibtex = official
                return ref

        # Fallback to basic extraction if no DOI or lookup failed
        author_part = ref_text.split(ref.year)[0] if ref.year else ref_text[:50]
        ref.authors = re.findall(r'([A-Z][a-z]+,?\s+[A-Z]\.)', author_part) or ["Unknown"]
        return ref

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', nargs='*', default=['refs'], help="PDF files or folders (default: ./refs)")
    parser.add_argument('-o', '--output', default='refs.bib', help="Output BibTeX file (default: refs.bib)")
    parser.add_argument('--email', help="Contact email for User-Agent header (optional)")
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    extractor = ReferenceExtractor(verbose=args.verbose, email=args.email)
    all_refs = []

    for path_str in args.inputs:
        path = Path(path_str).resolve()
        files = path.glob('**/*.pdf') if path.is_dir() else [path]
        for f in files:
            if f.suffix.lower() != '.pdf':
                if args.verbose: print(f"Skipping non-PDF: {f.name}")
                continue
            if args.verbose: print(f"Processing PDF: {f.name}")
            all_refs.extend(extractor.extract_from_pdf(f))
            time.sleep(0.5)  # be kind to the API

    # Deduplicate by DOI or title
    unique_refs = []
    seen = set()
    for r in all_refs:
        key = (r.doi or r.title or r.generate_cite_key()).strip().lower()
        if key in seen: continue
        seen.add(key)
        unique_refs.append(r)

    with open(args.output, 'w', encoding='utf-8') as f:
        f.write("\n\n".join(r.to_bibtex() for r in unique_refs))

    print(f"\nCompleted! {len(unique_refs)} entries written to {args.output}")
    if not HAS_REQUESTS:
        print("Tip: install requests to enable DOI verification: pip install requests")

if __name__ == "__main__":
    main()