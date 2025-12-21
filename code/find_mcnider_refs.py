#!/usr/bin/env python3
"""Find all references to McNider & England (1995) stability functions paper.

Searches across PDF, MD, DOCX, TEX, TXT, BIB, and other formats for citations
to the paper "Stability functions based upon shear functions" in Boundary-Layer
Meteorology.

Usage:
    python find_mcnider_refs.py [directory] [--verbose]
    python find_mcnider_refs.py . --formats pdf,md,tex,bib
"""

import re
import argparse
from pathlib import Path
from typing import List, Dict, Tuple
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

# Try to import DOCX library
try:
    import docx
    HAS_DOCX = True
except ImportError:
    HAS_DOCX = False


class McNiderRefFinder:
    """Search for references to McNider & England stability functions paper."""
    
    # Search patterns for the paper
    PATTERNS = [
        # Direct title matches
        r'Stability\s+functions?\s+based\s+upon\s+shear\s+functions?',
        
        # Author-year citations
        r'McNider\s*(?:&|and)\s*England\s*[\(\[]?\s*1995\s*[\)\]]?',
        r'England\s*(?:&|and)\s*McNider\s*[\(\[]?\s*1995\s*[\)\]]?',
        
        # BibTeX key patterns
        r'mcnider(?:1995|_1995|1995[a-z]?)',
        r'england(?:1995|_1995|1995[a-z]?)',
        
        # Journal citation patterns
        r'Boundary[- ]Layer\s+Meteorol(?:ogy|\.)?\s*,?\s*7[234]',  # Volume 72-74
        
        # Combined patterns
        r'McNider.*England.*1995.*(?:stability|shear)',
        r'England.*McNider.*1995.*(?:stability|shear)',
    ]
    
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.results: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    
    def search_text(self, text: str, file_path: Path) -> List[Tuple[int, str]]:
        """Search text for McNider paper references.
        
        Returns:
            List of (line_number, matched_context) tuples
        """
        matches = []
        lines = text.split('\n')
        
        for line_num, line in enumerate(lines, 1):
            for pattern in self.PATTERNS:
                if re.search(pattern, line, re.IGNORECASE):
                    # Extract context (40 chars before and after match)
                    match = re.search(pattern, line, re.IGNORECASE)
                    if match:
                        start = max(0, match.start() - 40)
                        end = min(len(line), match.end() + 40)
                        context = line[start:end].strip()
                        matches.append((line_num, context))
                        break  # Only record first pattern match per line
        
        return matches
    
    def search_pdf(self, pdf_path: Path) -> List[Tuple[int, str]]:
        """Search PDF file for references."""
        if self.verbose:
            print(f"  Searching PDF: {pdf_path.name}")
        
        try:
            if HAS_FITZ:
                doc = fitz.open(pdf_path)
                text = "\n".join(page.get_text() for page in doc)
                doc.close()
            elif HAS_PYPDF2:
                with open(pdf_path, 'rb') as f:
                    reader = PyPDF2.PdfReader(f)
                    text = "\n".join(page.extract_text() for page in reader.pages)
            else:
                return []
            
            return self.search_text(text, pdf_path)
            
        except Exception as e:
            if self.verbose:
                print(f"    Error reading PDF: {e}")
            return []
    
    def search_docx(self, docx_path: Path) -> List[Tuple[int, str]]:
        """Search DOCX file for references."""
        if not HAS_DOCX:
            return []
        
        if self.verbose:
            print(f"  Searching DOCX: {docx_path.name}")
        
        try:
            doc = docx.Document(docx_path)
            text = "\n".join(para.text for para in doc.paragraphs)
            return self.search_text(text, docx_path)
            
        except Exception as e:
            if self.verbose:
                print(f"    Error reading DOCX: {e}")
            return []
    
    def search_text_file(self, file_path: Path) -> List[Tuple[int, str]]:
        """Search text-based file (MD, TEX, TXT, BIB, etc.)."""
        if self.verbose:
            print(f"  Searching text: {file_path.name}")
        
        try:
            text = file_path.read_text(encoding='utf-8', errors='ignore')
            return self.search_text(text, file_path)
            
        except Exception as e:
            if self.verbose:
                print(f"    Error reading file: {e}")
            return []
    
    def search_file(self, file_path: Path) -> List[Tuple[int, str]]:
        """Search a single file based on its extension."""
        suffix = file_path.suffix.lower()
        
        if suffix == '.pdf':
            return self.search_pdf(file_path)
        elif suffix == '.docx':
            return self.search_docx(file_path)
        elif suffix in ['.md', '.tex', '.txt', '.bib', '.rst', '.latex', '.markdown']:
            return self.search_text_file(file_path)
        else:
            # Try as text file anyway
            return self.search_text_file(file_path)
    
    def search_directory(self, directory: Path, formats: List[str], recursive: bool = True):
        """Search directory for McNider paper references.
        
        Args:
            directory: Root directory to search
            formats: List of file extensions to search (e.g., ['pdf', 'md'])
            recursive: Whether to search subdirectories
        """
        print(f"Searching in: {directory}")
        print(f"Formats: {', '.join(formats)}")
        print("-" * 70)
        
        # Build file pattern
        pattern = "**/*" if recursive else "*"
        
        for ext in formats:
            ext = ext.strip().lower()
            if not ext.startswith('.'):
                ext = f'.{ext}'
            
            for file_path in directory.glob(f"{pattern}{ext}"):
                if not file_path.is_file():
                    continue
                
                matches = self.search_file(file_path)
                
                if matches:
                    self.results[str(file_path)] = matches
                    print(f"\n✓ Found {len(matches)} reference(s) in: {file_path}")
                    
                    for line_num, context in matches:
                        print(f"  Line {line_num}: ...{context}...")
    
    def generate_report(self, output_file: Path = None):
        """Generate summary report of findings."""
        report_lines = [
            "=" * 70,
            "McNider & England (1995) Reference Search Report",
            "=" * 70,
            "",
            f"Total files with references: {len(self.results)}",
            f"Total references found: {sum(len(matches) for matches in self.results.values())}",
            "",
            "Files by type:",
        ]
        
        # Count by file type
        type_counts = defaultdict(int)
        for file_path in self.results.keys():
            suffix = Path(file_path).suffix.lower()
            type_counts[suffix] += 1
        
        for suffix, count in sorted(type_counts.items()):
            report_lines.append(f"  {suffix or '(no extension)'}: {count} file(s)")
        
        report_lines.extend(["", "=" * 70, "Detailed Results", "=" * 70])
        
        # Detailed results by file
        for file_path, matches in sorted(self.results.items()):
            report_lines.append(f"\n{file_path}")
            report_lines.append("-" * 70)
            
            for line_num, context in matches:
                report_lines.append(f"  Line {line_num}: {context}")
        
        report = "\n".join(report_lines)
        
        # Print to console
        print("\n" + report)
        
        # Write to file if requested
        if output_file:
            output_file.write_text(report, encoding='utf-8')
            print(f"\nReport saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Find references to McNider & England (1995) stability functions paper'
    )
    parser.add_argument(
        'directory',
        nargs='?',
        default='.',
        help='Directory to search (default: current directory)'
    )
    parser.add_argument(
        '-f', '--formats',
        default='pdf,md,tex,bib,txt,docx,rst',
        help='Comma-separated list of file extensions to search (default: pdf,md,tex,bib,txt,docx,rst)'
    )
    parser.add_argument(
        '-r', '--recursive',
        action='store_true',
        default=True,
        help='Search subdirectories recursively (default: True)'
    )
    parser.add_argument(
        '--no-recursive',
        dest='recursive',
        action='store_false',
        help='Do not search subdirectories'
    )
    parser.add_argument(
        '-o', '--output',
        help='Save report to file'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Check dependencies
    if 'pdf' in args.formats.lower() and not (HAS_FITZ or HAS_PYPDF2):
        print("Warning: No PDF library found. Install one of:")
        print("  pip install pymupdf  (recommended)")
        print("  pip install PyPDF2   (fallback)")
        print("\nSkipping PDF files...")
        args.formats = args.formats.replace('pdf', '').replace(',,', ',')
    
    if 'docx' in args.formats.lower() and not HAS_DOCX:
        print("Warning: python-docx not found. Install with:")
        print("  pip install python-docx")
        print("\nSkipping DOCX files...")
        args.formats = args.formats.replace('docx', '').replace(',,', ',')
    
    # Parse directory
    directory = Path(args.directory).resolve()
    if not directory.exists():
        print(f"Error: Directory not found: {directory}")
        return 1
    
    # Parse formats
    formats = [f.strip() for f in args.formats.split(',') if f.strip()]
    if not formats:
        print("Error: No valid file formats specified")
        return 1
    
    # Run search
    finder = McNiderRefFinder(verbose=args.verbose)
    finder.search_directory(directory, formats, recursive=args.recursive)
    
    # Generate report
    output_file = Path(args.output) if args.output else None
    finder.generate_report(output_file)
    
    # Return status
    return 0 if finder.results else 1


if __name__ == '__main__':
    exit(main())
