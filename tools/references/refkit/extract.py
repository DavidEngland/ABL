from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, List

from .models import BibEntry
from .normalize import clean_doi

try:
    import fitz  # type: ignore

    HAS_FITZ = True
except ImportError:
    HAS_FITZ = False

try:
    import PyPDF2  # type: ignore

    HAS_PYPDF2 = True
except ImportError:
    HAS_PYPDF2 = False


def _first_year(text: str) -> str:
    match = re.search(r"\b(19|20)\d{2}\b", text)
    return match.group(0) if match else ""


def _extract_doi(text: str) -> str:
    match = re.search(r"(10\.\d{4,9}/[-._;()/:A-Za-z0-9]+)", text)
    return clean_doi(match.group(1)) if match else ""


def _parse_tex_bibitem(text: str, source: str) -> List[BibEntry]:
    entries: List[BibEntry] = []
    pattern = re.compile(
        r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}\s*(.+?)(?=\\bibitem|\\end\{thebibliography\}|$)",
        re.DOTALL,
    )
    for match in pattern.finditer(text):
        key = match.group(1).strip()
        body = re.sub(r"\s+", " ", match.group(2)).strip()
        fields = {
            "title": body[:180],
            "year": _first_year(body),
        }
        doi = _extract_doi(body)
        if doi:
            fields["doi"] = doi
        entries.append(BibEntry(entry_type="article", cite_key=key, fields=fields, source=source))
    return entries


def _pdf_text(path: Path, max_pages: int = 3) -> str:
    if HAS_FITZ:
        doc = fitz.open(path)
        chunks = []
        for i in range(min(max_pages, len(doc))):
            chunks.append(doc[i].get_text())
        doc.close()
        return "\n".join(chunks)

    if HAS_PYPDF2:
        chunks = []
        with open(path, "rb") as handle:
            reader = PyPDF2.PdfReader(handle)
            for i in range(min(max_pages, len(reader.pages))):
                chunks.append(reader.pages[i].extract_text() or "")
        return "\n".join(chunks)

    raise RuntimeError("Install pymupdf or PyPDF2 for PDF extraction")


def extract_from_pdf(path: Path) -> List[BibEntry]:
    text = _pdf_text(path)
    title = ""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    for line in lines[:20]:
        if 15 < len(line) < 220:
            title = line
            break

    doi = _extract_doi(text)
    year = _first_year(text)
    key_seed = re.sub(r"[^a-z0-9]", "", (path.stem.lower() or "ref"))[:24]
    cite_key = f"{key_seed}{year}" if year else key_seed

    fields = {"title": title or path.stem}
    if doi:
        fields["doi"] = doi
    if year:
        fields["year"] = year

    return [BibEntry(entry_type="article", cite_key=cite_key, fields=fields, source=str(path))]


def extract_from_markdown(path: Path) -> List[BibEntry]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    dois = sorted(set(clean_doi(m.group(1)) for m in re.finditer(r"(10\.\d{4,9}/[-._;()/:A-Za-z0-9]+)", text)))
    entries: List[BibEntry] = []
    for idx, doi in enumerate(dois, start=1):
        entries.append(BibEntry(entry_type="article", cite_key=f"mdref{idx}", fields={"doi": doi, "title": doi}, source=str(path)))
    return entries


def extract_from_path(path: Path) -> List[BibEntry]:
    suffix = path.suffix.lower()
    if suffix == ".tex":
        return _parse_tex_bibitem(path.read_text(encoding="utf-8", errors="ignore"), str(path))
    if suffix == ".pdf":
        return extract_from_pdf(path)
    if suffix in {".md", ".markdown", ".html", ".htm"}:
        return extract_from_markdown(path)
    return []


def collect_files(inputs: Iterable[str], recursive: bool, formats: List[str]) -> List[Path]:
    wanted = {f".{f.lower().lstrip('.')}" for f in formats}
    files: List[Path] = []

    for raw in inputs:
        path = Path(raw)
        if path.is_file() and path.suffix.lower() in wanted:
            files.append(path)
            continue
        if path.is_dir():
            for suffix in wanted:
                pattern = f"**/*{suffix}" if recursive else f"*{suffix}"
                files.extend(path.glob(pattern))

    # Stable order for deterministic outputs.
    return sorted(set(files), key=lambda p: str(p))
