from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

from .models import BibEntry


ENTRY_PATTERN = re.compile(r"@(\w+)\s*\{\s*([^,]+)\s*,(.*?)\n\s*\}\s*", re.DOTALL)
FIELD_PATTERN = re.compile(r"(\w+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|\"[^\"]*\")\s*,?", re.DOTALL)


def _strip_wrappers(value: str) -> str:
    value = value.strip()
    if value.startswith("{") and value.endswith("}"):
        return value[1:-1].strip()
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1].strip()
    return value


def parse_bib_file(path: Path) -> List[BibEntry]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    entries: List[BibEntry] = []

    for match in ENTRY_PATTERN.finditer(text):
        entry_type = match.group(1).strip().lower()
        cite_key = match.group(2).strip()
        body = match.group(3)

        fields: Dict[str, str] = {}
        for fmatch in FIELD_PATTERN.finditer(body):
            key = fmatch.group(1).strip().lower()
            value = _strip_wrappers(fmatch.group(2))
            fields[key] = value

        entries.append(BibEntry(entry_type=entry_type, cite_key=cite_key, fields=fields, source=str(path)))

    return entries


def write_bib_file(path: Path, entries: List[BibEntry]) -> None:
    lines: List[str] = []
    for entry in entries:
        lines.append(f"@{entry.entry_type}{{{entry.cite_key},")
        for key in sorted(entry.fields.keys()):
            value = entry.fields[key].replace("\n", " ").strip()
            lines.append(f"  {key} = {{{value}}},")
        lines.append("}")
        lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
