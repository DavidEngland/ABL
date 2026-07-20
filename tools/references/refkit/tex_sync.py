from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Set

from .models import BibEntry


CITE_PATTERN = re.compile(r"\\cite[t|p|alp|year|author]*\*?(?:\[[^\]]*\])?(?:\[[^\]]*\])?\{([^}]+)\}")


def extract_citation_keys(tex_path: Path) -> Set[str]:
    text = tex_path.read_text(encoding="utf-8", errors="ignore")
    keys: Set[str] = set()

    for match in CITE_PATTERN.finditer(text):
        raw_keys = match.group(1)
        for key in raw_keys.split(","):
            clean = key.strip()
            if clean:
                keys.add(clean)

    return keys


def compare_tex_to_bib(tex_files: List[Path], entries: List[BibEntry]) -> Dict[str, List[str]]:
    cited: Set[str] = set()
    for path in tex_files:
        cited |= extract_citation_keys(path)

    bib_keys = {e.cite_key for e in entries}

    missing_in_bib = sorted(cited - bib_keys)
    unused_in_bib = sorted(bib_keys - cited)

    return {
        "missing_in_bib": missing_in_bib,
        "unused_in_bib": unused_in_bib,
        "cited_total": sorted(cited),
        "bib_total": sorted(bib_keys),
    }
