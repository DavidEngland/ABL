from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class BibEntry:
    entry_type: str
    cite_key: str
    fields: Dict[str, str] = field(default_factory=dict)
    source: Optional[str] = None

    def get(self, key: str) -> str:
        return self.fields.get(key.lower(), "")

    def set_if_missing(self, key: str, value: str) -> bool:
        key = key.lower()
        if not value:
            return False
        if self.fields.get(key, "").strip():
            return False
        self.fields[key] = value.strip()
        return True


@dataclass
class VerificationIssue:
    cite_key: str
    doi: str
    field: str
    local_value: str
    remote_value: str
    severity: str  # info | warning


@dataclass
class VerificationResult:
    cite_key: str
    doi: str
    status: str  # verified | conflict | unresolved | no-doi
    issues: List[VerificationIssue] = field(default_factory=list)
    improved_fields: List[str] = field(default_factory=list)


FIELD_MAP = {
    "title": "title",
    "journal": "container-title",
    "author": "author",
    "year": "issued",
    "volume": "volume",
    "number": "issue",
    "pages": "page",
    "publisher": "publisher",
}
