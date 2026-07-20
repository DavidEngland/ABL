from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Optional

from .normalize import clean_doi


class CrossrefClient:
    def __init__(self, email: str = "", cache_path: Optional[Path] = None, timeout: int = 20):
        self.email = email.strip()
        self.timeout = timeout
        self.cache_path = cache_path
        self._cache: Dict[str, Dict] = {}
        self._loaded = False

    def _load_cache(self) -> None:
        if self._loaded:
            return
        self._loaded = True
        if not self.cache_path or not self.cache_path.exists():
            return
        try:
            self._cache = json.loads(self.cache_path.read_text(encoding="utf-8"))
        except Exception:
            self._cache = {}

    def _save_cache(self) -> None:
        if not self.cache_path:
            return
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache_path.write_text(json.dumps(self._cache, indent=2, sort_keys=True), encoding="utf-8")

    def get_work(self, doi: str) -> Optional[Dict]:
        doi = clean_doi(doi)
        if not doi:
            return None
        self._load_cache()
        if doi in self._cache:
            return self._cache[doi]

        encoded = urllib.parse.quote(doi, safe="")
        url = f"https://api.crossref.org/works/{encoded}"
        ua = "ABL-refs-toolkit/0.1"
        if self.email:
            ua += f" (mailto:{self.email})"
        req = urllib.request.Request(url, headers={"User-Agent": ua, "Accept": "application/json"})

        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as response:
                payload = json.loads(response.read().decode("utf-8", errors="ignore"))
                message = payload.get("message", {})
        except Exception:
            message = None

        if message:
            self._cache[doi] = message
            self._save_cache()
        return message


def crossref_to_bib_fields(message: Dict) -> Dict[str, str]:
    fields: Dict[str, str] = {}

    title = (message.get("title") or [""])[0]
    if title:
        fields["title"] = title

    container = (message.get("container-title") or [""])[0]
    if container:
        fields["journal"] = container

    authors = []
    for a in message.get("author", []) or []:
        family = (a.get("family") or "").strip()
        given = (a.get("given") or "").strip()
        if family and given:
            authors.append(f"{family}, {given}")
        elif family:
            authors.append(family)
    if authors:
        fields["author"] = " and ".join(authors)

    year_parts = (((message.get("issued") or {}).get("date-parts") or [[None]])[0])
    if year_parts and year_parts[0]:
        fields["year"] = str(year_parts[0])

    for src, dst in [("volume", "volume"), ("issue", "number"), ("page", "pages"), ("publisher", "publisher")]:
        value = message.get(src)
        if value:
            fields[dst] = str(value)

    doi = message.get("DOI")
    if doi:
        fields["doi"] = doi

    return fields
