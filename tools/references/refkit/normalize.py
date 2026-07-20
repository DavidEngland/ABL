from __future__ import annotations

import re
from typing import Iterable, List


def normalize_space(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "")).strip()


def normalize_compare(value: str) -> str:
    # Normalize for tolerant field matching.
    value = normalize_space(value).lower()
    value = value.replace("{", "").replace("}", "")
    value = value.replace("–", "-").replace("—", "-")
    value = re.sub(r"[\.,;:]", "", value)
    return value


def clean_doi(value: str) -> str:
    if not value:
        return ""
    value = value.strip()
    value = re.sub(r"^https?://(dx\.)?doi\.org/", "", value, flags=re.IGNORECASE)
    value = re.sub(r"^doi:\s*", "", value, flags=re.IGNORECASE)
    value = value.strip().rstrip(".,;)\"")
    return value


def split_authors(value: str) -> List[str]:
    if not value:
        return []
    return [normalize_space(v) for v in re.split(r"\s+and\s+", value, flags=re.IGNORECASE) if normalize_space(v)]


def author_family_names(authors: Iterable[str]) -> List[str]:
    names = []
    for author in authors:
        if "," in author:
            family = author.split(",", 1)[0]
        else:
            family = author.split()[-1] if author.split() else ""
        family = re.sub(r"[^a-zA-Z]", "", family).lower()
        if family:
            names.append(family)
    return names
