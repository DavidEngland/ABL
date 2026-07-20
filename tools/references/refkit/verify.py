from __future__ import annotations

from dataclasses import asdict
from typing import Dict, List, Tuple

from .crossref import CrossrefClient, crossref_to_bib_fields
from .models import BibEntry, VerificationIssue, VerificationResult
from .normalize import clean_doi, normalize_compare, split_authors, author_family_names


def _authors_match(local_authors: str, remote_authors: str) -> bool:
    local_families = author_family_names(split_authors(local_authors))
    remote_families = author_family_names(split_authors(remote_authors))
    if not local_families or not remote_families:
        return False
    # Strong enough for verification report: first author and majority overlap.
    overlap = set(local_families) & set(remote_families)
    return local_families[0] == remote_families[0] and len(overlap) >= min(2, len(remote_families))


def verify_entries(
    entries: List[BibEntry],
    client: CrossrefClient,
    fill_missing: bool,
) -> Tuple[List[BibEntry], List[VerificationResult]]:
    results: List[VerificationResult] = []

    for entry in entries:
        doi = clean_doi(entry.get("doi"))
        if not doi:
            results.append(VerificationResult(cite_key=entry.cite_key, doi="", status="no-doi"))
            continue

        message = client.get_work(doi)
        if not message:
            results.append(VerificationResult(cite_key=entry.cite_key, doi=doi, status="unresolved"))
            continue

        remote = crossref_to_bib_fields(message)
        issues: List[VerificationIssue] = []
        improved_fields: List[str] = []

        for field, remote_value in remote.items():
            local_value = entry.get(field)
            if not local_value:
                if fill_missing and entry.set_if_missing(field, remote_value):
                    improved_fields.append(field)
                continue

            if field == "author":
                same = _authors_match(local_value, remote_value)
            else:
                same = normalize_compare(local_value) == normalize_compare(remote_value)

            if not same:
                issues.append(
                    VerificationIssue(
                        cite_key=entry.cite_key,
                        doi=doi,
                        field=field,
                        local_value=local_value,
                        remote_value=remote_value,
                        severity="warning",
                    )
                )

        status = "verified"
        if issues:
            status = "conflict"
        elif not remote:
            status = "unresolved"

        results.append(
            VerificationResult(
                cite_key=entry.cite_key,
                doi=doi,
                status=status,
                issues=issues,
                improved_fields=improved_fields,
            )
        )

    return entries, results


def results_to_json(results: List[VerificationResult]) -> List[Dict]:
    payload: List[Dict] = []
    for result in results:
        item = asdict(result)
        item["issues"] = [asdict(issue) for issue in result.issues]
        payload.append(item)
    return payload
