from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from .models import VerificationResult


def summarize(results: List[VerificationResult]) -> Dict[str, int]:
    stats = {"verified": 0, "conflict": 0, "unresolved": 0, "no-doi": 0}
    for result in results:
        stats[result.status] = stats.get(result.status, 0) + 1
    return stats


def write_markdown_report(path: Path, results: List[VerificationResult], title: str) -> None:
    stats = summarize(results)
    lines: List[str] = []

    lines.append(f"# {title}")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Verified: {stats.get('verified', 0)}")
    lines.append(f"- Conflicts: {stats.get('conflict', 0)}")
    lines.append(f"- Unresolved DOI lookups: {stats.get('unresolved', 0)}")
    lines.append(f"- Entries without DOI: {stats.get('no-doi', 0)}")
    lines.append("")

    conflict_rows = [r for r in results if r.status == "conflict"]
    unresolved_rows = [r for r in results if r.status in {"unresolved", "no-doi"}]

    lines.append("## Conflicts")
    lines.append("")
    if not conflict_rows:
        lines.append("No DOI metadata conflicts found.")
    else:
        for row in conflict_rows:
            lines.append(f"### {row.cite_key} ({row.doi})")
            lines.append("")
            for issue in row.issues:
                lines.append(f"- Field `{issue.field}`")
                lines.append(f"  - local: {issue.local_value}")
                lines.append(f"  - crossref: {issue.remote_value}")
            if row.improved_fields:
                lines.append(f"- Filled missing fields: {', '.join(row.improved_fields)}")
            lines.append("")

    lines.append("## Unresolved")
    lines.append("")
    if not unresolved_rows:
        lines.append("No unresolved entries.")
    else:
        for row in unresolved_rows:
            if row.status == "no-doi":
                lines.append(f"- {row.cite_key}: missing DOI")
            else:
                lines.append(f"- {row.cite_key}: DOI lookup failed for {row.doi}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
