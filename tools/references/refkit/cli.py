from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List

from .bibtex_io import parse_bib_file, write_bib_file
from .crossref import CrossrefClient
from .extract import collect_files, extract_from_path
from .models import BibEntry
from .reporting import write_markdown_report
from .tex_sync import compare_tex_to_bib
from .verify import results_to_json, verify_entries


def _write_json(path: Path, payload: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def cmd_extract(args: argparse.Namespace) -> int:
    files = collect_files(args.inputs, recursive=args.recursive, formats=args.formats.split(","))
    entries: List[BibEntry] = []

    for path in files:
        if path.suffix.lower() == ".bib":
            entries.extend(parse_bib_file(path))
        else:
            entries.extend(extract_from_path(path))

    write_bib_file(Path(args.output_bib), entries)
    print(f"Wrote {len(entries)} entries to {args.output_bib}")
    return 0


def _load_entries_from_paths(paths: List[str]) -> List[BibEntry]:
    entries: List[BibEntry] = []
    for raw in paths:
        path = Path(raw)
        if path.suffix.lower() != ".bib":
            continue
        entries.extend(parse_bib_file(path))
    return entries


def cmd_verify_doi(args: argparse.Namespace) -> int:
    entries = _load_entries_from_paths(args.bib)
    client = CrossrefClient(email=args.email or "", cache_path=Path(args.cache))
    _, results = verify_entries(entries, client, fill_missing=False)

    write_markdown_report(Path(args.report_md), results, title="DOI Verification Report")
    _write_json(Path(args.report_json), results_to_json(results))
    print(f"Verified {len(results)} entries; report: {args.report_md}")
    return 0


def cmd_improve(args: argparse.Namespace) -> int:
    entries = _load_entries_from_paths(args.bib)
    client = CrossrefClient(email=args.email or "", cache_path=Path(args.cache))
    improved_entries, results = verify_entries(entries, client, fill_missing=True)

    write_bib_file(Path(args.output_bib), improved_entries)
    write_markdown_report(Path(args.report_md), results, title="Bibliography Improvement Report")
    _write_json(Path(args.report_json), results_to_json(results))
    print(f"Improved {len(improved_entries)} entries; output: {args.output_bib}")
    return 0


def cmd_sync_paper(args: argparse.Namespace) -> int:
    entries = _load_entries_from_paths(args.bib)
    tex_files = [Path(p) for p in args.tex]
    result = compare_tex_to_bib(tex_files, entries)

    md_lines = [
        "# Paper Citation Sync Report",
        "",
        f"- Missing in BibTeX: {len(result['missing_in_bib'])}",
        f"- Unused in BibTeX: {len(result['unused_in_bib'])}",
        "",
        "## Missing in BibTeX",
        "",
    ]
    if result["missing_in_bib"]:
        md_lines.extend([f"- {k}" for k in result["missing_in_bib"]])
    else:
        md_lines.append("None")

    md_lines.extend(["", "## Unused in BibTeX", ""])
    if result["unused_in_bib"]:
        md_lines.extend([f"- {k}" for k in result["unused_in_bib"]])
    else:
        md_lines.append("None")

    report_path = Path(args.report_md)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(md_lines).rstrip() + "\n", encoding="utf-8")
    _write_json(Path(args.report_json), result)

    print(f"Synced {len(tex_files)} TeX files against {len(entries)} BibTeX entries")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="References toolkit: extract, verify, improve, and sync BibTeX")
    sub = parser.add_subparsers(dest="command", required=True)

    p_extract = sub.add_parser("extract", help="Extract references from PDF/TeX/Bib/Markdown into BibTeX")
    p_extract.add_argument("inputs", nargs="+", help="Input files or directories")
    p_extract.add_argument("--formats", default="pdf,tex,bib,md,html", help="Comma-separated formats to include")
    p_extract.add_argument("--recursive", action="store_true", help="Recursively scan directories")
    p_extract.add_argument("--output-bib", required=True, help="Output BibTeX file")
    p_extract.set_defaults(func=cmd_extract)

    p_verify = sub.add_parser("verify-doi", help="Verify BibTeX fields against Crossref DOI metadata")
    p_verify.add_argument("--bib", nargs="+", required=True, help="Input .bib file(s)")
    p_verify.add_argument("--email", default="", help="Contact email for Crossref user-agent")
    p_verify.add_argument("--cache", default="tools/references/output/crossref_cache.json", help="Crossref cache JSON path")
    p_verify.add_argument("--report-md", default="tools/references/output/doi_verification_report.md")
    p_verify.add_argument("--report-json", default="tools/references/output/doi_verification_report.json")
    p_verify.set_defaults(func=cmd_verify_doi)

    p_improve = sub.add_parser("improve", help="Balanced improvement: fill missing fields, report conflicts")
    p_improve.add_argument("--bib", nargs="+", required=True, help="Input .bib file(s)")
    p_improve.add_argument("--email", default="", help="Contact email for Crossref user-agent")
    p_improve.add_argument("--cache", default="tools/references/output/crossref_cache.json", help="Crossref cache JSON path")
    p_improve.add_argument("--output-bib", required=True, help="Improved .bib output file")
    p_improve.add_argument("--report-md", default="tools/references/output/improve_report.md")
    p_improve.add_argument("--report-json", default="tools/references/output/improve_report.json")
    p_improve.set_defaults(func=cmd_improve)

    p_sync = sub.add_parser("sync-paper", help="Compare TeX citation keys to BibTeX keys")
    p_sync.add_argument("--tex", nargs="+", required=True, help="Input .tex files")
    p_sync.add_argument("--bib", nargs="+", required=True, help="Input .bib file(s)")
    p_sync.add_argument("--report-md", default="tools/references/output/sync_paper_report.md")
    p_sync.add_argument("--report-json", default="tools/references/output/sync_paper_report.json")
    p_sync.set_defaults(func=cmd_sync_paper)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
