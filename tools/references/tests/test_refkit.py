import tempfile
import unittest
from pathlib import Path

from refkit.normalize import clean_doi
from refkit.bibtex_io import parse_bib_file
from refkit.tex_sync import compare_tex_to_bib


class RefkitTests(unittest.TestCase):
    def test_clean_doi(self):
        self.assertEqual(clean_doi("https://doi.org/10.1000/xyz.1."), "10.1000/xyz.1")
        self.assertEqual(clean_doi("DOI: 10.5555/abc"), "10.5555/abc")

    def test_tex_bib_sync(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tex = root / "paper.tex"
            bib = root / "refs.bib"

            tex.write_text("""
            \\cite{known2020,missing2022}
            """, encoding="utf-8")
            bib.write_text("""
            @article{known2020,
              title = {Known},
              year = {2020},
            }
            @article{unused2019,
              title = {Unused},
              year = {2019},
            }
            """, encoding="utf-8")

            entries = parse_bib_file(bib)
            result = compare_tex_to_bib([tex], entries)

            self.assertEqual(result["missing_in_bib"], ["missing2022"])
            self.assertEqual(result["unused_in_bib"], ["unused2019"])


if __name__ == "__main__":
    unittest.main()
