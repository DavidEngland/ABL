# Extract from all PDFs in refs/
python code/refs_to_bib.py refs/*.pdf --output ABL_refs.bib -v

# Or process everything recursively
python code/refs_to_bib.py refs/ drafts/ --recursive --deduplicate --output complete_bibliography.bib