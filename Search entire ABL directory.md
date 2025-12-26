# Search entire ABL directory
cd /Users/davidengland/Documents/GitHub/ABL
python code/find_mcnider_refs.py . --output mcnider_refs_report.txt

# Or just drafts
python code/find_mcnider_refs.py refs/ --verbose

# Or specific formats
python code/find_mcnider_refs.py . --formats pdf,md,bib