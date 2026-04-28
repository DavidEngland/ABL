# UAH Timesheet Quickstart

This workflow keeps payroll happy and minimizes your effort.

## 1) Submit using UAH's official paper timesheet

Use the UAH Payroll paper timesheet form for submission. Keep your local CSV files as your work log and copy totals/notes into the official form.

## 2) Keep your daily log in CSV

Use `social/timesheet_template.csv` to track:

- start/end time
- lunch
- regular hours
- project/task summary

Your identity fields are already prefilled:

- Name: David England
- UAH A-Number: A00084992
- Email: David.England@UAH.Edu

## 3) Auto-generate activity evidence from repo history

From repo root:

```bash
bash social/timesheet_from_git.sh 2026-04-06 2026-04-10
```

This creates:

- `social/timesheet_activity_2026-04-06_to_2026-04-10.csv`
- `social/timesheet_mtime_files_2026-04-06_to_2026-04-10.csv`

Use these as source material to fill your official UAH timesheet quickly and consistently.

## 4) Suggested weekly routine

1. Run the script for the pay-week.
2. Review generated CSV summaries.
3. Enter actual hours on your official UAH form.
4. Attach or archive the CSVs for your personal records.

## 5) Fill official UAH correction CSV automatically

If UAH gives you a correction-form CSV export/template, you can auto-fill it while keeping the official structure intact:

```bash
bash social/fill_uah_correction_csv.sh \
	"/Users/davidengland/Downloads/Timesheet Correction - 01.csv" \
	"/Users/davidengland/Downloads/Timesheet Correction - 01.prefill.csv" \
	"8,8,8,8,8,0,0,8,8,8,8,8,0,0" \
	"Research, writing, and coding tasks in ABL repo"
```

Tip: Set optional metadata once per terminal session if needed:

```bash
export UAH_POSITION="<position number>"
export UAH_TITLE="<job title>"
export UAH_DEPARTMENT="<department>"
export UAH_HOME_LABOR="<home labor>"
```