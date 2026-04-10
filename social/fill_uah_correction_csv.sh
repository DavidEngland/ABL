#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  social/fill_uah_correction_csv.sh INPUT_CSV OUTPUT_CSV HOURS_14 [COMMENT]

Arguments:
  INPUT_CSV   Path to official UAH correction CSV export/template.
  OUTPUT_CSV  Path for the filled output CSV.
  HOURS_14    Comma-separated 14 daily regular-hour values.
              Example: 8,8,8,8,8,0,0,8,8,8,8,8,0,0
  COMMENT     Optional comment text for the Comment: row.

Environment variables (optional):
  UAH_NAME         Default: David England
  UAH_A_NUMBER     Default: A00084992
  UAH_POSITION     Default: (blank)
  UAH_TITLE        Default: (blank)
  UAH_DEPARTMENT   Default: (blank)
  UAH_HOME_LABOR   Default: (blank)
EOF
}

if [[ $# -lt 3 ]]; then
  usage
  exit 1
fi

INPUT_CSV="$1"
OUTPUT_CSV="$2"
HOURS_14="$3"
COMMENT="${4:-}"

UAH_NAME="${UAH_NAME:-David England}"
UAH_A_NUMBER="${UAH_A_NUMBER:-A00084992}"
UAH_POSITION="${UAH_POSITION:-}"
UAH_TITLE="${UAH_TITLE:-}"
UAH_DEPARTMENT="${UAH_DEPARTMENT:-}"
UAH_HOME_LABOR="${UAH_HOME_LABOR:-}"

if [[ ! -f "$INPUT_CSV" ]]; then
  echo "Input CSV not found: $INPUT_CSV" >&2
  exit 1
fi

awk -F',' \
  -v OFS=',' \
  -v name="$UAH_NAME" \
  -v anum="$UAH_A_NUMBER" \
  -v position="$UAH_POSITION" \
  -v title="$UAH_TITLE" \
  -v department="$UAH_DEPARTMENT" \
  -v home_labor="$UAH_HOME_LABOR" \
  -v hours_csv="$HOURS_14" \
  -v comment="$COMMENT" '
BEGIN {
  n = split(hours_csv, h, ",")
  if (n != 14) {
    print "HOURS_14 must contain exactly 14 comma-separated values" > "/dev/stderr"
    exit 2
  }

  total = 0
  week1 = 0
  week2 = 0
  for (i = 1; i <= 14; i++) {
    val = h[i] + 0
    h[i] = sprintf("%.2f", val)
    total += val
    if (i <= 7) {
      week1 += val
    } else {
      week2 += val
    }
  }

  total_s = sprintf("%.2f", total)
  week1_s = sprintf("%.2f", week1)
  week2_s = sprintf("%.2f", week2)
}

# Name/A#/Position/Title line.
NR == 6 {
  $2 = name
  $6 = anum
  if (position != "") {
    $10 = position
  }
  if (title != "") {
    $14 = title
  }
}

# Department / Home Labor lines.
NR == 8 {
  if (department != "") {
    $2 = department
  }
}
NR == 10 {
  if (home_labor != "") {
    $2 = home_labor
  }
}

# Regular Hours detail row (14 values + total).
NR == 16 {
  for (i = 1; i <= 14; i++) {
    $i + 2 = h[i]
  }
  $17 = total_s
}

# Total Hours row mirrors Regular Hours when leave is 0.
NR == 34 {
  for (i = 1; i <= 14; i++) {
    $i + 2 = h[i]
  }
  $17 = total_s
}

# Week 1 / Week 2 totals row.
NR == 35 {
  $8 = week1_s
  $15 = week2_s
}

# Comment row.
NR == 38 {
  if (comment != "") {
    $2 = comment
  }
}

{ print }
' "$INPUT_CSV" > "$OUTPUT_CSV"

echo "Wrote: $OUTPUT_CSV"
