#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: social/timesheet_from_git.sh YYYY-MM-DD YYYY-MM-DD

Generates two CSV files in social/ for the date range:
1) timesheet_activity_<start>_to_<end>.csv      (daily git activity summary)
2) timesheet_mtime_files_<start>_to_<end>.csv   (files with mtimes in range)

Example:
  social/timesheet_from_git.sh 2026-04-06 2026-04-10
EOF
}

if [[ $# -ne 2 ]]; then
  usage
  exit 1
fi

START_DATE="$1"
END_DATE="$2"

if ! date -j -f "%Y-%m-%d" "$START_DATE" "+%Y-%m-%d" >/dev/null 2>&1; then
  echo "Invalid start date: $START_DATE (expected YYYY-MM-DD)" >&2
  exit 1
fi

if ! date -j -f "%Y-%m-%d" "$END_DATE" "+%Y-%m-%d" >/dev/null 2>&1; then
  echo "Invalid end date: $END_DATE (expected YYYY-MM-DD)" >&2
  exit 1
fi

if [[ "$START_DATE" > "$END_DATE" ]]; then
  echo "Start date must be <= end date" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

ACTIVITY_CSV="$SCRIPT_DIR/timesheet_activity_${START_DATE}_to_${END_DATE}.csv"
MTIME_CSV="$SCRIPT_DIR/timesheet_mtime_files_${START_DATE}_to_${END_DATE}.csv"

# Include the full end day by using end + 1 day as an exclusive upper bound.
END_PLUS_ONE="$(date -j -v+1d -f "%Y-%m-%d" "$END_DATE" "+%Y-%m-%d")"

{
  echo 'Date,Commit Count,Commit Subjects,Files Changed,Project/Task,Notes'

  current_date="$START_DATE"
  while [[ "$current_date" < "$END_PLUS_ONE" ]]; do
    commit_count="$(git -C "$REPO_ROOT" rev-list --count --since "$current_date 00:00:00" --until "$current_date 23:59:59" HEAD)"

    subjects="$(git -C "$REPO_ROOT" log --since "$current_date 00:00:00" --until "$current_date 23:59:59" --pretty=format:%s | paste -sd '; ' -)"
    files_changed="$(git -C "$REPO_ROOT" log --since "$current_date 00:00:00" --until "$current_date 23:59:59" --name-only --pretty=format: | sed '/^$/d' | sort -u | wc -l | tr -d ' ')"

    if [[ -z "$subjects" ]]; then
      subjects="No commit recorded"
    fi

    # Basic CSV escaping for quoted text columns.
    subjects_escaped="${subjects//\"/\"\"}"

    echo "$current_date,$commit_count,\"$subjects_escaped\",$files_changed,Research/Writing/Code,Fill in actual hours + key outcomes"

    current_date="$(date -j -v+1d -f "%Y-%m-%d" "$current_date" "+%Y-%m-%d")"
  done
} > "$ACTIVITY_CSV"

{
  echo 'Modified Date,File Path'
  find "$REPO_ROOT" -type f -not -path '*/.git/*' -newermt "$START_DATE" ! -newermt "$END_PLUS_ONE" \
    -exec stat -f '%Sm,%N' -t '%Y-%m-%d' {} \; \
    | sed "s#$REPO_ROOT/##" \
    | sort
} > "$MTIME_CSV"

echo "Wrote: $ACTIVITY_CSV"
echo "Wrote: $MTIME_CSV"
