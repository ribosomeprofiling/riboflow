#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."

hits=$(grep -rnE "^[[:space:]]+(prefix|indexed_prefix|label|suffix|bed_suffix|stats_label|route|args)[[:space:]]*=" \
         modules/ subworkflows/ workflows/ 2>/dev/null | grep -v "def ")

if [ -n "$hits" ]; then
    echo "ERROR: undeclared (session-global) assignment in a process body:" >&2
    echo "$hits" >&2
    echo >&2
    echo "Add 'def', and resolve the name from task.ext.* in the output: block." >&2
    exit 1
fi
echo "OK: no session-global script variables."
