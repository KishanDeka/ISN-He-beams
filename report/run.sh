#!/usr/bin/env bash
# Build report.pdf and clean up the LaTeX build files.
set -euo pipefail
cd "$(dirname "$0")"

# rerun until references and equation numbers are stable (at most 4 passes)
for _ in 1 2 3 4; do
    if ! pdflatex -interaction=nonstopmode -halt-on-error report.tex > /dev/null; then
        echo "pdflatex failed, see report.log" >&2
        grep -A3 '^!' report.log >&2 || true
        exit 1
    fi
    grep -q "Rerun to get" report.log || break
done
grep -E "Overfull|undefined" report.log >&2 || true

rm -f *.aux *.log *.out *.toc *.synctex.gz
echo "built report.pdf"
