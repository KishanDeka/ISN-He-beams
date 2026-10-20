#!/usr/bin/env bash
# Build report.pdf and clean up the LaTeX build files.
set -euo pipefail
cd "$(dirname "$0")"

# two passes so references and figure numbers resolve
for _ in 1 2; do
    if ! pdflatex -interaction=nonstopmode -halt-on-error report.tex > /dev/null; then
        echo "pdflatex failed, see report.log" >&2
        grep -A3 '^!' report.log >&2 || true
        exit 1
    fi
done

rm -f report.aux report.log report.out report.toc report.synctex.gz
echo "built report.pdf"
