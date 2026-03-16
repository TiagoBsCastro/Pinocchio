#!/usr/bin/env bash
set -euo pipefail

# Build script for A&A LaTeX paper
# - Compiles draft.tex into build/draft.pdf by default
# - Uses latexmk if available, otherwise falls back to pdflatex+bibtex
# - By default, cleans auxiliary files after a successful build (keeps the PDF)
# - Provides clean and dependency-check options

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

FILE="draft.tex"
OUTDIR="build"
OPEN_PDF=false
DO_CLEAN=false
DO_DEPS=false
QUIET=false
KEEP_AUX=false

usage() {
  cat <<'EOF'
Usage: ./build.sh [options]

Options:
  -f, --file FILE       TeX file to compile (default: draft.tex)
  -o, --outdir DIR      Output directory (default: build)
  -c, --clean           Clean build artifacts and exit
      --open            Open resulting PDF with xdg-open
      --deps            Check LaTeX toolchain dependencies and exit
      --keep-aux        Do not remove auxiliary files after build
  -q, --quiet           Less verbose compiler output (where supported)
  -h, --help            Show this help and exit

Examples:
  ./build.sh
  ./build.sh -f draft.tex --open
  ./build.sh --clean
EOF
}

has() { command -v "$1" >/dev/null 2>&1; }

einfo() { echo "[INFO] $*"; }
ewarn() { echo "[WARN] $*" >&2; }
eerr()  { echo "[ERR ] $*" >&2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f|--file) FILE="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -c|--clean) DO_CLEAN=true; shift;;
    --open) OPEN_PDF=true; shift;;
    --deps) DO_DEPS=true; shift;;
    --keep-aux) KEEP_AUX=true; shift;;
    -q|--quiet) QUIET=true; shift;;
    -h|--help) usage; exit 0;;
    *) eerr "Unknown option: $1"; usage; exit 2;;
  esac
done

BASENAME="$(basename "$FILE" .tex)"
PDF_OUT="$OUTDIR/$BASENAME.pdf"

check_deps() {
  local missing=()
  for tool in latexmk pdflatex bibtex xdg-open; do
    if ! has "$tool"; then
      missing+=("$tool")
    fi
  done
  if [[ ${#missing[@]} -eq 0 ]]; then
    einfo "All optional tools found: latexmk, pdflatex, bibtex, xdg-open"
    return 0
  else
    ewarn "Missing tools (script will fallback when possible): ${missing[*]}"
    return 1
  fi
}

clean() {
  if has latexmk; then
    einfo "Cleaning via latexmk (-C) in outdir '$OUTDIR'"
    latexmk -C -outdir="$OUTDIR" >/dev/null 2>&1 || true
  fi
  if [[ -d "$OUTDIR" ]]; then
    einfo "Removing directory '$OUTDIR'"
    rm -rf "$OUTDIR"
  fi
  einfo "Clean complete"
}

clean_aux() {
  # Remove auxiliary files but keep final PDFs
  if has latexmk; then
    einfo "Removing auxiliary files via latexmk (-c) in '$OUTDIR'"
    latexmk -c -outdir="$OUTDIR" >/dev/null 2>&1 || true
  else
    if [[ -d "$OUTDIR" ]]; then
      einfo "Removing auxiliary files from '$OUTDIR' (keeping PDFs)"
      rm -f "$OUTDIR"/*.aux "$OUTDIR"/*.bbl "$OUTDIR"/*.blg "$OUTDIR"/*.bcf \
            "$OUTDIR"/*.brf "$OUTDIR"/*.idx "$OUTDIR"/*.ilg "$OUTDIR"/*.ind \
            "$OUTDIR"/*.lof "$OUTDIR"/*.log "$OUTDIR"/*.lot "$OUTDIR"/*.nav \
            "$OUTDIR"/*.out "$OUTDIR"/*.snm "$OUTDIR"/*.synctex.gz "$OUTDIR"/*.toc \
            "$OUTDIR"/*.run.xml "$OUTDIR"/*.fdb_latexmk "$OUTDIR"/*.fls || true
    fi
  fi
}

compile_latexmk() {
  local args=(-pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir="$OUTDIR" -f)
  [[ "$QUIET" == true ]] && args+=(-silent)
  einfo "latexmk ${args[*]} $FILE"
  latexmk "${args[@]}" "$FILE"
}

compile_fallback() {
  mkdir -p "$OUTDIR"
  local pdargs=(-interaction=nonstopmode -halt-on-error -file-line-error -output-directory "$OUTDIR")
  einfo "pdflatex ${pdargs[*]} $FILE"
  pdflatex "${pdargs[@]}" "$FILE" >/dev/null

  # Decide whether to run bibtex: check for \bibliography or \addbibresource
  local do_bib=0
  if grep -q '\\bibliography{' "$FILE" || grep -q '\\addbibresource{' "$FILE"; then
    do_bib=1
  fi

  if [[ $do_bib -eq 1 ]] && has bibtex; then
    if [[ -f "$OUTDIR/$BASENAME.aux" ]]; then
      ( cd "$OUTDIR" && einfo "bibtex $BASENAME" && bibtex "$BASENAME" >/dev/null || true )
    fi
  else
    [[ $do_bib -eq 1 ]] && ewarn "BibTeX requested by source but 'bibtex' not found; continuing without bibliography"
  fi

  # Two more passes for references
  pdflatex "${pdargs[@]}" "$FILE" >/dev/null
  pdflatex "${pdargs[@]}" "$FILE" >/dev/null
}

main() {
  if [[ "$DO_DEPS" == true ]]; then
    check_deps; exit 0
  fi

  if [[ "$DO_CLEAN" == true ]]; then
    clean; exit 0
  fi

  if [[ ! -f "$FILE" ]]; then
    eerr "TeX file not found: $FILE"
    exit 3
  fi

  if has latexmk; then
    compile_latexmk
  else
    ewarn "'latexmk' not found, falling back to pdflatex+bibtex"
    compile_fallback
  fi

  if [[ -f "$PDF_OUT" ]]; then
    einfo "PDF generated: $PDF_OUT"
    if [[ "$OPEN_PDF" == true ]]; then
      if has xdg-open; then
        einfo "Opening PDF with xdg-open"
        xdg-open "$PDF_OUT" >/dev/null 2>&1 || true
      else
        ewarn "'xdg-open' not found; not opening PDF"
      fi
    fi
    # Post-build cleanup of auxiliary files unless user asked to keep them
    if [[ "$KEEP_AUX" == false ]]; then
      clean_aux
    else
      einfo "Keeping auxiliary files as requested (--keep-aux)"
    fi
  else
    ewarn "PDF not found at expected location: $PDF_OUT"
    exit 4
  fi
}

main "$@"
