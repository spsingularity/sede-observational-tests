#!/usr/bin/env bash
# Build the JCAP-class PDF (jcappub.sty) for Paper IV from SEDE_observational_tests.md.
#   makedoc (title/abstract -> YAML) -> strip the citeproc refs placeholder
#   -> pandoc --natbib (emits \citep) with template_jcap.tex -> xelatex + bibtex.
set -e
cd "$(dirname "$0")"
mkdir -p tex

python3 tools/makedoc.py SEDE_observational_tests.md .SEDE.obstests.md
trap 'rm -f .SEDE.obstests.md' EXIT
# drop the citeproc "## References / ::: {#refs}" placeholder — natbib \bibliography handles it
python3 - <<'PY'
import re
t=open('.SEDE.obstests.md',encoding='utf-8').read()
t=re.sub(r'\n##\s+References\s*\n+:::\s*\{#refs\}\s*\n:::\s*\n','\n',t)
open('.SEDE.obstests.md','w',encoding='utf-8').write(t)
PY

pandoc -f markdown-superscript-subscript .SEDE.obstests.md -o tex/SEDE_observational_tests.tex \
  --standalone \
  --shift-heading-level-by=-1 \
  --natbib \
  --template=tools/template_jcap.tex

# figure paths: from tex/ the repo results dir is ../../results ; prefer vector .pdf
perl -0pi -e 's#\{\.\./results/([^}]+?)\.png\}#{../../results/\1}#g' tex/SEDE_observational_tests.tex
# numeric journal (JHEP): no author-prominent form — normalise \citet -> \citep (silences natbib)
perl -0pi -e 's#\\citet\{#\\citep{#g' tex/SEDE_observational_tests.tex
# merge adjacent subscripts X_{a}_{b} -> X_{a,b} (mdmath can emit a double subscript, e.g.
# R-1_{\mathrm{cl}}_{\mathrm{stop}}); nesting-aware so \mathrm{..} content is handled.
perl -0pi -e 's/_\{((?:[^{}]|\{[^{}]*\})+)\}_\{((?:[^{}]|\{[^{}]*\})+)\}/_{$1,$2}/g' tex/SEDE_observational_tests.tex
# number display equations: pandoc emits \[ ... \]; convert to a numbered equation environment
perl -0pi -e 's/\\\[(.*?)\\\]/\\begin{equation}$1\\end{equation}/gs' tex/SEDE_observational_tests.tex

( cd tex && \
  pdflatex -interaction=nonstopmode SEDE_observational_tests.tex >SEDE_observational_tests.build.log 2>&1 ; \
  BIBINPUTS="..:$BIBINPUTS" bibtex SEDE_observational_tests      >>SEDE_observational_tests.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode SEDE_observational_tests.tex >>SEDE_observational_tests.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode SEDE_observational_tests.tex >>SEDE_observational_tests.build.log 2>&1 ) || true

if [ -f tex/SEDE_observational_tests.pdf ]; then
  echo "built paper/tex/SEDE_observational_tests.pdf"
  cp tex/SEDE_observational_tests.pdf SEDE_observational_tests.pdf       # keep the committed final PDF in sync
  grep -c "^!" tex/SEDE_observational_tests.build.log | awk '{print $1" LaTeX errors (see tex/SEDE_observational_tests.build.log)"}'
  grep -c "Warning--" tex/SEDE_observational_tests.build.log 2>/dev/null | awk '{print $1" bibtex warnings"}'
else
  echo "BUILD FAILED — see tex/SEDE_observational_tests.build.log"; grep -A2 '^!' tex/SEDE_observational_tests.build.log | head -20
fi
