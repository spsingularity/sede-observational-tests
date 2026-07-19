#!/usr/bin/env python3
"""
Build-time preprocessor: wrap the paper's *informal inline math notation* in $...$ so
xelatex renders it (the source SEDE.md keeps the readable plain notation, which also
displays fine on GitHub). We ONLY touch runs that contain an ASCII '_' or '^' — those are
the cases pandoc mis-renders (literal underscore / stray caret). Lone Greek (ρ, Δ) and
Unicode super/sub (², ₀) already render in STIX Two Text and are left untouched.

Protected and never altered: fenced code, inline `code`, existing $...$/$$...$$ math,
and markdown image/link syntax. Run:  python mdmath.py IN.md OUT.md
"""
import re, sys

GREEK = {
    'α':r'\alpha ','β':r'\beta ','γ':r'\gamma ','δ':r'\delta ','ε':r'\epsilon ',
    'ζ':r'\zeta ','η':r'\eta ','θ':r'\theta ','ι':r'\iota ','κ':r'\kappa ',
    'λ':r'\lambda ','μ':r'\mu ','ν':r'\nu ','ξ':r'\xi ','π':r'\pi ','ρ':r'\rho ',
    'σ':r'\sigma ','τ':r'\tau ','υ':r'\upsilon ','φ':r'\phi ','ϕ':r'\phi ',
    'χ':r'\chi ','ψ':r'\psi ','ω':r'\omega ','ϵ':r'\epsilon ','ϑ':r'\vartheta ',
    'Γ':r'\Gamma ','Δ':r'\Delta ','Θ':r'\Theta ','Λ':r'\Lambda ','Ξ':r'\Xi ',
    'Π':r'\Pi ','Σ':r'\Sigma ','Φ':r'\Phi ','Ψ':r'\Psi ','Ω':r'\Omega ',
}
SUP = {'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4','⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁻':'-','⁺':'+'}
SUB = {'₀':'0','₁':'1','₂':'2','₃':'3','₄':'4','₅':'5','₆':'6','₇':'7','₈':'8','₉':'9'}
OPS = {'∝':r'\propto ','≈':r'\approx ','≲':r'\lesssim ','≳':r'\gtrsim ','≡':r'\equiv ',
       '≤':r'\leq ','≥':r'\geq ','≪':r'\ll ','≫':r'\gg ','≠':r'\neq ','∈':r'\in ',
       '×':r'\times ','·':r'\cdot ','−':'-','∞':r'\infty ','√':r'\surd ','ℏ':r'\hbar ','ℓ':r'\ell ',
       '⊙':r'\odot ','⊕':r'\oplus ','⊗':r'\otimes ','∂':r'\partial ','∇':r'\nabla '}

# characters allowed inside an inline math run (besides spaced operators and brace groups)
MATHCHARS = r"A-Za-z0-9_\^{}()\[\]/+\.'" + "".join(re.escape(c) for c in
            list(GREEK) + list(SUP) + list(SUB) + list(OPS))
# a run: starts with a math char, continues with math chars, {..} groups (spaces ok inside),
# or spaced binary operators
RUN = re.compile(
    r"[A-Za-zΑ-Ωα-ω0-9(\[]"                        # start
    # continue with math chars, {..} groups, or spaced binary ops. NB: ASCII '-' is a prose
    # hyphen (the paper writes math minus as '−' U+2212), so it is NOT a connector.
    # NB: operator padding is [ \t]* (NOT \s*) so a run never spans a newline — inline $...$
    # math cannot cross a line break.
    r"(?:\{[^}]*\}|[" + MATHCHARS + r"]|[ \t]*[−=+/<>][ \t]*|[ \t]*[∝≈≲≳≡≤≥≪≫≠∈][ \t]*)*"
)

def to_latex(s):
    # multi-letter ASCII uppercase = acronym/label, set upright roman so it is not
    # italicised as a product of single-letter variables (SEDE, CDM, DIC, TTTEEE, II).
    # Runs first, on the raw text, before any backslash commands are introduced; single
    # uppercase letters (H, T, S, …) are left as italic math variables.
    s = re.sub(r"[A-Z]{2,}", lambda m: r"\mathrm{" + m.group() + "}", s)
    # unicode super/sub -> ^{..}/_{..}
    s = re.sub("[" + "".join(SUP) + "]+", lambda m: "^{" + "".join(SUP[c] for c in m.group()) + "}", s)
    s = re.sub("[" + "".join(SUB) + "]+", lambda m: "_{" + "".join(SUB[c] for c in m.group()) + "}", s)
    # group sub/superscripts on the RAW unicode (before Greek map) so multi-char scripts
    # like ^φφ or _grav stay in one group. Latin subscript words -> roman.
    s = re.sub(r"_([A-Za-z][A-Za-z0-9]*)", lambda m: r"_{\mathrm{" + m.group(1) + "}}", s)
    s = re.sub(r"_([0-9])(?![{}])", r"_{\1}", s)
    s = re.sub(r"\^(\{[^}]*\}|[A-Za-zΑ-Ωα-ω0-9]+)",
               lambda m: "^" + m.group(1) if m.group(1).startswith("{") else "^{" + m.group(1) + "}", s)
    # now translate the remaining unicode letters/operators to LaTeX commands
    for u, t in {**GREEK, **OPS}.items():
        s = s.replace(u, t)
    return s

# any of these Unicode glyphs in a run is an unambiguous math signal (Greek letter,
# unicode super/sub, or a mapped operator) -> the run is real math, wrap it.
UNIMATH = frozenset(set(GREEK) | set(SUP) | set(SUB) | set(OPS))

def convert_run(m):
    s = m.group(0)
    # wrap a run only if it carries a math signal: an ASCII '_'/'^' (which pandoc mis-renders)
    # or a Unicode math glyph. A pure-ASCII prose run (no signal) renders fine -> leave it.
    if "_" not in s and "^" not in s and not any(c in UNIMATH for c in s):
        return s
    # peel leading prose words a spaced operator dragged in: "crosses −1" -> lead="crosses ",
    # math="−1"; "from ≈0" -> lead="from ". Only multi-letter lowercase words (real English)
    # are peeled -- single-letter leads like "w = −1" are variables and stay in the math.
    lead = ""
    while True:
        mm = re.match(r"[a-z]{2,}\s+", s)
        if not mm: break
        lead += mm.group(0); s = s[mm.end():]
    if not s or ("_" not in s and "^" not in s and not any(c in UNIMATH for c in s)):
        return lead + s          # nothing math-worthy left after peeling
    # don't drag trailing prose punctuation / unbalanced brackets into the math
    trail = ""
    while s:
        c = s[-1]
        if c in " \t": s = s[:-1]                               # trailing whitespace -> drop
        elif c in ".,;:": trail = c + trail; s = s[:-1]         # punctuation -> keep as prose
        elif c in "−=+/<>∝≈≲≳≡≤≥≪≫≠∈": trail = c + trail; s = s[:-1]  # dangling op -> out of math
        elif c in "_^": s = s[:-1]                              # dangling sub/sup marker -> drop (safety)
        elif c == ")" and s.count("(") < s.count(")"): trail = c + trail; s = s[:-1]
        elif c == "]" and s.count("[") < s.count("]"): trail = c + trail; s = s[:-1]
        elif c == "}" and s.count("{") < s.count("}"): trail = c + trail; s = s[:-1]
        else: break
    s = s.rstrip()
    if not s: return m.group(0)
    # .strip(): a Greek substitution leaves a trailing space (Λ -> "\Lambda "); a space just
    # inside the closing $ makes pandoc reject the whole thing as math. Must not happen.
    return lead + "$" + to_latex(s).strip() + "$" + trail

def protect(text):
    spans, out, i = [], [], 0
    pat = re.compile(r"```.*?```|``.*?``|`[^`]*`|\$\$.*?\$\$|\$[^$\n]*\$|!\?\[[^\]]*\]\([^)]*\)|\[[^\]]*\]\([^)]*\)",
                     re.S)
    for mt in pat.finditer(text):
        out.append(text[i:mt.start()])
        out.append(f"\x00{len(spans)}\x00"); spans.append(mt.group(0)); i = mt.end()
    out.append(text[i:])
    return "".join(out), spans

def restore(text, spans):
    return re.sub(r"\x00(\d+)\x00", lambda m: spans[int(m.group(1))], text)

def main(inp, outp):
    text = open(inp, encoding="utf-8").read()
    body, spans = protect(text)
    body = RUN.sub(convert_run, body)
    open(outp, "w", encoding="utf-8").write(restore(body, spans))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
