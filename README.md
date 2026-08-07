# Ten Advances in Mathematics and Theoretical Computer Science — explainers

Self-contained HTML explainers for the ten results surveyed in OpenAI's report
*"Ten Advances in Mathematics and Theoretical Computer Science"* (249 pages, ten
independent chapters).

**Source paper:** <https://openai.com/index/ten-advances-in-mathematics-and-theoretical-computer-science/>
The PDF itself is not redistributed here — download it from the link above. Every
chapter of this repo cites the page range it covers.

## What this is

Each chapter of the report states a recent theorem and sketches its proof at a level
aimed at people already inside the relevant subfield. These explainers rewrite that
material for a different reader: **someone with a PhD in *some* field, but not
*this* field.** A probabilist reading the operator-algebras chapter, a systems
person reading the lattice chapter.

Concretely, each explainer is one HTML file you can open offline, with no network
access, and read end to end.

## What's actually new here (vs. reading the paper)

- **Four-layer structure per chapter** — plain-language statement → the objects
  defined from scratch → the proof → what the result does *not* say. The paper
  assumes the first two layers; these files build them.
- **Every provable step is proved.** Anything provable from graduate-core
  background is written out in full. Genuinely deep external theorems are not
  hand-waved either — they are stated precisely in a boxed "Import *N*" with a
  citation, so the reader always knows exactly which black boxes the argument rests
  on. Chapter 6, for instance, needs only one such import.
- **Figure-dense, not wall-of-text.** 13 to 38 hand-built SVG diagrams per chapter
  (256 across the ten finished chapters). Every
  abstract definition or lemma gets a picture; no long prose run goes without a
  visual. The figures are generated from data by the scripts in `figgen/`, so the
  geometry in a diagram is computed, not sketched by eye.
- **Numeric validation sections.** Each chapter recomputes its key quantities with
  real numbers and ships the runnable listing in-page — e.g. the sphere-packing
  chapter runs the Cohn–Elkies rate out to dimension 6×10⁶, the CVP chapter
  executes its entire 3SAT→lattice construction on a toy formula in F₁₆ (3120
  unknowns, 8949 equations, rank computed exactly) and re-proves the paper's
  parameter inequalities at N = 100 in exact integer arithmetic, and the
  compactness chapter computes the W(2) generalized-quadrangle ("doily") directly,
  and the permanent chapter re-runs its entire coefficient-independence construction
  at n = 13 — all 2,401 entries of a 49×49 Jacobian rebuilt by brute-force permanents,
  factored, and inverted in exact arithmetic — plus its block-cancellation identity in
  exact ℤ[ζ] arithmetic. The codes chapter verifies its imported Kravchuk/Weyl
  representation-theory formulas in exact rational arithmetic (transition weights
  summing to exactly 1 across random Young diagrams of both parities) and certifies
  the paper's kissing-number exponent 0.39661 from explicit five-parameter tuples.
- **Shared primers written once and reused.** Property (T) / expanders (ch. 3 §3,
  reused by ch. 4) · lattices with Blichfeldt and Minkowski proved (ch. 8 §2,
  reused by ch. 7) · LP certificates (ch. 1, reused by ch. 2).
- **Mechanically checked.** Each file passes a battery of automated checks — LaTeX
  well-formedness, SVG validity, link integrity, proof-box completeness, element
  overlap and frame containment in a headless browser — and every formula is
  verified against the *rendered* PDF page rather than extracted text, which
  flattens sub- and superscripts.
- **Offline by construction.** MathJax is inlined into each file (~2.3 MB each), so
  the math renders with no CDN and no scripts fetched from anywhere.
- Dark-mode aware throughout.

Built with the `rigorous-explainer` workflow.

## Read them

The chapters are live on GitHub Pages — click a title below. (GitHub itself won't
render a raw `.html` file, so read them at these links rather than through the file
list.) They are also self-contained: download one and it works offline.

## Status — chapters complete

All ten chapters are finished and verified. The cross-chapter index page is
still to come.

| Ch | Topic | PDF pp. | Figures | Status |
|----|-------|---------|---------|--------|
| 1 | [Sphere packing (Cohn–Elkies rate)](https://az9713.github.io/ten-math-cs-advances/ch01-sphere-packing.html) | 3–28 | 32 | ✅ done |
| 2 | [Binary & spherical codes](https://az9713.github.io/ten-math-cs-advances/ch02-binary-spherical-codes.html) | 29–78 | 34 | ✅ done |
| 3 | [A non-sofic group](https://az9713.github.io/ten-math-cs-advances/ch03-nonsofic-group.html) | 79–95 | 24 | ✅ done |
| 4 | [Connes rigidity counterexample](https://az9713.github.io/ten-math-cs-advances/ch04-connes-rigidity.html) | 96–113 | 22 | ✅ done |
| 5 | [Permanent circuit/formula bounds](https://az9713.github.io/ten-math-cs-advances/ch05-permanent-bounds.html) | 114–153 | 29 | ✅ done |
| 6 | [Quantum parallel repetition](https://az9713.github.io/ten-math-cs-advances/ch06-quantum-parallel-repetition.html) | 154–182 | 29 | ✅ done |
| 7 | [CVP hardness to factor n^(1/400)](https://az9713.github.io/ten-math-cs-advances/ch07-cvp-hardness.html) | 183–218 | 38 | ✅ done |
| 8 | [Ehrhart volume conjecture](https://az9713.github.io/ten-math-cs-advances/ch08-ehrhart-volume.html) | 219–228 | 13 | ✅ done |
| 9 | [Multicolor Ramsey R_k(3)](https://az9713.github.io/ten-math-cs-advances/ch09-multicolor-ramsey.html) | 229–235 | 15 | ✅ done |
| 10 | [Compactness & degeneracy](https://az9713.github.io/ten-math-cs-advances/ch10-compactness-degeneracy.html) | 236–249 | 20 | ✅ done |

The build spec these were written against is also live:
[implementation-plan.html](https://az9713.github.io/ten-math-cs-advances/implementation-plan.html).

**Still missing:** the planned index artifact — one page showing all ten results at a
glance with significance, difficulty ratings, and a map of which machinery is shared
between chapters.

No problem sets in any chapter; that was a deliberate scope decision.

## Layout

```
ch0*.html               the deliverables — open any one directly in a browser
implementation-plan.html  the full spec these were built against
figgen/
  chNN_figs.py          figure data + the machine checks for that chapter
  chNN_build.py         assembler: splices part files with generated SVGs
  chNN_figlib.py        per-chapter SVG helpers
  inline_mathjax.py     bundles MathJax into a finished chapter (sha384-verified)
  strip_mathjax.py      inverse; produces a checkable copy, since the 2.1 MB
                        bundle defeats the static checkers
.ignore/chNN_parts/     the prose, one part file per section (build input)
.ignore/chNN_figs/      generated SVGs (build input, regenerable)
```

Rebuilding a chapter is `python figgen/chNN_build.py` followed by
`python figgen/inline_mathjax.py chNN-*.html`. Both are idempotent. On Windows set
`PYTHONIOENCODING=utf-8` first, or the ﬁ ligatures in the source text crash cp1252.

Not tracked: the source PDF, its extracted text, page renders, and audit
screenshots — all regenerable from the paper.
