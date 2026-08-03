# HANDOFF — resume point for ten-proofs explainer project

**Read this first each new session.** Not a git repo; disk state is the only state.

## Project goal
Build PhD-bridging explainer artifacts for `ten-proofs-oai.pdf` — "Ten Advances in
Mathematics and Theoretical Computer Science" (OpenAI, **249 pages**, ten independent
chapters). Audience: a reader with a PhD in *some* field but not the chapter's field.
Full spec: **`implementation-plan.html`** (repo root). Built with the
**rigorous-explainer skill** (`~/.claude/skills/rigorous-explainer/`).

Run mode: **fully autonomous** (Simon dropped the calibration gate). Deliver each
chapter to him as it completes; log judgment calls here.

## Template standard (calibrated on ch9 per Simon's feedback — FOLLOW THIS)
- 4-layer plan structure + full proofs of everything provable at graduate-core level;
  deep external theorems stated precisely as "Import N" with citations (advisory
  check_proofs flags on those boxes are intended).
- **MathJax embedded INLINE** (offline rendering — Simon saw raw LaTeX when it was a
  CDN link; his viewer blocks external scripts). Build step:
  `python figgen/inline_mathjax.py FILE.html` (bundle `.ignore/mathjax-3.2.2-tex-mml-svg.js`,
  sha384 verified). Checkers cannot parse the 2.1 MB bundle, so run ALL checks on a
  stripped copy: `python figgen/strip_mathjax.py FILE.html .ignore/check_FILE.html`.
- **High figure density** (Simon's second correction: "too many pure walls of words"):
  every abstract definition/lemma gets a diagram; no 10-line prose run without a visual.
  Range so far: 13 (ch8) to 38 (ch7).
- Literal dollar signs in prose: `<span class="tex2jax_ignore">&dollar;</span>`.
- Numeric validation section with real computed numbers + an in-page runnable listing.
- **NO problem sets** (decided; Simon can request them later as a follow-up pass).
- Dark-mode aware CSS (variables; `figure > svg` scoping). Every chapter reuses ch9's head.
- Long chapters: part-file build — write `.ignore/chNN_parts/pNN.html`, splice with
  `figgen/chNN_build.py` (replaces `<!--FIG:name-->` markers from `.ignore/chNN_figs/`).

## Current state (2026-08-02, end of session 7)
**Seven chapters DONE, verified, and delivered** — 1, 3, 4, 6, 8, 9, 10. Each is a
self-contained ~2.2–2.6 MB HTML with inline MathJax, all mechanical checks at 0 on its
stripped copy, formulas verified against rendered PDF pages, and every figure eyeballed
in a dark-mode screenshot audit. Build assets per chapter live in `figgen/`
(`chNN_figs.py` = figure data + machine checks; `chNN_build.py` = assembler) with part
files in `.ignore/chNN_parts/`.

**Ch 7 is content-complete but NOT hardened** — the live work item, see below.

Session 7 ended because the **Fable 5 usage limit was hit** mid-build; Simon then set
the default model to Opus 5. Nothing was lost — all 17 ch7 part files had landed.

## Next task
**Finish Chapter 7 (`ch07-cvp-hardness.html`), then 5 → 2, then the index page.**

Ch 7 state as of this handoff (CVP n^(1/400) hardness, PDF pp. 183–218, 36 pp):
- All 17 part files written: `.ignore/ch07_parts/p00–p16.html` (p16 = §16 scope & limits).
- All 38 figures generated in `.ignore/ch07_figs/` by `figgen/ch07_figs.py` +
  `ch07_svgs.py` + `ch07_figlib.py`.
- **Assembled and confirmed to splice cleanly**: `PYTHONIOENCODING=utf-8 python
  figgen/ch07_build.py` → 336,658 bytes, 17 parts, 37 figures, zero missing figs.
- `renders/p183–218.png` already rendered at 150 dpi for formula verification.

Remaining for ch7 — the standard hardening tail (roughly the last third of a build):
1. Run the hardening loop below on a stripped copy → `.ignore/check_ch07.html`.
2. Formula-by-formula verification against `renders/p183–218.png`.
3. Dark-mode figure audit of all 38 figures — **use segmented screenshots** (gotcha
   below); fix label collisions, clipped viewBoxes, wrong Unicode glyphs.
4. `python figgen/inline_mathjax.py ch07-cvp-hardness.html` (→ ~2.4 MB); re-verify
   rendering by screenshot after inlining.
5. Update this file, then message the coordinator (SendMessage to "main") with the path.

Then **Ch 5** (permanent circuit/formula bounds, pp. 114–153, 40 pp) and **Ch 2**
(binary & spherical codes, pp. 29–78, 50 pp — the longest; reuses the LP-certificate
spine from ch1). Finally the **index artifact**: one master page with all ten results at
a glance, significance, difficulty ratings, shared-machinery map, links (per
`implementation-plan.html` §6).

Update this file BEFORE starting each chapter, not after.

## PAUSED — 2026-08-02
Simon stopped work here: **7-day usage hit 71%** with 3 days left before the quota
resets. Nothing is broken; ch7 is exactly where "Next task" describes (parts + figures
built and spliced, hardening tail not yet run). **Resume date: 2026-08-06** (email
reminder set). Pick up at step 1 of the ch7 hardening list above.

Reminder: Google Calendar event "Resume ten-proofs explainer — finish Chapter 7
hardening", **Thu 2026-08-06, 9:00 AM PT**, with a 0-minute *email* reminder (Google
sends the mail at 9:00) plus a popup. Event id `jm4j349lgprhnnlssiujec0t24` on
simon.y.szeto@gmail.com. Delete/move it there if the resume date changes.

## Chapter status
| Ch | Topic | PDF pages | Status |
|----|-------|-----------|--------|
| 1 | Sphere packing (Cohn–Elkies rate) | 3–28 | **DONE** — `ch01-sphere-packing.html` (32 figs, 46 checks, d=6×10⁶ numerical run) |
| 2 | Binary & spherical codes | 29–78 | not started (longest; LP spine from ch1) |
| 3 | Non-sofic group | 79–95 | **DONE** — `ch03-nonsofic-group.html` (24 figs, 30 checks; §3 = property-(T)/expander primer) |
| 4 | Connes rigidity counterexample | 96–113 | **DONE** — `ch04-connes-rigidity.html` (22 figs, ~25 exact checks; reuses ch3 §3) |
| 5 | Permanent circuit/formula bounds | 114–153 | not started |
| 6 | Quantum parallel repetition | 154–182 | **DONE** — `ch06-quantum-parallel-repetition.html` (29 figs, 34 checks, 1 import only) |
| 7 | CVP n^(1/400) hardness | 183–218 | **PARTS DONE, NEEDS HARDENING** — see Next task |
| 8 | Ehrhart volume conjecture | 219–228 | **DONE** — `ch08-ehrhart-volume.html` (13 figs; §2 = lattice primer, reuse in ch7) |
| 9 | Multicolor Ramsey Rk(3) | 229–235 | **DONE** — `ch09-multicolor-ramsey.html` (15 figs; the calibration chapter) |
| 10 | Compactness & degeneracy | 236–249 | **DONE** — `ch10-compactness-degeneracy.html` (20 figs; computed doily W(2)) |

Shared primers, written once and reusable: **property (T)/expanders** = ch3 §3 (used by
ch4) · **lattices** (Blichfeldt + Minkowski proved) = ch8 §2 (reuse in ch7) ·
**LP certificates** = ch1 (reuse in ch2) · **codes** = to be written in ch2.

## Hardening loop (per chapter, all on the stripped check copy)
autolink_sections (on the real file) → strip_mathjax → checktex, checklt, check_links,
check_svg, check_code, check_prose, check_proofs → check_overlap, check_frame,
verify_dom (headless Chrome) → tall screenshot via shoot.py (`--size 900x40000`), slice
into strips, eyeball EVERY figure → fix → re-run. Delete `.bak` files at the end.

## Environment gotchas
- Read tool CANNOT render this PDF visually; page images come from PyMuPDF →
  `renders/pNNN.png` (150 dpi). 157 pages rendered so far.
- `PYTHONIOENCODING=utf-8` always, or Windows cp1252 crashes on ﬁ ligatures.
- Harness PDF metadata falsely reported 23 pages; true count is 249.
- **Full-page 60k-px Chrome screenshots hit tile-paint glitches (random black rects).**
  Audit long pages via segmented shots: inject a negative `margin-top` on `body` to
  bring the region into view, and pass `--wait 20000`. Pattern: `.ignore/ch06_seg`.
- **Bash heredocs mangle backslashes** in quoted Python blocks (`\t` → literal TAB) —
  use the Edit tool or a `.py` file for anything containing backslashes.
- verify_dom's virtual-time budget is too short for inline-MathJax files — check the
  stripped copy instead.
- SVG text: no literal `_` or `^` (check_svg hard-fails) — use tspan baseline-shift.
  Unique SVG marker ids per page (prefix c8/c9/etc.).
- CSS counters number the figures while prose cites "Fig. N" — inserting a figure shifts
  every later number, so remap prose references BEFORE inserting.
- The 6 "stray-$" advisories in every chapter are the MathJax config's own delimiter
  arrays; check_proofs advisories on Import boxes are intended. Neither is a bug.

## Session-transient scratch (regenerate; the durable record is the chapter HTML)
- `figgen/chNN_figs.py` → figure data + machine checks; `figgen/ch07_svgs.py` → SVG
  emitters; `figgen/chNN_build.py` → splices `.ignore/chNN_parts/*.html` with
  `.ignore/chNN_figs/*.svg` into the chapter HTML. Idempotent; rerun any time.
- `figgen/inline_mathjax.py` + `strip_mathjax.py` → the inline/strip pair every chapter
  depends on. These two matter most; the rest are per-chapter.
- `.ignore/mathjax-3.2.2-tex-mml-svg.js` — **KEEP**. Needed by inline_mathjax.py; if
  lost, re-fetch from cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-svg.js
  (sha384-msWocAZtTDh+149KSxjbRTGVfGbDjiff/j0JX1iu/A/GxxTWSu6ozNemlfkag/8d).
- `.ignore/ch0*_*.png`, `.ignore/check_ch0*.html` — audit artifacts, regenerable.
- Pyright warnings on `figgen/*.py` (unused vars, dynamic dict/tuple typing) are IDE
  noise — the scripts are actually run and only real outputs are quoted.

## How to work
- Deliverables are local self-contained HTML; deliver each to Simon as it completes.
- "Done" = formulas verified against the *rendered* PDF page (not extracted text, which
  flattens sub/superscripts) AND a clean dark-mode figure audit.
- Budget: one big chapter per session is realistic (ch1 and ch6 each consumed a full
  context). If context runs low, leave a clean resume point here rather than a
  half-built chapter — that protocol has now survived four interruptions without loss.
