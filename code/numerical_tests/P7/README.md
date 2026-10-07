# P7: numerical test of Section 9 (the per-scale count)

This is a test on small analogues of the construction. It is **not part of the proof**. The flow tracer and the
checks were written from scratch for this test on 2026-10-06, against the research draft of Sections 3 and 9, and
they use no other code. The statement numbers below are those of the paper. All geometry uses exact rationals;
only comparisons with irrational numbers use mpmath at 60 digits.

Results recorded on 2026-10-06 are in `results/original/`. Only summaries of the large campaign records are
included (the six campaign files have 44 MB).

## What is tested

For each configuration the master flow F⁰ is traced exactly, and for several scales s the scale flow F_s
(Definition 3.21) is derived from it. Then:

| claim | statement of the paper | what is checked |
|---|---|---|
| (a) | Theorem 9.1(a) | \|Z(s)\| ≤ (εs + b_W)(N(s) + Λ\*(s)); `a_ratio` = left side / right side must be ≤ 1. The steps of its proof are checked at sampled heights: Step 1 (Z ⊂ ℝ × V_s), Step 3 (Σ_Z c_Z ≤ Sh), Step 4 (Sh ≤ L_T + Λ, Theorem 6.9(f); L_T ≤ N(s); the loss bounds of Proposition 6.10), and the identities of Theorem 6.9(b),(c). |
| (b) | Step 2 of the proof of Theorem 9.1 | No path of F_s makes an alive contact with relint bot(Z), Z ∈ Z(s). At every alive contact, the width test of Theorem 4.11 with w(s, Z) of Lemma 9.3(e) passes (ratio ≤ 1). |
| (c) | Theorem 9.2(a) | I ≤ J ≤ \|Z(s)\| β̂((1−ε)s), with I = ∫(1 − ω₀ − E)₊ and J = Σ_Z ∫ c_Z, and the resulting lower bound for \|Z(s)\|. Its Step 1 uses [R, Lemma 3.5] (`lem35`), which is also checked. |
| (d) | Lemma 9.6(c), Lemma 9.7 | {T\* ≤ s} (Definition 3.22) equals the set of x whose path is T_s-terminated in F_s, and equals the set characterised in Lemma 9.6(b),(c); the symmetric differences must have measure exactly 0. T_s is nested in s. L_T ≤ N(s). |

Self-checks of the tracer: an independent point-by-point tracer is compared with the flow tracer at sampled x.
Every M-terminated piece must have an R1 winner with lexicographically smaller (g, x). Conservation holds on F⁰.

Negative control: counting all squares inside V_s instead of Z(s) (`negctl_ratio_Zall`) must be able to exceed
the bound of (a), and it does.

### The scaled analogue

The parameters of the paper (y₀ ≈ 4·10¹⁰, angles ≈ 10⁻⁶, δ = 10⁻⁵) cannot be simulated. The test uses
α(s) = c₀ s^{−1/2} and β̂(t) = c₁ t^{−3/4} with the following parameter ranges in the records: k = 10–95,
δ = 0.01–0.3, c₀ = 0.01–0.6, c₁ = 0.004–0.06, y₀ = 1–3, ε = 0.2–0.76. The numerical constants of the paper
that rely on the small angles are replaced by their exact counterparts (see `p7x_meas.py`):

- b := β̂((1−ε)y₀) bounds a(Z) for Z ∈ Z(s) (it replaces 10⁻⁴), and h := 1 + b replaces 1.0001;
- w(s, Z) = (sec α(s) − 1) q_y + δ + sin a + 1 − cos a is used in place of 0.50001 α(s)² q_y + … of Lemma 9.3(e);
- w₀ is recomputed as sup_s (sec α(s) − 1)(s + h) + δ + sin b + 1 − cos b + 10⁻¹², rounded up;
- V_s = [(1−ε)s − h, s + h], so b_W := 2h replaces 2.0002.

## Files

| file | content |
|---|---|
| `p7x_core.py` | exact flow tracer (ray steps, priority D > W > contact > H, R1, R2) and an independent point tracer |
| `p7x_meas.py` | F_s, T\*, Z(s), N(s), Λ₀, line quantities, width test |
| `p7x_eval.py` | evaluation of one configuration: all checks (a)–(d) (the field names are listed at its top) |
| `p7x_configs.py` | configuration generators |
| `p7x_run.py` | campaign driver (families `blocked`, `blocked2`, `drift`, `random`, `round2`; see its header) |
| `p7x_targeted.py`, `p7x_targeted2.py` | hand-built degenerate cases: death exactly at a contact (T1), exits and contacts exactly at height s + 2 (T2), a ramp exactly at the boundary of H (T3), floor/wall contacts (T4), tightness of the width test (T5), the ceiling version (T6). In the run of 2026-10-06 the T1 and T4 configurations of `p7x_targeted.py` were rejected as invalid (overlapping squares). `p7x_targeted2.py` repeats them with corrected configurations. |
| `p7x_t7.py` | tightness of the width test from the gap side (108 configurations) |
| `p7x_maxratio.py` | deterministic attempt to push the ratio of (a) towards 1 (blocked columns with many rows) |
| `p7x_merge.py` | merge-heavy configurations: every M piece and the R1 winners are checked by the point tracer |
| `p7x_example.py` | re-traces one recorded case where F⁰ reaches a square of Z(s) but F_s does not, and prints exact data |

On a clean Linux machine (2026-10-07, `../rerun_linux_2026-10-07/`) every row was reproduced; the server ran 45 configurations of campaign D in the 90-second limit, identical to the first 45 recorded ones.
| `p7x_summary.py` | aggregates the JSON-lines outputs |
| `results/original/summary_all.txt` | output of `p7x_summary.py` over campaigns A–E, R2 and the targeted, maxratio and t7 records |
| `results/original/camp*.log`, `merge71.log`, `merge72.log`, `targeted.log`, `maxratio.log` | console outputs |
| `results/original/records_targeted_merge.zip` | the records `targeted.jsonl`, `targeted2.jsonl`, `maxratio.jsonl`, `t7.jsonl`, `merge_71.jsonl`, `merge_72.jsonl` |
| `results/original/example_case.jsonl` | the one record of campaign B used by `p7x_example.py` (copied unchanged) |
| `results/original/example_out.txt` | output of `p7x_example.py` |
| `results/original/resource_log.txt` | launch times, memory, and stop times of the runs of 2026-10-06 |

The scripts write into a subfolder `out/`, which they create. The six campaign records (`campA.jsonl` …
`campE.jsonl`, `campR2.jsonl`; 3,014 records, 44 MB) are not included. `summary_all.txt` summarises them,
and `p7x_run.py` regenerates them (see "Reproducibility").

## Quick run

Python 3.12 with `mpmath` (numpy is not needed). From this folder:

```
python p7x_run.py blocked2 505 90 quick_D.jsonl
python p7x_summary.py quick_D.jsonl
python p7x_t7.py
python p7x_targeted2.py
python p7x_targeted.py
python p7x_example.py
```

Measured on 2026-10-07 on a Windows 11 laptop (Python 3.12.10, one process at a time), in two sessions. The
times depend on the load of the machine.

| command | time | peak memory | result |
|---|---|---|---|
| `p7x_run.py blocked2 505 90 quick_D.jsonl` | 91 s (time limit 90 s) | 38 MB | 73 and 107 configurations in the two sessions (the number depends on the speed of the machine); identical to the first records of campaign D |
| `p7x_summary.py quick_D.jsonl` | < 1 s | 18 MB | 292 and 428 cases; every counter `viol_*` is 0; largest ratio of (a) 0.9636 |
| `p7x_t7.py` | 2–3 s | 22 MB | `T7 configs 108 violations 0`, best width ratio 0.99953; identical to the record |
| `p7x_targeted2.py` | < 1 s | 22 MB | identical to the record |
| `p7x_targeted.py` | 11–14 s | 27 MB | identical to the record (`T6 ceiling/reflection checks 16 mismatches []`) |
| `p7x_example.py` | < 1 s | 22 MB | identical to `results/original/example_out.txt`, except the label of the first line |

Longer re-runs made for this release: `p7x_maxratio.py 100` (86 s, 47 MB) gives all 56 cases identical to the
record, with best ratio 0.98607. `p7x_merge.py 71 90` (91 s, 36 MB) gives its first 84 configurations
identical to the record. The memory guard of `p7x_run.py` works only on Windows and is inactive elsewhere.

## Results recorded on 2026-10-06

From `summary_all.txt` (campaigns A–E, R2, targeted, maxratio, t7): 3,203 valid configurations (556,622
squares) and 15,450 cases (configuration, s). Z(s) is nonempty in 4,991 cases, with 538,030 elements in all.
**No violation** of (a)–(d) or of any pointwise check.

- (a): the largest ratio is 0.98607, from `p7x_maxratio.py` (blocked columns, m = 10 rows). The ratio comes
  close to 1 for blocked columns with many rows, so the inequality is nearly sharp in these analogues.
- (b): 200,967 width tests at alive contacts, largest ratio 0.99953 (`p7x_t7.py`). No alive contact with a square
  of Z(s). In 614 cases F⁰ reaches a square of Z(s) while F_s does not, so the truncation of F_s does real
  work (`example_out.txt` shows one such case with exact coordinates).
- (c): I > 0 in 163 cases; the largest ratio I/J is 0.98969 and the largest ratio of the lower bound to |Z(s)|
  is 0.0784.
- (d): the symmetric differences are exactly 0 in all cases. T_s is nested in s in all 1,585 configurations where
  this was recorded.
- Tracer: 62,098 sampled x agree with the point tracer, with no anomaly. In the merge runs (seeds 71, 72; 332
  configurations; `merge71.log`, `merge72.log`), all 1,107 M-terminated pieces and 14,995 winner entries pass
  the R1 checks, and the evaluation finds no violation.
- Negative control: counting all squares in V_s exceeds the bound of (a) in 2,362 cases (up to 7.94), so the
  check can fail.
- (C3) is not used in Section 9 (Remark 9.14). In 652 cases with Z(s) ≠ ∅, (C3) is violated, and 474 cases have a
  square of Z(s) with tilt ≥ α(s). All checks pass there too. 2,467 cases have ω₀ < 1.

## Reproducibility

On 2026-10-07 the records were compared field by field with re-runs of the scripts in this folder. Run time and
memory fields were ignored, and the renamed key (last section) was mapped.

| record | re-run | result |
|---|---|---|
| campaign D (`blocked2 505`), E (`drift 606`), R2 (`round2 404`) | first 107, 98, 49 configurations | identical |
| campaigns A (`blocked 101`), B (`drift 202`), C (`random 303`) | first 91, 116, 134 configurations | identical in every recorded field. The current `p7x_eval.py` adds three fields (`zt_over_alpha`, `extra_constraint_holds`, `d_Ts_nested`). |
| `t7`, `targeted`, `targeted2`, `maxratio` | full runs | identical |
| `merge_71` | first 84 configurations | identical |
| `example_out.txt` | `p7x_example.py` on `example_case.jsonl` | identical except the label of the first line |
| `summary_all.txt` | `p7x_summary.py` over the original records | identical |

Campaigns A, B and C were started a few minutes before the last change of `p7x_eval.py` (see
`resource_log.txt`). They therefore ran with an earlier version that did not record the (C3) flag, the tilt
ratio and the nesting of T_s. The configurations and every recorded quantity are reproduced exactly. But the
nesting check of Lemma 9.6(c) and the (C3) statistics above come only from the other runs.

The campaigns stop at a time limit, and a configuration that hits the per-configuration time limit of
`p7x_eval.py` would change the random stream after it. This did not happen in the re-runs above.

## Names in the recorded files

- `d_lemB_symdiff` (in the records and `example_case.jsonl`) is the quantity now called `d_Tchar_symdiff`.
  "Lemma B" was the research draft's name for what is now Lemma 9.6. All recorded values are 0.
- "DEFS P2 2nd inequality" in the first line of `example_out.txt` is the inequality of Remark 9.4. "DEFS" was the
  definitions file of the research draft.
- `lem35` is [R, Lemma 3.5]. `P4star` is the shadow inequality of Theorem 6.9 (P4\*). `extra_constraint_holds`
  says whether (C3) holds at s.

## Limitations

- These are small analogues with large angles and gaps. They test the composition of the definitions and the
  conclusions of Theorems 4.11 (P2) and 6.9 (P4\*), not their proofs, and not the constants of the real range.
- Sets of measure zero in x are seen only through sampled point traces.
