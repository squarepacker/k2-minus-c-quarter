# Numerical tests for Section 6 (P4\*: the shadow inequality)

Exact-arithmetic tests of Theorem 6.9 (P4\*), Lemma 6.7 (the exceptional set) and Proposition 6.10 (the wall term), on small configurations: the scaled analogues have container side k = 4 to 11 (20 in the wall sweep), δ from 0.05 to 4k (44; 80 in the wall sweep) and α_F up to almost π/4; the families marked `_paper` use the constants of the paper, δ = 10⁻⁵ and α_F = 2 atan(7.5·10⁻⁷) ≈ 1.5·10⁻⁶, with k = 6 to 8. These are tests, not part of the proof; "no counterexample" is the outcome of a search and is provisional.

Squares have rational centres and rational t = tan(φ/2), so every computation is exact (`fractions`). Two independent tracers are used: a beam tracer (x-intervals with affine data, R1 resolved by a sweep) and a per-point tracer (one floor point at a time, R1 by an explicit backward search of all live arrivals); they are compared at sampled points.

## Files

| File | Content |
|---|---|
| `core.py` | the two tracers; the quantities of Theorem 6.9 at a height y (identities, Ov, B, U, margin, losses); the exceptional set; truncation to the scale flows F_s |
| `generators.py` | configuration families (valleys, zigzags, deaths, near 45°, walls, floor, jams, towers, rows, drift to a wall, steep DAG configurations) |
| `run_tests.py` | driver: all checks on each configuration; `unit` runs hand-computed checks (including the configuration of Proposition 6.17) |
| `negative_control.py` | rule R1 switched off: the final inequality must then fail (Theorem 6.9(d),(f)) |
| `ties.py` | hand-made exact ties of the priority rules (degenerate configurations) |
| `wall_sweep.py` | sharpness of the wall bound of Proposition 6.10 on an explicit configuration |
| `exceptional_types.py` | which type of exceptional height (Lemma 6.7(iv)) actually breaks the inequality |
| `summarize.py` | aggregates JSONL results per family |
| `results/original/` | the recorded outputs (see below), unchanged |

Requirements: Python 3.12, standard library only. Outputs go to `out/` (for `run_tests.py`, to the path given on the command line).

## What is checked

At sampled heights y (random heights, exceptional heights and vertex heights together with their neighbours at distance 10⁻¹² or 10⁻¹³, heights of event points) for F⁰ and for several scale flows F_s per configuration (about 9 on average; the thresholds and heights of the F_s include exceptional heights of F⁰ and tilts of passable squares):
- the identities of Theorem 6.9(b), (c), (e) and Sh = L_T + Λ − U; Ov = 0 (Theorem 6.9(d)); B = 0 and the inequality of Theorem 6.9(f) at every non-exceptional height (at exceptional heights the inequality may fail; failures there are counted, not flagged);
- L_T^{F_s}(y) ≤ N(s) (Theorem 6.9(f)) and the bounds of Proposition 6.10 (L_W^{F_s}(y) ≤ 2y tan α(s), the D, M, E losses, Ov_gap^{F_s} ≤ Ov_gap^{F⁰});
- Exc(F_s) ⊆ Exc(F⁰) ∩ (0, h_s) (Lemma 6.7(iii)) and the types of exceptional heights (Lemma 6.7(iv));
- slopes of the crossing maps ≥ 1 (Lemma 6.1(iii)); the partition of (0, k) into pieces; agreement of the beam tracer with the per-point tracer for F⁰;
- with R1 switched off (a negative control inside every configuration with merges), the heights with Ov > 0 are counted; `negative_control.py` uses configurations built so that the final inequality must then fail.

Not tested: squares at exactly 45° (impossible with rational t), the actual scale k of order 10¹², and an independent per-point check of the scale flows F_s (the per-point tracer checks F⁰ only).

## Quick checks

Measured on the author's Windows PC on 2026-10-07 (other programs were running at the same time):

| Command | Time, peak memory | Result |
|---|---|---|
| `python run_tests.py unit` | 0.8 s, 20 MB | 11 hand-computed checks, all PASS |
| `python summarize.py results/original/res_scaled.jsonl results/original/res_merge.jsonl results/original/res_defs.jsonl results/original/res_dag.jsonl results/original/res_manyvar.jsonl results/original/res_ties.jsonl` | 0.2 s, 12 MB | the printed summary is identical to `results/original/summary_all.txt` |
| `python summarize.py results/original/res_dag.jsonl results/original/res_manyvar.jsonl results/original/res_ties.jsonl` | 0.2 s, 12 MB | identical to `results/original/summary_round2.txt` |
| `python negative_control.py` | 6 s, 18 MB | identical to `results/original/res_negctl.json` (40 records) |
| `python ties.py` | 4 s, 19 MB | identical to `results/original/res_ties.jsonl` (7 records) |
| `python wall_sweep.py` | 1 s, 18 MB | identical to `results/original/res_wallsweep.json` (48 records) |
| `python exceptional_types.py 60` | 60 s, 24 MB | 666 configurations; the same pattern as the recorded 420-second run (below) |
| `python run_tests.py dag,near45,dag 6 4004 120 out/quick_dag.jsonl 60 manyvar` | 15 s, 25 MB | the first 6 records of `res_dag.jsonl`, identical |
| `python run_tests.py valley,deaths,floor,rows,zigzag,floor_paper 6 5005 120 out/quick_manyvar.jsonl 60 manyvar` | 27 s, 26 MB | the first 6 records of `res_manyvar.jsonl`, identical |
| `python run_tests.py valley,zigzag,rows,near45,deaths,floor,jam,tower,wall,wallmax 10 1001 100 out/quick_scaled.jsonl 30` | 25 s, 23 MB | 10 configurations, 0 violations |
| `python run_tests.py valley_paper,zigzag_paper,rows_paper,deaths_paper,floor_paper,jam_paper,tower_paper 7 2002 100 out/quick_paper.jsonl 30` | 109 s, 37 MB | 3 configurations within the budget (the paper's constants are slow), 0 violations |

"Identical" means equal record by record, ignoring the fields `time` and `rss_mb`, with `_defs` read as `_paper` (see below). `run_tests.py` appends to its output file. A `summarize.py` over the quick-run outputs gave 40 configurations, 0 violations, 0 errors and 0 point mismatches. On a clean Linux machine (2026-10-07, `../rerun_linux_2026-10-07/`) every row was reproduced, and `exceptional_types.py 420`, the recorded budget, ran 4,731 configurations with the same pattern as the recorded run (no failure at the heights of constant-death type; failures, with least margin −4, only at the heights of sides of axis-parallel squares).

## Recorded results (`results/original/`)

Totals of `summary_all.txt` (all six result files): 1,461 configurations (1,267 scaled, 194 with the paper's constants), 14,230 squares, 165,485 heights checked for F⁰, 12,652 scale flows F_s with 633,460 heights, 99,489 point checks with 0 mismatches, 316 configurations with merges; **0 violations, 0 errors**. At exceptional heights the inequality failed 773 times for F⁰ and 5,794 times for F_s, which Theorem 6.9(f) allows. Inside `run_tests.py`, the negative control (R1 off) produced Ov > 0 at 2,802 of 5,691 heights; in these random configurations the final inequality itself did not fail.

- `negative_control.py`: with R1 off, the final inequality fails at non-exceptional heights in 20 of the 40 configurations (780 heights, Ov up to 0.042); with R1 on, in none (Ov = 0, all identities hold).
- `wall_sweep.py`: the largest ratio L_W/(2y tan α_s) is 0.948623 (k = 20, a = 0.01), equal to the hand-computed value within 2.2·10⁻¹⁶.
- `exceptional_types.py` (420 s, 2,938 configurations): at the 5,096 exceptional heights of constant-death type the inequality never failed; at the 8,164 heights of bottom or top sides of axis-parallel squares it failed 2,818 times (smallest margin −4).
- `ties.py`: 7 configurations, 0 violations.
- Outside the hypotheses of the paper: in the steep `dag` family (α_F close to π/4, δ ≥ 0.5) the gap-crosser density reaches 2.106, so three paths can cross the same point (M_g ≥ 3). This is no contradiction: Theorem 6.9 does not use M_g ≤ 2, and Lemma 7.7 assumes δ ≤ ½ cos 2α_max. With the paper's constants the density stayed ≤ 2.0.

| File | Produced by (command from `memlog.txt`, with the present script names; to rerun, write `_paper` for `_defs`) | Run (2026-10-06) |
|---|---|---|
| `smoke.jsonl` | `run_tests.py valley,zigzag,deaths,near45,wall,floor,jam,tower 8 11 400 smoke.jsonl 30` | 18:16 |
| `smoke2.jsonl` | `run_tests.py rows,rows_defs,wallmax,valley_defs,deaths_defs,jam_defs 6 21 400 smoke2.jsonl 40` | 18:18 |
| `res_scaled.jsonl` | `run_tests.py valley,zigzag,rows,near45,deaths,floor,jam,tower,wall,wallmax 100000 1001 2100 res_scaled.jsonl 60` | 18:19–18:54 |
| `res_defs.jsonl` | `run_tests.py valley_defs,zigzag_defs,rows_defs,deaths_defs,floor_defs,jam_defs,tower_defs 100000 2002 2100 res_defs.jsonl 60` | 18:19–18:55 |
| `res_merge.jsonl` | `run_tests.py valley,zigzag,valley_defs,wallmax,rows,wall 100000 3003 2100 res_merge.jsonl 60` | 18:19–18:54 |
| `res_negctl.json` | `negative_control.py` | 18:22 |
| `res_ties.jsonl` | `ties.py` | 18:56 |
| `res_dag.jsonl` | `run_tests.py dag,near45,dag 100000 4004 900 res_dag.jsonl 60 manyvar` | 18:56–19:11 |
| `res_manyvar.jsonl` | `run_tests.py valley,deaths,floor,rows,zigzag,floor_defs 100000 5005 900 res_manyvar.jsonl 60 manyvar` | 18:56–19:11 |
| `res_wallsweep.json` | `wall_sweep.py` | 18:56 |
| `res_exclab.json` | `exceptional_types.py` (420 s) | 18:57–19:04 |
| `summary_round2.txt` | `summarize.py res_dag.jsonl res_manyvar.jsonl res_ties.jsonl` | 19:14 |
| `summary_all.txt` | `summarize.py` of the six files `res_scaled`, `res_merge`, `res_defs`, `res_dag`, `res_manyvar`, `res_ties` | 19:14 |
| `memlog.txt` | run log with the free memory at the start of every run | |

**Outputs from earlier versions of the scripts.** `run_tests.py` and `generators.py` were last changed at 18:23 on 2026-10-06. By the file times, `smoke.jsonl`, `smoke2.jsonl`, `res_scaled.jsonl`, `res_defs.jsonl` and `res_merge.jsonl` (runs started 18:16–18:19) and `res_negctl.json` (18:22) were produced with earlier versions, which were not kept; `summary_all.txt` therefore combines records of both versions. With the published scripts, `res_negctl.json` is reproduced exactly; the recorded smoke run is not (its records differ from the first configuration on: the current driver checks more heights and more scale flows, and the random configurations differ from then on); the long runs were not repeated. All later outputs were produced by the published versions.

**Names in the recorded files** (left unchanged there):
- `_defs` in family names (`valley_defs`, …) and in the file name `res_defs.jsonl` = `_paper` in the scripts: the constants of the paper; "defs" referred to the definitions draft.
- `memlog.txt` gives the working names of the scripts: `p4sx_run.py` = `run_tests.py`, `p4sx_negctl.py` = `negative_control.py`, `p4sx_ties.py` = `ties.py`, `p4sx_wallsweep.py` = `wall_sweep.py`, `p4sx_exclab.py` = `exceptional_types.py`, `p4sx_summary.py` = `summarize.py` (and `p4sx_core.py` = `core.py`, `p4sx_gen.py` = `generators.py`).

## Changes made for publication

The logic of the tests is unchanged. Changes: file names and imports; docstrings and comments use the numbering of the paper; the family suffix `_defs` was renamed `_paper` (also in the driver's test for it); outputs of the helper scripts go to `out/`, and `run_tests.py` creates the folder of its output file; `exceptional_types.py` takes its time budget (default 420 s, as recorded) as an optional argument.
