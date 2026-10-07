# Numerical tests for Section 4 (P1 and P2)

Exact-arithmetic tests of the statements of Section 4 of the paper: the height identity (Theorem 4.2, P1), its corollary (Corollary 4.4), quantization at bottom sides (Lemma 4.8, Theorem 4.11, P2), the reflection lemma (Lemma 4.7), the auxiliary inequality (Lemma 4.6) and a scaled analogue of Theorem 4.13 (P2 for the squares of Z(s)).

These are tests on small, scaled configurations (container side k between 8 and 134), not part of the proof. "No counterexample" below is the outcome of a search and is therefore provisional.

All geometry is exact: squares have rational centres and rational t = tan(φ/2), so all positions, gaps and contact points are fractions; only the thresholds α and a(Z) are evaluated with mpmath (60 digits). The flow is traced without rule R1 (every path that reaches an entry candidate goes on if a(Y) < α_T). A path of the flow with R1 is an initial part of the path of the same floor point without R1, so every first contact of the flow with R1 is among the contacts checked here.

## Files

| File | Content |
|---|---|
| `tracer.py` | exact tracer; independent re-verification of each contact (Cyrus–Beck clipping); the quantities of Corollary 4.4, Lemma 4.8 and Theorem 4.11 |
| `common.py` | output folder `out/`, memory self-check, exact builders of columns of squares |
| `run_families.py` | adversarial families (exact contacts, squares on the floor, overhangs, valleys, piles, towers) and a random-cluster fuzzer |
| `consequence.py` | scaled analogue of Theorem 4.13(b)–(d): forward sampling and exhaustive backward search |
| `side_contacts.py` | non-vacuity check for `consequence.py`: forbidden squares placed where live paths reach them |
| `merge_points.py` | exact merge points (the situation in which R1 decides) |
| `reflection.py` | Lemma 4.7(a)–(c) on 3,000 random rational squares |
| `thresholds.py` | the threshold a\* of Lemma 4.6 and other constants used in Section 4 |
| `witnesses.py` | explicit near-extremal configurations |
| `aggregate.py` | aggregates the `*_summary.json` files of a folder |
| `results/original/` | the recorded outputs (see below), unchanged |

Requirements: Python 3.12 and `mpmath`. Each script writes into `out/` next to the scripts.

## Quick checks

Measured on the author's Windows PC on 2026-10-07 (other programs were running at the same time):

| Command | Time, peak memory | Result |
|---|---|---|
| `python thresholds.py` | 5 s, 26 MB | identical to `results/original/thresholds_out.json` |
| `python reflection.py` | 3 s, 19 MB | identical: 3,000 squares, 0 mismatches |
| `python merge_points.py` | 8 s, 22 MB | identical: 42 exact merge points, 0 failures |
| `python witnesses.py` | 3 s, 21 MB | identical to `results/original/witnesses_out.json` |
| `python aggregate.py results/original` | 1 s, 20 MB | `out/aggregate_out.json` is byte-for-byte identical to the recorded `aggregate_out.json` |
| `python run_families.py 60 smoke quick` | 42 s, 25 MB | 216 configurations, 12,062 contacts, 0 failures (the time budget allows the families exact, floorZ and overhang) |
| `python run_families.py 120 quick_tvp quick tower,valley,pile` | 111 s, 32 MB | 332 configurations, 15,521 contacts, 0 failures |
| `python run_families.py 60 fuzz quick fuzz` | 46 s, 29 MB | 741 configurations, 13,395 contacts, 0 failures (same seed 463 as the recorded fuzz run) |
| `python consequence.py 90 conseq_quick quick` | 57 s, 25 MB | 6 configurations, 0 violations, backward search 14,688 points / 0 found, sanity 87/87 |
| `python side_contacts.py 100` | 102 s, 22 MB | the first 136 of the 288 cases, identical to the recorded ones; 0 violations |

"Identical" means equal after the renaming of keys described below, ignoring the field `elapsed`. The side contacts with forbidden squares (the point of `side_contacts.py`) occur only among the recorded cases 209–288 (c0 = 3/5), which the 100-second run does not reach. The full run `python side_contacts.py` took 2 minutes when it was recorded; at the speed of the quick check it would take about 3.5 minutes. It was re-run on a clean Linux machine on 2026-10-07 (61 s): all 288 cases are identical to the record (`../rerun_linux_2026-10-07/`), and so are the outputs of the other quick checks that are claimed identical above.

## What was tested, and in which form

- **Theorem 4.2 (P1).** The identity (4.3) is checked exactly at every recorded contact, and each contact is re-verified by an independent routine (Lemma 4.1 (F2), (F4)).
- **Corollary 4.4 and Lemma 4.8 were tested in a weaker form than stated in the paper.** The paper has the term (sec α − 1)·q_y. The scripts use α²·q_y/(2 cos α), which is larger, since sec α − 1 = α²/(2 cos α)·(1 − α²/12 + O(α⁴)). So the lower bound q_y − m ≥ −α² q_y/(2 cos α) of Corollary 4.4 and the width K + δ + f(a(Z)) of Lemma 4.8 were checked with K = α² q_y/(2 cos α) for α = α_T (the flow threshold) and for α = α_eff (the largest inclination actually traversed), and the inclusions (4.6), (4.7) with the α_eff version of this K. The recorded ratios (`ratio_tilt_T`, `ratio_tilt_eff`, `r_prover`, `r_g`, `r_eff`) are relative to this larger constant. The paper's form with sec α − 1 was not tested directly.
- **The form with 0.50001 α² q_y** (Corollary 4.4, second part, and Theorem 4.11(a)) was tested directly only when α_T ≤ 10⁻³ (keys `cor_defs`, `r_defs` in the recorded files), that is, for the parameter sets `a1e-3_d1e-5` and `theory_a1.5e-6_d1e-5`, not for the larger α_T of the other sets or of the fuzzer. This is the range of Theorem 4.11, which assumes α(s) ≤ 10⁻³.
- **Theorem 4.11(b).** The auxiliary width 0.50001 α² q_y + δ + 1.0001 a(Z) is only recorded, not required, because the paper asserts it only for a(Z) ≤ a\*. `thresholds.py` computes a\* = 2.000133357782697…·10⁻⁴ exactly (Lemma 4.6(c)) and builds configurations in which the conclusion with the auxiliary width fails for larger a(Z) (from about a = 4.69·10⁻³ at δ = 10⁻⁵).
- **Theorem 4.13, scaled analogue** (`consequence.py`, `side_contacts.py`): α(s) = c0·s^(−1/2) with c0 ∈ {3/10, 9/20, 3/5} and s ∈ {4, 9, 16, 25}; the constant w0 is built with α²/(2 cos α) in place of 0.50001 α², since α(s) exceeds 10⁻³ here. This is an analogue at small scale, not the parameters of the paper (k of order 10¹²).
- Not tested: the ceiling flow (the reflected piles, family `pile_refl`, produced no contacts, so that part was vacuous), squares at exactly 45° (impossible with rational t), the actual scale of the paper.

## Recorded results (`results/original/`)

Outcome of all recorded runs: no failure of any checked statement, no anomaly.

- `run_families.py` (4 runs): 29,098 configurations, 4,595,287 traced paths, 739,761 first contacts checked, of which 727,596 on closed bottom sides (31,167 at their endpoints); 0 failures of the independent re-verification, 0 failures of the checked inequalities, 0 anomalies. Largest ratios: `r_prover` 0.99997743, `r_g` 0.99999999, `r_eff` 0.9999999987, `ratio_tilt_eff` 0.9999999999992, `ratio_gap_delta` 0.9999982; `ratio_gap` reaches 1 (equality is allowed there).
- `consequence.py` (4 runs): 142 configurations with 1,404 squares of the Z(s)-analogue and 12,754 squares excluded by the bound of Lemma 4.8 (the two sets may overlap), 277,544 traced paths: no path met the closed bottom side of such a square or passed through one; the backward search from 622,608 points of bottom sides found no path; the backward search reproduced 4,746 of 4,746 forward contacts (sanity check).
- `side_contacts.py`: 288 cases (144 with a forbidden square, 144 controls); forbidden squares were met 46 times, always on a side and at the termination point of the path, never on the bottom side, never passed; the control squares were reached on their bottom sides 8,320 times (forward contacts plus backward-search hits).
- `witnesses.py`: at a vertex contact the ratio of Lemma 4.8 is 1 − 2.17·10⁻⁵ (case E1), so the width cannot be reduced by a constant factor; the cases with alternating tilts (`E3_alt_...`) did not reach their target.

| File | Produced by (command from `memlog.txt`; exact arguments partly reconstructed) | Run (2026-10-06) |
|---|---|---|
| `smoke_summary.json`, `smoke_viol.jsonl` | `run_families.py 60 smoke quick` | 17:27–17:28 |
| `thresholds_out.json` | `thresholds.py` | 17:31 |
| `conseq_quick_summary.json` | `consequence.py 90 conseq_quick quick` | 17:31 |
| `witnesses_out.json` | `witnesses.py` | 17:32 |
| `conseq_summary.json` | `consequence.py 1500 conseq` (the budget ended within c0 = 3/10) | 17:32–17:57 |
| `main_summary.json`, `main_viol.jsonl` | `run_families.py 1740 main` | 17:34–18:03 |
| `conseq_c06_summary.json` | `consequence.py` with c0 = 3/5, budget 1300 s | 17:35–17:56 |
| `reflection_out.json` | `reflection.py` | 17:57 |
| `main2_summary.json`, `main2_viol.jsonl` | `run_families.py 1800 main2 full tower,valley,pile` | 17:57–18:10 |
| `fuzz_summary.json`, `fuzz_viol.jsonl` | `run_families.py 1500 fuzz full fuzz` | 17:58–18:23 |
| `conseq2_c045_summary.json` | `consequence.py` with c0 = 9/20, budget 1450 s | 18:04–18:28 |
| `side_contacts_out.json` | `side_contacts.py 700` | 18:11–18:13 |
| `merge_points_out.json` | `merge_points.py` | 18:13 |
| `aggregate_out.json` | `aggregate.py` | 18:28 |
| `memlog.txt` | run log: start, end and free memory of every run | |

The `*_viol.jsonl` files are empty: no failure was written. Five further empty files, created only as a side effect of an import in the earlier versions of the scripts, are not included.

**Outputs from earlier versions of the scripts.** By the file times, these outputs were produced before the last change of a script they used: `smoke_summary.json`, `thresholds_out.json`, `conseq_quick_summary.json`, `witnesses_out.json`, `conseq_summary.json` (all before the last change of `tracer.py` at 17:34, of the column builders and of `run_families.py` at 17:46, and, for the consequence runs, of `consequence.py` at 17:50), `main_summary.json` (before the change of `run_families.py` at 17:46; it has no valley family) and `conseq_c06_summary.json` (before the changes at 17:46 and 17:50). The earlier versions were not kept. With the published scripts, `thresholds_out.json` and `witnesses_out.json` are reproduced exactly; the recorded smoke run (tower families only, one parameter set) and the recorded quick consequence run (47 instead of 62 squares of the Z(s)-analogue, other configurations) are not; the long runs were not repeated. The first attempt of the main run (17:32) stopped after 2 seconds with exit code 1; the run log calls the restart a "relaunch after 0/0 fix", and `tracer.py` was last changed at 17:34, just before the restart.

**Names in the recorded files** (left unchanged there; the scripts use the new names):
- `prover` in keys (`r_prover`, `w_prover`, `worst_r_prover`, `P2_prover`, `prover_holds`, `relint_r_prover`, `one_minus_r_prover`) = `alphaT` in the scripts: the width of Lemma 4.8 with α = α_T, as stated in the research draft of Section 4.
- `defs` in keys (`cor_defs`, `r_defs`, `P2_defs`) = `050001` in the scripts: the form with 0.50001 α² q_y; "defs" referred to the definitions draft.
- `memlog.txt` names the scripts by working labels: "extreme witnesses" = `witnesses.py`, "thresh" = `thresholds.py`, "conseq" = `consequence.py`, "refl lemma3.3" = `reflection.py` (3.3 was the draft number of Lemma 4.7), "side/non-vacuity" = `side_contacts.py`, "merge exact" = `merge_points.py`, "smoke", "main run", "main2", "fuzz" = `run_families.py`.

## Changes made for publication

The logic of the tests is unchanged. Changes: file names; the column builders and the memory self-check moved from the driver to `common.py` (the self-check also works outside Windows); outputs go to `out/`; docstrings and comments use the numbering of the paper; the keys `prover` → `alphaT` and `defs` → `050001` were renamed (values unchanged); `aggregate.py` takes the folder as an argument and carries over every `worst_*` key, which gives the same output on the recorded summaries.
