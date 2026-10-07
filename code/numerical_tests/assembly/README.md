# Assembly: end-to-end numerical test of the assembled chain (Sections 9–10)

This is a test on small analogues of the construction. It is **not part of the proof**. For one packing at a
time it traces the master flow, the scale flows and the ceiling flows. It then evaluates every link of the
chain that leads from the per-scale count (Section 9) to the integration inequality (Theorem 10.3) with
**measured** quantities. The code was written from scratch for this test on 2026-10-06, against the research
draft of these sections, and uses no other tracer. The statement numbers below are those of the paper. It uses
float64 arithmetic: an inequality counts as violated when it fails by more than 10⁻⁹. An independent point
tracer cross-checks the flow tracer.

Results recorded on 2026-10-06 are in `results/original/`. Only a summary of the three campaign records is
included (they have 27 MB).

## What is tested

| link | statement of the paper | quantities (labels in `summary_ABC.txt`) |
|---|---|---|
| (1) | Theorem 9.1 and the steps of its proof; Theorem 9.2(a) | `L1_r_paper` = \|Z(s)\| / ((εs + b_W)(N(s) + Λ\*(s))) with Λ\* as in the paper; `L1_r_meas` uses the measured wall loss; `L1_r_sharp` uses the integrated measured losses (the sharpest form). Step 1: Z ⊂ ℝ × V_s (`VIOL_outV`). Lemma 9.3(c): a(Z) < β̂((1−ε)s) (`VIOL_tiltv`). Step 2: no alive contact on Z (`VIOL_p2v`). Step 3: Sh ≥ Σ c_Z (`L1_step3min`). Theorem 6.9(f) pointwise (`L1_p4min`), Lemma 9.7: L_T ≤ N (`L1_LTminusN`), Proposition 6.10 (`L1_lam_excess`, `L1_wall_ratio`). Theorem 9.2(a): `L1_lb1`, `L1_lb2`, `L1_lbZ`. |
| (2) | Lemma 9.6(c); Step 8 of the proof of Theorem 10.3 (first use of Tonelli) | `L2_symd` = measure of {T_s-terminated} Δ {T\* ≤ s}; N(s) nondecreasing (`VIOL_L2_mono`); `L2_T1_exact`, `L2_T1_grid` (Step 8) |
| (3) | Step 10 (second use of Tonelli, the constant invF) and the per-family chain of Steps 6–10 | `L3_T2_ratio` (with f = (1 − ω₀ − E)₊ as in the paper), `L3_T2s_ratio` (stress form f\* = the sum of the short chords), `L3_T2_swap_rel`, `L3_LBint_ratio`, `L3_P7int_ratio_*`, `L3_chain_ratio_*` |
| (4) | Steps 2–5 and 12 (summation) | `L4_ratio4_paper` / `L4_ratio4_sharp` = (1 − ω₀) ℓ(H) / (sum of the four contributions), with the paper's terms or with the measured kinds (i) and (ii); `L4_eq10_1_ratio` (the bound (10.1) for ℓ(H)); floor and ceiling separately (`L4_floor_*`, `L4_ceil_*`); consistency of the reflection (`L4_*_reflect_diff`) |

All ratios must be ≤ 1. The constants of the lemmas are replaced by the values measured on the same packing:
`K_A_ii_max`, `K_Aprime_meas_max`, `K_CLam_meas_max`, `K_Gamma_meas_max`, the actual E and the actual wall
loss. Both families (floor and ceiling flows) are evaluated.

### The scaled analogue

k = 16–32, δ = 0.02–0.2, c₀ = 0.08–0.25, y₀ ∈ {2, 3, 4}, several (c₁, ε, ω₀) per packing (`asmx_run.py`). The
widths use (sec α(s) − 1) q_y exactly. The number 1.0001 (the vertical extent of a square of Z(s)) is replaced by
h = cos b\* + sin b\* with b\* = c₁((1−ε)y₀)^{−3/4}, and b_W = 2.0002 by 2h. w₀ is recomputed accordingly, and
Ω_W is computed by quadrature.

## Files

| file | content |
|---|---|
| `asmx_core.py` | piecewise-affine flow tracer (floor interval split into pieces of constant combinatorial type), F_s, T\*, line evaluation |
| `asmx_chain.py` | parameters, line structure (kinds (i)–(iii), Y_b, Y_t), per-family flow data, links (1)–(4) |
| `asmx_gen.py` | configuration generators (`rows`, `nested`, `columns`, `lshape`, `valley`, `floortouch`, `random`, `hill`, `nestedEG`), including `gen_hill_v1`, the first version of `hill` |
| `asmx_run.py` | campaign driver (usage and the original commands at its top), with two negative controls per configuration |
| `asmx_summary.py` | aggregates campaign outputs (produces the format of `summary_ABC.txt`) |
| `asmx_targeted.py` | T1: negative control "T before M"; T2: independent point tracer (Cyrus–Beck) cross-check; T3: stress of Theorem 9.1 on a fine s-grid; T4: stress of the wall loss with `hill` configurations |
| `asmx_t2b.py` | T2b: point-tracer cross-check of merges and R1 winners on merge-heavy configurations |
| `asmx_flagged.py`, `asmx_replay.py`, `asmx_replay2.py`, `asmx_replaycheck.py` | replay tools: list the cases with diagnostic `r_shadow` > 1, regenerate them from the campaign seed, and re-integrate with 3000 points |
| `results/original/summary_ABC.txt` | summary of campaigns A, B, C (`asmx_summary.py`) |
| `results/original/targeted_T1.json` … `targeted_T4.json`, `targeted_T2b.json` | outputs of the targeted checks |
| `results/original/flagged.json`, `replay_campA.json`, `replay_campC.json`, `replay2_campA_0_1.json`, `replay2_campA_1_14.json` | the flagged cases and their replays |
| `results/original/resource_log.txt` | start and end times, memory, and the low-memory waits of the runs of 2026-10-06 |

The scripts write into a subfolder `out/`, which they create. The campaign records (`campA.jsonl`,
`campB.jsonl`, `campC.jsonl`; 1,346 records, 27 MB) are not included. `asmx_run.py` regenerates them (see
"Reproducibility").

## Quick run

Python 3.12 with `numpy`. From this folder:

```
python asmx_run.py quickC 13 40 150 valley,hill,columns,floortouch,rows,lshape 20,24,28
python asmx_summary.py quickC.jsonl
python asmx_targeted.py T2
python asmx_t2b.py
python asmx_replay.py campC 13 valley,hill,columns,floortouch,rows,lshape 20,24,28 results/original/flagged.json
```

Measured on 2026-10-07 on a Windows 11 laptop (Python 3.12.10, numpy 2.5.3, one process at a time):

| command | time | peak memory | result |
|---|---|---|---|
| `asmx_run.py quickC 13 40 150 …` | 30 s | 146 MB | 40 configurations, identical to the first 40 records of campaign C |
| `asmx_summary.py quickC.jsonl` | < 1 s | 15 MB | 238 variants, 9,520 cases; every `VIOL_*` counter is 0; the negative controls fire (`neg_w0_r_gt1` 5, `diag_nobW_gt1` 13) |
| `asmx_targeted.py T2` | 2 s | 30 MB | identical to `results/original/targeted_T2.json` (16,000 of 16,000 points agree) |
| `asmx_t2b.py` | 2 s | 29 MB | identical to `results/original/targeted_T2b.json` |
| `asmx_replay.py campC 13 …` | 9 s | 35 MB | identical to `results/original/replay_campC.json` |

On a clean Linux machine (2026-10-07, `../rerun_linux_2026-10-07/`) every row was reproduced up to floating-point last digits: in `targeted_T2b.json` one value differs by 2.6·10⁻¹² (relative), and in the 40 records of campaign C some values near 0 (minimal margins of order 10⁻¹⁵) differ by at most 3.6·10⁻¹⁵.

The memory guard of `asmx_run.py` (it waits until 1.5 GB are free) works only on Windows and is skipped elsewhere.

## The campaigns of 2026-10-06

| campaign | command (n_configs was large; each campaign stopped at its time limit) | valid configurations |
|---|---|---|
| A | `asmx_run.py campA 11 100000 1200 --hill-v1` (all eight families of `GENS`, k ∈ {16, 20, 24}) | 563 |
| B | `asmx_run.py campB 12 100000 1200 nestedEG 24,28,32 --hill-v1` | 420 |
| C | `asmx_run.py campC 13 100000 1000 valley,hill,columns,floortouch,rows,lshape 20,24,28` | 363 |

The seeds were not written down in 2026-10-06. They were recovered on 2026-10-07 by replaying the first
records: for seeds 11, 12, 13 the first four configurations reproduce k and the number of squares. The
families and k values can be read off the records, and the time limits off `resource_log.txt`.

## Results recorded on 2026-10-06

From `summary_ABC.txt`: 1,346 packings, 7,610 parameter variants, 304,400 cases (variant, family, s), with
Z(s) nonempty in 52,901. **No violation in any link**, with one diagnostic artefact explained below.

| link | worst ratio (≤ 1 required) |
|---|---|
| (1) Theorem 9.1 with Λ\* of the paper / with measured wall loss / sharpest form | 0.8996 / 0.9147 / 0.9147 |
| (1) Theorem 9.2(a): lower bound / \|Z(s)\| | 0.0012 |
| (1) Step 1, Lemma 9.3(c), Step 2: squares outside V_s, too tilted, or alively contacted | 0, 0, 0 |
| (1) Proposition 6.10: wall loss / 2(s + 2) tan α(s) | 0.170 |
| (2) Lemma 9.6(c): symmetric difference | ≤ 4.6·10⁻¹⁴ (equality) |
| (2) Step 8 | 0.574 |
| (3) Step 10, with f / with the stress form f\* | 0.862 / 0.9005 |
| (3) per-family chain, with f / with f\* | 2.8·10⁻⁴ / 0.159 |
| (4) Step 12 with the paper's terms / with measured kinds (i), (ii) | 0.0131 / 0.4605 |
| (4) bound (10.1) for ℓ(H) | 0.891 |
| (4) floor and ceiling separately | 2.8·10⁻⁴ |

- Tools: total measure error ≤ 2.1·10⁻¹⁴, overlap Ov ≤ 2.4·10⁻¹⁵. The independent point tracer agrees at 16,000 of
  16,000 points (T2), at 18,000 of 18,000 winner points and at 610 of 610 merge points (T2b).
- Measured constants (the values that replace the constants of the lemmas, much smaller than in the paper because
  the analogues are small): A_ii ≤ 0.0053, A′ ≤ 0.066, C_Λ ≤ 0.139, Γ ≤ 0.00095.
- **Diagnostic `r_shadow`** (`VIOL_L1_r_shadow 18`). This is \|Z(s)\| divided by the integral of the shadow Sh over
  V_s, computed with 24 sample heights. It exceeded 1 in 18 (variant, side) cases. The replays recompute these
  integrals with 3000 midpoints (`replay_campA.json`, `replay_campC.json`, `replay2_campA_*.json`). The largest
  value is then 1.0000049 (campA-131), and there the same quadrature gives Σ_Z ∫ c_Z = 27.99986 instead of the
  exact 28. So the excess is quadrature error, not a violation. The pointwise inequality Sh(y) ≥ Σ_Z c_Z(y),
  which is what Step 3 uses, holds at every sample (smallest margin −7·10⁻¹⁵, rounding).
- Negative controls (the checks can fail):
  - Dropping the fractional condition of H (w₀ → 0) gives 76,606 alive contacts on squares of Z(s) and ratios up
    to 10.5.
  - Deciding T before M breaks {T\* ≤ s} = T_s in 18 of 445 packings with merges. In T1 this happens in 6 of 12
    runs, exactly by the merge measure at the T_max squares; with M before T it holds in all 12.
  - Without b_W the ratio of (1) exceeds 1 in 3,136 cases (up to 4.23). With half of b_W it reaches 1.30 (T3).
- Lines of kind (iii) with f > 0 occur only in thin strips in these analogues. The chain was therefore also
  stressed with f\*.
- The wall factor 2(s + 2) tan α(s) was reached only up to 17%, so its sharpness is not shown.

## Reproducibility

On 2026-10-07 the records were compared field by field with re-runs of the scripts in this folder, ignoring run
times and memory.

| record | re-run | result |
|---|---|---|
| campaign C | first 40 configurations | identical |
| campaign A | first 8 configurations, with `--hill-v1` | identical |
| campaign A | the same, with the current `hill` generator | **different** in the configurations that use `hill` |
| campaign B | first 6 configurations, with or without `--hill-v1` | identical |
| `targeted_T1` … `T4`, `targeted_T2b` | full runs | identical |
| `replay_campC.json` | `asmx_replay.py` | identical |
| `replay2_campA_0_1.json`, `replay2_campA_1_14.json` | `asmx_replay2.py` (stages `genpart` 0–523 and `eval`) | identical |
| `summary_ABC.txt` | `asmx_summary.py` over the three original records | identical (one label renamed, see below) |

**Campaigns A and B ran with the first version of the generator `hill` (`gen_hill_v1`).** The generator was
rewritten after they started. To reproduce them exactly, use `asmx_run.py … --hill-v1`; the replay tools do this
for them automatically. Without it, campaign A cannot be reproduced exactly: every configuration that uses a
`hill` design differs. Campaign B uses only the family `nestedEG`, which does not call `hill`, so its records come
out the same with either version. Campaign C started after the rewrite.

Configurations are generated from one random stream per campaign. A configuration whose flow tracing hit the
60-second limit of `FamilyFlows` would differ on a slower machine. This did not happen in the re-runs above.

## Names in the recorded files

- `L4_41_ratio` in `summary_ABC.txt` is `L4_eq10_1_ratio` of the current `asmx_summary.py`: "41" referred to
  equation (4.1) of the research draft, which is (10.1) of the paper.
- `l35_min` is the slack of [R, Lemma 3.5]. `P7` refers to Section 9 of the paper (P7: the per-scale count), and
  `worstP7` to Theorem 9.1.

## Limitations

- The real parameter range (y₀ ≈ 4·10¹⁰, angles ≈ 10⁻⁶, δ = 10⁻⁵) cannot be reproduced. The test checks how the
  links compose and that the definitions are consistent with one another.
- The constants that enter from Sections 7 and 8 (C_Λ, A′) and from [R] (Lemmas 4.10, 4.13, 4.15) are replaced by
  measured values. Their bounds are not tested here (see the folders P6 and P8).
- Float64 arithmetic with tolerance 10⁻⁹, cross-checked by an independent point tracer.
