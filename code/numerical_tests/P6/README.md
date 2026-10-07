# P6: numerical test of Section 8 (the path version of [R, Lemma 4.14])

This is a test on small analogues of the construction. It is **not part of the proof**. The programs were
written from scratch for this test on 2026-10-06, against the research draft of Sections 3 and 8: they
implement the definitions independently and do not use any other code. The statement numbers below are those of
the paper. Coordinates are exact rationals (rotations by Pythagorean unit vectors) and all
comparisons are exact; only angles θ are evaluated with mpmath at 40 digits.

Results recorded on 2026-10-06 are in `results/original/`. Re-runs made on 2026-10-07 for this release are in
`results/rerun_2026-10-07/`.

## What is tested

The master flow F⁰ (Definitions 3.9, 3.10, 3.14 and 3.20; M is decided before T) and the ceiling flow
(Definition 3.24) are traced on random and hand-built packings. Then the following checks run for every pair
e = (X, Y) (Definition 8.2) of both families:

| label in the outputs | statement of the paper | what is checked |
|---|---|---|
| `L4` | Lemma 8.5 | w_X is injective on P_e, and \|I_e\| ≥ m_e |
| `L5` | Lemma 8.6 | (i) the hull endpoints lie in top(X) and bot(Y), and 0 ≤ t(w) ≤ δ; (ii) the open set O_e meets no square (exact polygon test against every nearby square) |
| `L6` | Lemma 8.7, Remark 8.8 | (a) the end points c(x′) lie on bot(Y); the slit between the two end segments meets no square (exact); at five abscissae of J′_e a vertical ray from b(x′), the highest point of X on that line, first meets Y, on bot(Y), at distance exactly t(w)/χ_e; g^v_e ≤ δ/χ_e (and < 1.0000016 δ with the parameters of the paper); (c) see the next section; `exc`: the excluded end interval is at most cos φ_X · Δ_e |
| `L7` | Lemma 8.9 (a), (c) | no square is a member of a floor pair and of a ceiling pair; the slits S_e of both families are pairwise disjoint (exact test of every two slits whose x-ranges overlap) |
| `chain` | Step 1 of the proof of Theorem 8.3 | Σ_{i≤j} min(θ(Z_{i−1}, Z_i), β̄) ≥ min(β̄, a(Z_j)) at every R1-passed entry (tolerance 10⁻³⁰ for the 40-digit angles) |
| `pw` | Definition 3.14 | an independent point-by-point tracer agrees with the flow tracer at sampled points, and every R1 winner has the lexicographically smallest (g, x) |

The tracer also checks itself: the terminal measures add up to k, entries are processed in an acyclic order, no
path reaches the top side first, branches have slope at least 1, and paths that merge arrive from different
directions. Any failure is recorded as a violation.

Not tested here: Lemmas 8.10–8.12, Corollary 8.13 and the value of A′, and the results quoted from [R].

### The end term of Lemma 8.7(c): tested only in a weaker form

Lemma 8.7(c) states |J′_e| ≥ cos φ_X · (m_e − Δ_e) ≥ cos β̄ · m_e − 1.0000016 δβ̄. (The draft against which the test was written, and the first version of the paper text, had 1.0000017, a slightly larger constant; the proof gives 1.0000016.)

- The first inequality is checked exactly (`L6_sharp`, margin `L6_sharp_min`).
- The second inequality was tested only in the weaker form |J′_e| ≥ cos β̄ · m_e − **2.1** δβ̄ (`L6_weak_form`,
  margin `L6_weak_min`, in units of δβ̄), for every parameter set. The paper states it with **1.0000016** in
  place of 2.1, and that constant was not tested directly.

Arithmetic on the recorded margins (not a separate test): the two forms differ by exactly 1.0999984 δβ̄. On the
records with the parameters of the paper (δ = 10⁻⁵, α_max = 1.5·10⁻⁶; campaigns 2, 6, 7 and the 288 targeted
cases), the smallest recorded weak-form margin is 1.1407 δβ̄ in the campaigns and 1.1134 δβ̄ in the targeted
cases. So the form with 1.0000016 would have held in these cases with margin at least 0.0134 δβ̄.

In the scaled analogues, β̄ is replaced by the tilt threshold α_max of the flow.

## Files

| file | content |
|---|---|
| `p6x_core.py` | exact implementation of the flows (ray steps, priority D > W > contact > H, R1, R2, M before T; ceiling flow by reflection) and of the independent point tracer |
| `p6x_checks.py` | the checks L4–L7, chain and the R1 audit |
| `p6x_configs.py` | generators of exact configurations |
| `p6x_run.py` | campaign driver: generates configurations, traces floor and ceiling flows, runs all checks, writes one JSON line per configuration |
| `p6x_targeted.py` | 288 hand-built cases aimed at the bound of Remark 8.8 (δ = 10⁻⁵, α_max = 1.5·10⁻⁶), and 60 merge configurations whose cap square is a tilt terminator, run with M before T, with T before M, and with R1 replaced by "leftmost wins" |
| `p6x_negctl.py` | negative controls: the checks must fire when R1 or R2 is removed or when the constructions of Lemmas 8.6/8.7 are broken on purpose |
| `p6x_ties.py` | hand-built ties of positive measure (t_S = t_D, cumulative gap exactly δ, squares stacked up to h, squares touching the walls) |
| `p6x_summary.py` | adds up the JSON-lines output of `p6x_run.py` |
| `results/original/summary_all.txt` | totals over the seven campaigns (2026-10-06) |
| `results/original/campaigns_raw.zip` | the seven campaign records `camp_*.jsonl` (6,016 configurations) and the console logs of campaigns 1–4 (unchanged) |
| `results/original/targeted_out.json`, `negctl_out.json` | outputs of `p6x_targeted.py` and `p6x_negctl.py` (2026-10-06) |
| `results/rerun_2026-10-07/ties_console.txt` | console output of `p6x_ties.py`. Its output was not saved in 2026-10-06, so this is the re-run. |
| `results/rerun_2026-10-07/negctl_out.json` | `p6x_negctl.py` re-run with the scripts of this folder (see "Reproducibility") |

The scripts write into a subfolder `out/`, which they create.

## Requirements and quick run

Python 3.12 with `mpmath` (numpy is not needed). Run the commands from this folder:

```
python p6x_run.py 3 25 120 quick_scaled_3.jsonl scaled
python p6x_summary.py quick_scaled_3.jsonl
python p6x_targeted.py
python p6x_ties.py
python p6x_negctl.py 160 40
```

Measured on 2026-10-07 on a Windows 11 laptop (Python 3.12.10, one process at a time), in two sessions. The
times depend on the load of the machine.

| command | time | peak memory | result |
|---|---|---|---|
| `p6x_run.py 3 25 120 ...` | 7–14 s | 22 MB | 25 configurations, `nviol 0` on every line. The records equal the first 25 of campaign 3 (see "Reproducibility"). |
| `p6x_summary.py quick_scaled_3.jsonl` | < 1 s | 12 MB | totals of the 25 configurations, `violations: 0` |
| `p6x_targeted.py` | 15–38 s | 24 MB | `sharp cases 288 invalid 0 with viol 0`, max exc ratio 0.99990, min weak-form margin 1.11343; merge rows 60, 0 violations in all three variants. The output is identical to `results/original/targeted_out.json`. |
| `p6x_ties.py` | < 1 s | 22 MB | `total violations 0` (five cases; identical to `results/rerun_2026-10-07/ties_console.txt`) |
| `p6x_negctl.py 160 40` | 18–41 s | 22 MB | base 0/40 configurations with a violation; the checks fire without R1 or R2. Identical to `results/rerun_2026-10-07/negctl_out.json`. |

On a clean Linux machine (2026-10-07, `../rerun_linux_2026-10-07/`) every row was reproduced; in `targeted_out.json` three floating-point values differ in the last digit (relative 1.1·10⁻¹⁶).

Usage of the campaign driver: `python p6x_run.py <seed> <n_configs> <time_limit_s> <out.jsonl> [scaled|tiny]`.
Configuration i uses `random.Random(seed*100000 + i)`. Seeds below 4 use k ∈ {6, 8, 10, 12}; seeds 4 and above
use k ∈ {10, 12, 16, 20}. `scaled` cycles through (δ, α_max) ∈ {(0.2, 0.1), (0.2, 0.03), (0.1, 0.05),
(0.05, 0.03), (0.02, 0.01)}. `tiny` uses δ = 10⁻⁵ and α_max = 1.5·10⁻⁶, the values of the paper.

## The campaigns of 2026-10-06

Each campaign ran with a large n_configs and stopped at its time limit. Campaign 1 was stopped by hand.

| campaign | seed, mode | k | configurations | time limit |
|---|---|---|---|---|
| `camp_scaled_1` | 1, scaled | 6–12 | 429 | stopped by hand after about 5 min |
| `camp_tiny_2` | 2, tiny | 6–12 | 1,748 | 2400 s |
| `camp_scaled_3` | 3, scaled | 6–12 | 2,596 | 2400 s |
| `camp_scaled_4` | 4, scaled | 10–20 | 1,193 | 2700 s |
| `camp_scaled_5` | 5, scaled | 10–20 | 18 | 45 s |
| `camp_tiny_6` | 6, tiny | 10–20 | 16 | 45 s |
| `camp_tiny_7` | 7, tiny | 10–20 | 16 | 45 s |

Configuration kinds (about 500 each, 1,000 for `merge` and `bridge`): random piles, lattices, jams, columns of
stacked tilted squares, valleys, layered valleys that force merges (`merge`), fans, L-shaped strips, tilts near
45° (`label45`; exact 45° is not reachable with rational rotations, so π/4 ± 10⁻⁴ is used in the targeted
cases), and squares bridging two supports over a hole (`bridge`, which produces holes in I_e). 4,236
configurations use the scaled parameters and 1,780 those of the paper.

## Results recorded on 2026-10-06

No violation in any check (`results/original/summary_all.txt`, `targeted_out.json`):

- 6,016 configurations, 325,530 squares, 153,814 pairs (71,741 with X a square).
- L4: no overlap of exit coordinates; |I_e| > m_e strictly for 66,877 pairs.
- L5: 153,814 hulls, 435,449 exact square tests.
- L6: 153,814 pairs, 768,495 vertical ray tests, 383 pairs with holes in I_e. The largest ratio `exc` (bound 1)
  is 0.99886 in the campaigns and 0.99990 in the targeted cases, so the bound of Remark 8.8 is nearly sharp.
  The smallest `L6_sharp` margin is 0 (equality cases).
- L7: 153,699 slits, 70,391 exact pairwise tests (69,052 of them in campaigns 4–7; see below).
- chain: 156,968 entries; smallest margin −1.4·10⁻⁴², an equality case within the 40-digit arithmetic.
- Point tracer: 622,967 sampled points (8,869 at merges), no disagreement.
- Merges: 1,918 families with merge measure M > 0.
- Negative controls (40 configurations each), configurations with a violation: base 0, without R1 14, without R2
  30, exceptional end not trimmed 26, hull widened 40, T before M 0. T before M also gives 0 in the 60 targeted
  merge configurations.
- `p6x_ties.py`: no violation. The script also defines a sixth construction (`contact_exactly_at_h`). Its own
  guard skips it because the squares would overlap, so a contact at exactly height h is not exercised.

The fields `cat_*` in `summary_all.txt` (counts of special situations) exist only in the records of campaigns 5–7.
`merge_meas` and `merge_at_T_meas` exist only from campaign 3 on.

## Reproducibility

The records were compared field by field with re-runs made on 2026-10-07 with the scripts of this folder. Run
times were ignored, and the renamed key (last section) was mapped to its old name.

| record | re-run | result |
|---|---|---|
| campaigns 5, 6, 7 | all 18, 16, 16 configurations | identical |
| campaign 4 | first 19 configurations | identical in every recorded field. The current scripts add the counters `cat_*`. |
| campaign 3 | first 25 configurations (the quick run) | identical except `L7_tests` and the added `cat_*` counters |
| campaigns 1, 2 | first 25 and 20 configurations | **not reproducible**: different configurations (N and W differ) |
| `targeted_out.json` | full run | identical |
| `negctl_out.json` | full run | same qualitative result, **different counts** (e.g. "exceptional end not trimmed": 31/40 instead of 26/40) |
| `summary_all.txt` | `p6x_summary.py` over the seven records, in alphabetical order of the file names | every number equal |

Reasons, from the timestamps of the scripts and of the runs:

- Campaigns 1 and 2 and the negative controls ran with **earlier versions of the generator and of the tracer**.
  These were changed a few minutes after the campaigns started ("wall snap" and merge statistics). They cannot be
  reproduced exactly with the current scripts. Their records remain valid tests of the same statements, and
  campaign 2 is the largest campaign with the parameters of the paper.
- Campaigns 1–3 ran with an earlier version of the L7 check, which tested fewer pairs of slits (`L7_tests` adds up
  to 8, 1,317 and 14 in these campaigns). The exact test of every two slits whose x-ranges overlap was used
  from campaign 4 on.

## Names in the recorded files

- `L6_defs_min`, `L6_defs_form` (campaign records, `negctl_out.json`) and `defs` (`targeted_out.json`) are the
  weak form of Lemma 8.7(c) described above. In the scripts they are now `L6_weak_min`, `L6_weak_form` and
  `weak`. "defs" referred to the definitions file of the research draft. `p6x_summary.py` skips this field in
  the old records, and its minimum is given above.

## Limitations

- These are small analogues (k ≤ 20). The scaled parameters (δ up to 0.2, α_max up to 0.1) lie outside the
  range of the paper. The parameters of the paper are used in campaigns 2, 6, 7 and in the targeted cases.
- Sampled points test the point tracer, so sets of measure zero are seen only through the exact checks.
- The tests check the conclusions of the lemmas on the traced flows, not their proofs.
