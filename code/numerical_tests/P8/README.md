# P8: numerical test of the constants (Sections 10–12)

This is an independent recomputation of the integration inequality, the explicit bound Ψ_π(k) and the
constants of Sections 10–12. It is **not part of the proof**. It was written from scratch from the formulas
of the research draft of these sections. It does not import or run any other code, in particular not
`code/k14_constants_check.py`, which is the script of the paper itself (Section 12.7). Arithmetic: mpmath at
50 digits, and 50-digit interval arithmetic (`mpmath.iv`, outward rounding) for every decisive inequality.

Results recorded on 2026-10-06 are in `results/original/`, unchanged.

## What is tested

| script | content |
|---|---|
| `p8x_core.py` | Definitions 3.1 (α, β̂, b\*, w₀, h₀), 3.2 (constraints (C0)–(C7), (Q1), (Q2)) and 3.6 (A, A′, C_Λ, invF, Γ, Ω_W, K_W, B̄). Theorem 10.3: the bounds (10.1) for ℓ(H) and (10.2) for Ω_W, and Ψ_π(k) of (10.3). Proposition 11.2: D₁ and the derivative of Ψ_π(k) − c\* k^{1/4}. Constants with K = 9 and K = 13, exact (`exact9`, `exact13`) and rounded up as in Sections 11.1 and 11.5 (`round9`, `round13`). |
| `p8x_main.py` | Table 1. For the parameter sets of Theorem A (Section 11.3), Theorem A13 (Section 11.5), Table 3 and Table 4: all constraints in interval arithmetic, conditions (i)–(iii) of Proposition 11.2, and the lower bounds for W. Independent quadratures of (10.1), (10.2), Γ (Lemma 10.5) and invF (Step 10 of the proof of Theorem 10.3). The exact ℓ(H) via the Hurwitz zeta function. 300 random small exact examples of (10.1). The asymptotic constants (Section 11.4). The integer claims of Corollary A (Section 11.6, Remark 11.7). The comparison with the bounds of [R] quoted in Section 1. 623 checks in all. |
| `p8x_scan.py` | "Ψ_π(k) ≥ c\* k^{1/4} for all k ≥ k₀" (Proposition 11.2): a log grid from k₀ to 10⁴⁰, all integers from k₀ − 50 to k₀ + 3000, and a rigorous interval cover of [k₀, 10¹⁰⁰]. Beyond 10¹⁰⁰ the monotonicity argument of Proposition 11.2 applies. |
| `p8x_asym.py` | The asymptotic constant: (1) closed form of c∞; (2) an independent 2-D maximisation of the limit functional (Lemma 11.4); (3) the family π_n of Proposition 11.5; (4) Proposition 11.6: (4a) sup of ρ(π) at fixed y₀, (4b) random admissible (π, k); (5) c₀ = 0.4472 against (C1); (6) an exploratory probe of the thresholds k₀. |
| `p8x_asym_fix.py` | Follow-up to (4a); see "The line (4a) 'False'" below. |
| `p8x_sens.py` | Slack of Theorems A and A13 against changes of C_Λ (Corollary 7.2), A′ (Corollary 8.13), A ([R, Lemmas 4.13 and 4.15]), b_W (Theorem 9.1), or an extra loss C′W/δ in Theorem 9.1, at fixed parameters (compare Section 12.6). |

## Results recorded on 2026-10-06

- `p8x_main_out.txt`: `TOTAL checks=623 fails=0 discrepancies=0 unsafe_roundings=11`. Every constraint holds in
  interval arithmetic for every parameter set. At π\*, k₀ = 4.62·10¹²: Ψ = 146.62429334 ≥ 0.1·k₀^{1/4} =
  146.60895355 (margin 0.01534). The quadratures agree with the closed forms. The exact sum over the full cells
  of H (Hurwitz zeta function), 5647.48557864514, is at least the integer-part version of (10.1), which is at
  least the bound (10.1) itself, 5647.48557862566. The 300 small exact examples of (10.1) pass. The lower
  bounds for W in Table 4 equal ⌈Ψ⌉.
- `p8x_scan_out.txt`: for Theorems A and A13 and the eight other pairs (c\*, k₀), no negative value on the grids.
  Ψ_π(k)/k^{1/4} is strictly increasing on the integers checked. The rigorous cover passes for every pair
  (155 pieces for Theorem A, 298 for Theorem A13). For fixed π\*, the crossing point of Ψ_π(k) = 0.1 k^{1/4} is
  k = 4.61667·10¹².
- `p8x_asym_out.txt`: c∞ = 0.16965241122 (exact constants; the closed form and the numerical maximisation agree
  to 4·10⁻⁴²), 0.16965029404 with the rounded constants. c∞^(13) = 0.14115933870 and 0.14115898336. In (4b),
  17,289 and 17,382 random admissible (π, k) give no violation of Ψ_π(k) < c∞ k^{1/4}.
- `p8x_asym_fix_out.txt`: see below.
- `p8x_sens_out.txt`: the relative slack in B̄ is 1.05·10⁻⁴ for Theorem A and 1.7·10⁻⁵ for Theorem A13. The
  threshold 4.62·10¹² survives C_Λ ≤ 1.0150, A′ ≤ 10.60905, A ≤ 10.60892 and b_W ≤ 53.09, each changed alone.
- `memlog.txt`: the memory checks and start/end times of the original runs.

### The 'stated' values are the research draft's values

In `p8x_main.py` and in its outputs, "stated" is the value **as displayed in the research draft of Sections
10–12**, against which the test was written. It is not always the value displayed in the paper. Each stated
value is compared with the 50-digit value: `OK` means correct and rounded in the safe direction, `DISC` means
off by more than one unit in the last digit (none occurred), and `UNSF` means correct but rounded to nearest
in the unsafe direction.

The 11 `UNSF` lines are draft values. **The paper displays all of them rounded in the safe direction**, following
its conventions in Section 12.1:

| quantity | draft (UNSF) | 50-digit value | paper |
|---|---|---|---|
| invF (Theorem A) | 0.99998819 (lower) | 0.9999881851… | ≥ 0.99998818 (Table 2) |
| bound (10.2) for Ω_W | 0.379415 (upper) | 0.3794151849… | ≤ 0.37941519 |
| 0.1·k₀^{1/4}, Theorem A | 146.60895 (upper) | 146.6089535… | ≤ 146.60896 |
| D₁, Theorem A | 0.0086430 (lower) | 0.0086429841… | ≥ 0.0086429 |
| 0.1·k₀^{1/4}, Theorem A13 | 215.83155 (upper) | 215.8315519… | ≤ 215.83156 |
| D₁, Theorem A13 | 0.0055791 (lower) | 0.0055790722… | ≥ 0.0055790 |
| B̄, Table 4, K = 9, k = 10³⁰ | 31.7249 (upper) | 31.7249406… | ≤ 31.72495 |
| B̄, Table 4, K = 13, k = 10¹³, 10¹⁴, 10¹⁶, 10³⁰ | 43.8461, 40.5842, 38.9055, 38.1280 (upper) | 43.846126…, 40.584226…, 38.905502…, 38.128048… | ≤ 43.84613, 40.58423, 38.90551, 38.12805 |

The decisive inequalities themselves are checked in interval arithmetic, so these roundings never affected a
conclusion. Every displayed value of the paper is also checked in its stated direction by
`code/k14_constants_check.py`. The pairs (c\*, k₀) with c\* = 0.05 and 0.08 in the outputs are from the research
draft and are not in Table 3 of the paper.

### The line (4a) 'False' in `p8x_asym_out.txt` is a withdrawn artefact

In the full run, (4a) reports at y₀ = 10³⁰ "sup rho ~ 0.1696602786 ... < c_inf ... : False", apparently
contradicting Proposition 11.6(b). **This line is withdrawn.** It comes from this script itself: its objective
is evaluated at 25 digits, and at ε = 1.28·10⁻²² the difference 1 − (1 − ε)^{3/4} cancels catastrophically.
This gives invF = 1.0000484 > 1, which is impossible because invF < 1 (Lemma 10.2). The paper is not
affected. `p8x_asym_fix.py` confirms this:

- at the reported point, invF = 0.99999998437 and ρ = 0.169650292527 < c∞ at 50 and at 80 digits, and in 50-digit
  interval arithmetic (only the 25-digit evaluation gives invF > 1);
- the maximisation, redone with a cancellation-free invF at 60 digits, gives sup ρ < c∞ at y₀ = 10¹², 10²⁰,
  10³⁰ and 10⁴⁰, confirmed with 120-digit intervals.

The note at the end of `p8x_asym_out.txt` records this. The line was kept in the record, unchanged.

## Quick run

Python 3.12 with `mpmath` (numpy is not needed). From this folder:

```
python p8x_main.py --quick
python p8x_scan.py --quick
python p8x_asym.py --quick
python p8x_asym_fix.py
python p8x_sens.py
```

The outputs go to `out/`. Measured on 2026-10-07 on a Windows 11 laptop (Python 3.12.10, mpmath 1.4.1, one
process at a time), in two sessions. The times depend on the load of the machine.

| command | time | peak memory | result, compared with `results/original/` |
|---|---|---|---|
| `p8x_main.py --quick` | 3–8 s | 24 MB | `TOTAL checks=623 fails=0 discrepancies=0 unsafe_roundings=11`. Every number equals the record except the number of random examples of (10.1) (30 instead of 300) and the run time. |
| `p8x_scan.py --quick` | 1–5 s | 20 MB | Theorems A and A13 only: crossing points, minima on the integers and the grid, and both rigorous covers (155 and 298 pieces, PASS) equal the record. |
| `p8x_asym.py --quick` | 5–24 s | 20 MB | (1)–(3), (5) and the rows y₀ = 10¹², 10²⁰ of (4a) equal the record. (4b) uses 2,000 samples, with no violation. The (6) probe ends at a slightly different optimum, with a positive interval margin as in the record. |
| `p8x_asym_fix.py` | 3–10 s | 20 MB | output identical to the record |
| `p8x_sens.py` | 1 s | 19 MB | every number equal to the record (labels now use the paper's numbering) |

The full modes (`p8x_main.py`, `p8x_scan.py`, `p8x_asym.py` without `--quick`) took 192 s, 665 s and 707 s
in 2026-10-06. They were re-run on a clean Linux machine on 2026-10-07 (`../rerun_linux_2026-10-07/`; 43 s, 46 s and 50 s there): the JSON outputs equal the records after mapping the renamed labels, and every number with four or more decimals in the text outputs equals the record. `--quick` reduces only sample sizes and the number of cases,
as described at the top of each script.

## Names in the recorded files

The records in `results/original/` use the numbering of the research draft. In the scripts here the labels use
the numbering of the paper:

| in the records | in the paper |
|---|---|
| `PROOF_P8`; `PROOF 4.1`, `PROOF 4.2`, `PROOF 2.2` | the research draft of Sections 10–12; its sections on the constants (now Section 12.2 with Table 1, Section 11.4, and the last paragraph of Section 12.2) |
| `DEFS formula` | an earlier formula for c∞ in the draft's definitions (c₀ = 1/√5), not used in the paper |
| `task: …`, `wall_task` | a value quoted in the assignment of this test (equal to the draft's value); `wall_task` is the wall term with the short value 1.358 |
| `Thm4.3`, `main` | Theorem A, explicit part (Section 11.3) |
| `Thm4.4`, `hand` | Theorem A13 (Section 11.5), the variant that does not use the computer-assisted Lemma 4.10 of [R] |
| `Lemma 4.2`, `lemma-4.2` | Proposition 11.2 |
| `(3.1)`, `(3.2)`, `(3.3)` | (10.1), (10.2), (10.3) |
| `Theorem 4.1(i)`, `Theorem 4.1(ii)` | Proposition 11.5, Proposition 11.6 |
| `corollary 4.5` | Corollary A |
| `v11` | [R], the previous paper (its bound quoted in Section 1) |
| `v1.1 0.033 ln k` | the constant of an earlier version of [R] |
| `(P5)`, `(P6)`, `(P7)` in `p8x_sens_out.txt` | Corollary 7.2 (C_Λ), Corollary 8.13 (A′), Theorem 9.1 (b_W) |
| `L4.13/L4.15` | [R, Lemmas 4.13 and 4.15] |
| `k0-minimality remark` | the remark at the start of Section 12 that the thresholds are not claimed to be the least ones |

## Limitations

- The test covers the integration and constants layer only: (10.3), Proposition 11.2, Theorems A and A13,
  Corollary A, and the asymptotic constant. The inputs from Sections 4–9 and from [R] (A, A′, C_Λ, b_W) are
  assumed.
- (6) in `p8x_asym.py` is an exploratory local optimisation. It finds no parameters that reach c\* = 0.1 at
  k = 4.61·10¹² or 4.60·10¹² (best ratios 0.0999793 and 0.0999478). This is not a proof that 4.62·10¹² is
  the least threshold, and the paper does not claim that it is.
