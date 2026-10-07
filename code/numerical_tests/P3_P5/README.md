# Numerical tests for Sections 5 and 7 (P3 and P5)

Exact-arithmetic tests of the geometric statements behind Proposition 5.10 (P3, end collisions) and Theorem 7.1 (P5, deaths, merges and gap overlap), on small configurations (container side k = 4 to 8; δ from 10⁻⁵ to about 0.5; most families outside the regime of the paper, some inside it). These are tests, not part of the proof; "no counterexample" is the outcome of a search and is provisional.

Squares have rational phases (s = tan(φ/2) rational), so the tracer of the flow F⁰ (rules D > W > contact > H, R1 and R2) works with exact fractions. For each run the analysis computes the gap regions of all paths exactly and checks:

- (a) at most two paths through a point, M_g ≤ 2 (Lemma 7.7); no point reached from the floor and from a square (Lemma 7.4); no overlap of two paths from the same source (Lemma 5.3);
- (b) at points reached by two paths: q = v₁ + a u + b n with 0 ≤ a ≤ w_e = δ tan θ / cos θ and 0 ≤ b ≤ δ, dist(X₁, X₂) < δ tan θ, and inequality (7.3) (Lemma 7.8);
- (c) deaths, (d) merges, (e) gap overlap: the bounds of Sections 7.3–7.5 that give Theorem 7.1(D), (M), (G);
- (f) end collisions: |E_{X,Y}| ≤ δ tan|Δ| for every source and target (Lemma 5.9(b), used for Proposition 5.10);
- conservation of measure D + W + E + M + T + H = k.

## Files

| File | Content |
|---|---|
| `tracer.py` | exact beam tracer of F⁰ |
| `analysis.py` | the checks (a)–(f) for one trace |
| `configs.py` | configuration families |
| `run_family.py` | runs one family (in chunks: `NEXT <index>` gives the start of the next chunk) |
| `summarize.py` | aggregates the run records `res_*.jsonl` of a folder per family |
| `degenerate.py` | targeted degenerate cases (ties of the rules, vertices on the floor, walls, contacts at the height bound) |
| `deaths.py` | near-extremal configurations for the death bound |
| `cross_check.py` | comparison with an independent per-point floating-point tracer |
| `mutations.py` | mutation tests: deliberately broken versions of the checks must be detected |
| `geom_search.py` | numerical search (differential evolution) for the boundary cases of Lemmas 7.5, 7.7 and 7.8 |
| `boundary_exact.py` | exact re-check (50 digits) of two boundary configurations found by `geom_search.py` |
| `results/original/` | the recorded outputs (see below) |

Requirements: Python 3.12 with `numpy` (`analysis.py`), `mpmath` and `scipy` (`geom_search.py`, `boundary_exact.py`). Outputs go to `out/`.

## Quick checks

Measured on the author's Windows PC on 2026-10-07 (other programs were running at the same time). Unpack `results/original/runs.zip` first, for example into `results/original/runs/`.

| Command | Time, peak memory | Result |
|---|---|---|
| `python summarize.py results/original/runs` | 1 s, 124 MB | `out/summary.json` is byte-for-byte identical to `results/original/summary.json` |
| `python degenerate.py` | 1 s, 28 MB | identical to `degen_results.json` (32 records; final line `BAD 2`, explained below) |
| `python deaths.py` | 1 s, 28 MB | identical to `death_results.json` (15 records) |
| `python mutations.py` | 7 s, 28 MB | identical to `mutation_results.json`: no flag without mutation; all four mutations detected |
| `python cross_check.py` | 6 s, 15 MB | identical to `cross_check_results.json` (18 configurations, largest deviation 0.0054) |
| `python boundary_exact.py` | 0.2 s, 19 MB | identical to `boundary_exact.json` |
| `python geom_search.py starmin 50 0.3` / `star 50 0.3` / `starmin 50 0.9` | 8–13 s, 73 MB each | identical to the recorded files of the same names |
| `python run_family.py valley 60` (likewise `smoke`, `bigvalley`, `mechA`, `three`, `fan`, `strip`) | 0.4–17 s, 30 MB | the whole family; identical to the recorded records (168, 1, 90, 48, 72, 144, 72) |
| `python run_family.py merge 50` | 42 s, 33 MB | the whole family (5,040 records), identical to the two recorded chunks |
| `python run_family.py near45 50` | 41 s, 29 MB | all 72 configurations; the 68 recorded records identical; the other 4 analysed for the first time (below) |
| `python run_family.py mrand 60` | 60 s, 36 MB | 1,598 records, identical to the first 1,598 recorded ones |
| `python run_family.py mrand 20 9000` | 20 s, 31 MB | identical to the 47 records of `res_mrand_9000.jsonl` |
| `python run_family.py rowjam 50` | 50 s, 37 MB | identical to the 143 records of `res_rowjam_0.jsonl` |
| `python run_family.py regime 40` / `tri 50` | 40 s / 50 s, 38 MB | identical to the 303 / 229 records of the first recorded chunk |

"Identical" means equal record by record, ignoring the timing fields (`sec`, `analysis_sec`) and with the renamed key of `boundary_exact.json` (below). On a clean Linux machine (2026-10-07, `../rerun_linux_2026-10-07/`) the same records were obtained, with two exceptions: some floating-point diagnostics differ in their last digits, and some configurations of `near45`, `tri`, `rowjam` and `mrand`, whose coordinates are computed in floating point, come out slightly different there (0 violations on them as well).

## Recorded results (`results/original/`)

`runs.zip` contains the 57 run-record files `res_<family>_<start>.jsonl` and 15 run logs `log_*.txt`, unchanged, except one log in which a path was redacted (below). The other files are stored as they were.

Over all families (`summary.json`): 16,677 run records, of which 12,517 analysed (4,160 configurations of the family `merge` were invalid and skipped); no error, no timeout; 5,166 runs with points reached by two paths, 2,874 with merges. **No violation of (a)–(f), no anomaly, conservation exact in every run.** Largest values: q_a/w_e 0.999002, b/δ 0.9999993, the ratio of (7.3) 0.9999999999, dist(X₁, X₂)/(δ tan θ) 0.987, |E_{X,Y}|/(δ tan|Δ|) 0.9999999998, gap overlap 0.874 of its bound, deaths 0.310 of their bound; at most one merge target per pair; M_g ≤ 2 everywhere (no point reached by three paths, also in the family `tri`, which looked for M_g = 3 near 45°).

- `degenerate.py`: in the configurations with a corner on the floor, the bound of Lemma 5.9(b) is attained (ratio exactly 1).
- `mutations.py`: halving w_e was detected in 14 of 25 configurations, switching R1 off in 25 of 25, "contact before death" in the tie configuration (H changes from 0 to 4/5), bounds evaluated with δ/2 in 25 of 25.
- `geom_search.py` and `boundary_exact.py`: the hypothesis δ ≤ ½ cos 2α_max of Lemma 7.5 is sharp. For θ = 0.6 the conclusion of Lemma 7.5 fails for the delay 0.4126682 > ½ cos θ = 0.4126678 (exact configuration with disjoint squares and clear rays), as in Remark 7.6. Inequality (7.3) fails for θ = 1.4 and δ = 0.45, far above ½ cos θ ≈ 0.085; in that configuration the ray from one source meets the other square, so it satisfies only the geometric hypotheses of Lemma 7.5, not those of a flow. In the regime of the paper (δ ≤ 10⁻⁵) the hypothesis holds with a wide margin.

**Facts to know about these records** (checked against the files):

1. **Largest q_a/w_e.** `summary.json` gives 0.999002 (family `rowjam`, run `rowjam_m1_s646`, δ = 1/100). An earlier summary of these tests quoted 0.9975, which is the largest value in the family `mrand` (run `mrand4716`). Both are below 1.
2. **169 rowjam runs from an earlier generator.** `res_rowjamU_0.jsonl` holds 169 runs (seeds 1–169) made with the first version of the rowjam generator (uniform gaps; its log is `old_log_rowjam_gapuniform.txt`). Their records carry the family name `rowjam`, so `summary.json` counts them in that family: 2,269 = 169 + 2,100, and seeds 1–169 occur twice. The published generator reproduces `res_rowjam_0.jsonl`, not these 169 records.
3. **near45: 68 of 72 configurations analysed.** The chunk that should have started at index 68 was run twice and ended both times with `NEXT 68` without a record. A chunk regenerates all configurations before its start index, and with the slow configuration generator this used up the 50-second limit before index 68 was reached (the traces themselves take 0.1–0.2 s). The published runner does the whole family in 41 s. It reproduces the 68 records, and the 4 missing configurations (δ = 1/10, seed 2) show M_g = 1 and no violation.
4. **mrand indices 8319–8999 were not run.** The last regular chunk stopped with `NEXT 8319`; the next run of the family was a 20-second memory probe starting at index 9000 (`res_mrand_9000.jsonl`, `memprobe_out.txt`). The run names are shifted by one (index 8319 = `mrand8320`).
5. **`BAD 2` of `degenerate.py` is a false alarm of the counting.** The two flagged records are the configuration D1 (gap exactly δ above an axis-parallel square A, for the two values of δ), whose expectation is `entries == 1`. But `entries` counts the beam pieces that enter a square, not the squares. The beam entering A is split at x = 1.2, the left edge of the square B above, so `entries` = 2. The measures show that B is never entered: D = 6 = k and H = 0, while a path entering B would cross the height bound inside B and end with H, as in D1b (gap slightly less than δ), where H = 4/5. A diagnostic run that printed every counted entry found two entries into A, for x in (1, 1.2) and (1.2, 2), and none into B. So the rule "death before contact" works as required; the expectation is wrong. It was left unchanged.

**Outputs from earlier versions of the scripts.** The scripts were last changed (2026-10-06) at 13:57 (`analysis.py`), 14:25 (`summarize.py`), 14:38 (`geom_search.py`), 14:39 (`boundary_exact.py`), 14:41 (`degenerate.py`), 14:42:53 (`configs.py`), 14:42:56 (`run_family.py`), 14:44 (`tracer.py`), 14:45 (`mutations.py`, `cross_check.py`) and 14:46 (`deaths.py`). By these times:
- the run records of `smoke`, `valley`, `bigvalley`, `merge`, `fan`, `three`, `mechA`, `strip`, `near45`, the 169 runs in `res_rowjamU_0.jsonl`, the rowjam chunks starting at 0–945 (seeds 1–1076) and the mrand chunks starting at 0–7508 were made with earlier versions of `run_family.py`, `configs.py` and `tracer.py`;
- the regime chunks starting at 0 and 303, the tri chunk starting at 0 and the rowjam chunk starting at 1076 with an earlier `tracer.py`;
- `degen_results.json` with earlier `tracer.py` and `configs.py`;
- `geom_two_all.json`, `geom_three_all.json`, `geom_starmin_0.01.json`, `geom_starmin_0.9.json` and `geom_starmin_1.4.json` with an earlier `geom_search.py`.

The earlier versions were not kept. Every recorded output from an earlier version that was rerun in the quick checks was reproduced exactly, except the 169 records of the earlier rowjam generator. The geometry searches `two` and `three` (13 and 41 minutes when recorded), `starmin` for θ = 0.01 and 1.4 and `star` for θ = 10⁻⁶, 10⁻³ and 1.4 were re-run on a clean Linux machine on 2026-10-07 (`../rerun_linux_2026-10-07/`) and equal the records; `three` produced 16 entries there with a larger time limit, and the 12 recorded ones are identical. Not rerun: the run chunks not listed in the quick checks.

**Names in the recorded files** (left unchanged there):
- `memlog.txt` gives the working names of the scripts: `xrun.py` = `run_family.py`, `xdegen.py` = `degenerate.py`, `xdeath.py` = `deaths.py`, `xcross.py` = `cross_check.py`, `mutate.py` = `mutations.py`, `geomsearch.py` = `geom_search.py`, with the modes `43` = `two` and `44` = `three`.
- The draft numbers 4.3 and 4.4 of Lemmas 7.5 and 7.7 also appear in the keys `lemma43` (`geom_two_all.json`), `lemma44` (`geom_three_all.json`) and `lemma43_fail` (`boundary_exact.json`; `two_sources_fail` in the script), and in the logs `log_g43.txt`, `log_g44.txt` in `runs.zip`.
- Renamed files (contents unchanged): `geom_43.json` → `geom_two_all.json`, `geom_44.json` → `geom_three_all.json`, `xcross_results.json` → `cross_check_results.json`.
- **Paths redacted:** in `old_log_rowjam_gapuniform.txt` (in `runs.zip`), one traceback line contained the absolute path of the working folder; that prefix was replaced by `<path>`. Nothing else was changed.
- The cache `cache_rowjam.pkl` of the recorded runs is not included; `run_family.py` builds its own cache in `out/`.

## Changes made for publication

The logic of the tests is unchanged. Changes: file and module names (the analysis module is `analysis.py`, so that its name does not clash with the function `analyze`, which `mutations.py` patches through the module); docstrings and comments use the numbering of the paper; the modes and keys of `geom_search.py` and `boundary_exact.py` named after draft lemma numbers were renamed; outputs go to `out/`; `summarize.py` takes the folder as an argument; a note on the false alarm was added to the docstring of `degenerate.py`.
