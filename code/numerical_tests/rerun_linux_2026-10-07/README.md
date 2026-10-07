# Re-run on a clean Linux machine (7 October 2026)

All quick checks listed in the folder READMEs, and the long runs that had not been repeated before, were run again on
a freshly installed Linux server from a copy of this repository. The results were then compared with the recorded
results and with the claims of the READMEs. This is a re-run of the same programs by the same author, on a different
operating system; it is not an independent test.

## Environment

- Amazon EC2 `r7i.large` (2 vCPU, 16 GiB), Ubuntu 24.04, a new server with nothing else installed.
- Python 3.12.3 with `mpmath` 1.4.1, `numpy` 2.5.3 and `scipy` 1.18.1 from PyPI (installed on 7 October 2026).
- Two lanes ran in parallel, one process each, with `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and `MKL_NUM_THREADS`
  set to 1 and `PYTHONHASHSEED=0`. Lane A ran the quick checks, lane B the long runs. The run took 30 minutes
  (00:55–01:25 UTC).

## Files

| file | content |
|---|---|
| `run_all.sh` | the script that was run; it expects the repository as `~/repo.zip` |
| `console_logs.zip` | console output of every command (`<lane>_<label>.log`), run time and peak memory (`.time`, GNU `time`), the order of the runs (`A_index.txt`, `B_index.txt`), and the output of `sha256sum -c SHA256SUMS` |
| `outputs.zip` | the `out/` folders written by the two lanes (`outputs/A/…`, `outputs/B/…`); a pickle cache file of `P3_P5` is left out |
| `compare_rerun.py` | the comparison; run `python compare_rerun.py` (standard library only) |
| `compare_out_2026-10-07.md` | its output, one row per claim |

`compare_rerun.py` can also compare the first records of campaign D of `P7/` and of campaign C of `assembly/` with the
recorded campaigns. These records are not in this repository (44 MB and 27 MB; see those READMEs), so without them
the two rows are reported as "not checked". The table `compare_out_2026-10-07.md` was made with them.

## Results

- `sha256sum -c SHA256SUMS` on the copy: 176 of 176 files OK.
- `code/k14_constants_check.py`: 538 checks, `ALL OK`.
- 73 commands; every one exited with code 0, and no log contains a Python traceback.
- 80 of 80 claims of the folder READMEs were reproduced, under the comparison rules below.

### Comparison rules

These are the rules stated in the folder READMEs:

- run-time and memory fields are not compared;
- keys that were renamed after the original runs are mapped (each folder README lists them);
- in `P6/`, the field `L7_tests` and the added `cat_*` counters are not compared;
- files written on Windows are compared up to line endings and a byte-order mark.

### Differences that the Windows re-runs did not show

1. **Last digits of floating-point values.** Some floating-point diagnostics differ in their last digits from the
   records made on Windows. The relative deviation is at most 1.2·10⁻¹¹. For values near 0, such as minimal margins
   of order 10⁻¹⁵ in `assembly/`, the absolute deviation is at most 3.6·10⁻¹⁵. The rational (exact) quantities are
   equal.
2. **Configurations built from floating-point coordinates.** Some generators of `P3_P5/` compute the positions of
   squares in floating point and then convert them exactly to rational numbers. Where the C library's
   trigonometric functions differ in the last bit, the rational coordinates differ, so the configuration is a
   slightly different one. This happened for 28 of 68 `near45`, 285 of 318 `tri`, 33 of 192 `rowjam` and 3 of
   1,212 `mrand` configurations. All checks pass on them (0 violations). Every other record equals the recorded one.
3. **Time budgets.** A run that stops at a time limit produces a different number of configurations on a different
   machine. For instance, the P1_P2 smoke run made 593 configurations here and 216 on the author's laptop. Where a
   README claims equality with the records, the same number of records is compared.

### Long runs not repeated before

All of these equal the records.

- **P8, full modes of `p8x_main.py`, `p8x_scan.py` and `p8x_asym.py`.** They took 43 s, 46 s and 50 s here, and 192 s,
  665 s and 707 s in the original runs.
  - The JSON outputs are equal after mapping the renamed labels.
  - Every number with four or more decimals in the text outputs equals the record (952, 197 and 98 numbers).
  - `p8x_main.py`: 623 checks, 0 fails, 0 discrepancies.
- **P3_P5, `geom_search.py`.**
  - `two`: all 32 entries identical.
  - `three`: with a time limit of 2700 s it produced 16 entries; the 12 recorded entries are identical.
  - `starmin` for θ = 0.01 and 1.4, and `star` for θ = 10⁻⁶, 10⁻³ and 1.4: identical.
- **P1_P2, `side_contacts.py` with its default budget (600 s; it finished in 61 s).** All 288 cases are identical to
  the record.
- **P4, `exceptional_types.py 420`.** It ran 4,731 configurations (2,938 in the record) and shows the same pattern as
  the record:
  - 8,209 exceptional heights of constant-death type, with no failure;
  - 13,262 heights of the sides of axis-parallel squares, with 4,550 failures (least margin −4);
  - the same example.

### Not repeated

These runs were not repeated:

- the long runs of `P1_P2/` (`main`, `main2`, the full fuzzer and the consequence runs);
- the full campaigns of `P4/`, `P6/`, `P7/` and `assembly/`, beyond the first records compared in the quick checks;
- the run chunks of `P3_P5/` that are not listed in its quick checks.
