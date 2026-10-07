# Numerical tests (not part of the proof)

These folders contain the programs and the recorded results of numerical tests of the statements proved in the paper. The tests were run on scaled analogues of the construction: small packings (k of order 10 to 100) with parameters much larger than in the theorem (δ, the tilt and line thresholds), so that the effects become visible; where the constants of cited lemmas enter, some tests replace them by quantities measured on the same packing. They test the logical chain of the arguments, not the actual parameter regime (k ≥ 4.62·10¹², inclinations of order 10⁻⁶). **A test that finds no counterexample is a search, not a proof.**

Each test was written by an AI system (Claude, Anthropic) in a session separate from the one that wrote the corresponding proof draft, from the definitions and statements, with its own code; the tests were run on 6 October 2026 against the research drafts from which the paper was written. The programs were prepared for publication on 7 October 2026 (paths, quick modes, comments referring to the paper) without changing the logic of any test; each folder's README lists the changes and says which recorded outputs come from earlier versions of its scripts and cannot be reproduced exactly.

| Folder | Paper | What is tested | README |
|---|---|---|---|
| `P1_P2/` | Section 4 | Height identity and its corollary; quantization at bottom sides (five forms); the consequence used in Section 9 (no path of a scale flow reaches a forbidden square); reflection | `P1_P2/README.md` |
| `P3_P5/` | Sections 5 and 7 | At most two paths through a point; the two-source lemma and its hypothesis δ ≤ ½ cos 2α_max (sharp); double points and their rectangles; deaths, merges, gap overlap; end collisions; mutation tests | `P3_P5/README.md` |
| `P4/` | Section 6 | The shadow inequality in all parts, exceptional heights, the wall term; negative control without rule R1 | `P4/README.md` |
| `P6/` | Section 8 | Density and injectivity, the hull slit, the vertical comparison (Lemma 8.7), disjoint slits; negative controls | `P6/README.md` |
| `P7/` | Section 9 | The per-scale count and the lower bound for \|Z(s)\|; the set of tilt-terminated points; negative controls | `P7/README.md` |
| `P8/` | Sections 10–12 | The constants and thresholds, the integration inequality carried to the explicit bound, interval covers of [k0, 10¹⁰⁰], the asymptotic constant, sensitivity | `P8/README.md` |
| `assembly/` | Sections 3–11 | End-to-end check of the assembled chain on a scaled analogue, every definition implemented as stated; negative controls | `assembly/README.md` |

## Layout of each folder

- the programs (Python 3.12; `mpmath`, and in some folders `numpy` or `scipy`);
- `README.md`: what is tested and in which form, how to run it (quick checks of at most a few minutes, with measured run times and memory), the recorded results, the changes made for publication, and the limitations;
- `results/original/`: the outputs of the original runs, verbatim (in one record of `P3_P5/` an absolute path was replaced by `<path>`); large sets of records are zipped. Some records use internal names of the research drafts; each README explains them;
- `P6/results/rerun_2026-10-07/`: outputs of two re-runs made when the programs were prepared for publication (see `P6/README.md`);
- `rerun_linux_2026-10-07/`: a re-run of all quick checks, and of the long runs not repeated before, on a clean Linux machine, with the comparison against the records (80 of 80 claims of the READMEs reproduced; see its README);
- `out/`: created by the programs when they are run.

## Running the quick checks

```
pip install mpmath numpy scipy
cd code/numerical_tests/P4
python run_tests.py unit
```

Each README lists its quick checks with the expected output and the measured run time and memory (each at most a few minutes and under 200 MB); run them from the folder of the test. On Windows, if `python` opens the Microsoft Store instead of running, call the installed interpreter directly. The full runs of 6 October 2026 took from minutes to several hours; their commands are given in the READMEs, and some of them were made with earlier versions of the programs and cannot be reproduced exactly (each README says which).

## Summary of the recorded results

No test found a violation of a statement in the form used by the paper. Several tests were not vacuous: their negative controls (for example the paper's rules with rule R1 removed, or a bound halved) produced violations, and some inequalities were nearly attained in targeted cases. Two tests checked a weaker form than the one stated in the paper (Section 4: the loss bound α²q_y/(2 cos α) instead of (sec α − 1)q_y; Section 8: the end term 2.1·δβ̄ instead of 1.0000016·δβ̄); see the READMEs. The reviews of the drafts are summarized in `reviews/REVIEWS.md`.
