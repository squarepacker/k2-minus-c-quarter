# Re-run of 7 October 2026: comparison with the README claims

80 of 80 claims reproduced. Fields not compared: L7_tests, analysis_sec, cat_Xtilt_gt_half_alpha, cat_Y_near45, cat_Y_terminator, cat_floor_pair, cat_foreign_within_0.01delta, cat_multi_piece, elapsed, rss_mb, sec, t, t_flow, t_total, time.

| folder | claim | reproduced | detail |
|---|---|---|---|
| all | SHA256SUMS of the uploaded repository | yes | 176 of 176 files OK |
| all | code/k14_constants_check.py | yes | 538 checks, last line ALL OK |
| all | every command exits with code 0 | yes | 73 commands; non-zero: none |
| all | no Python traceback in any log | yes | none |
| P1_P2 | thresholds_out.json equals the record | yes | identical (timing fields dropped) |
| P1_P2 | reflection_out.json equals the record | yes | identical (timing fields dropped) |
| P1_P2 | merge_points_out.json equals the record | yes | identical (timing fields dropped) |
| P1_P2 | witnesses_out.json equals the record | yes | identical (timing fields dropped) |
| P1_P2 | side_contacts_out.json equals the record (full run, not re-run before) | yes | identical (timing fields dropped) |
| P1_P2 | aggregate_out.json equals the record | yes | byte-identical up to line endings |
| P1_P2 | p12_smoke: 0 failures, 0 anomalies | yes | 593 configurations, 27815 contacts (time budget; the README's laptop counts are smaller) |
| P1_P2 | p12_tvp: 0 failures, 0 anomalies | yes | 583 configurations, 36737 contacts (time budget; the README's laptop counts are smaller) |
| P1_P2 | p12_fuzz: 0 failures, 0 anomalies | yes | 1511 configurations, 27471 contacts (time budget; the README's laptop counts are smaller) |
| P1_P2 | consequence.py quick: 0 violations, 0 found by the backward search, sanity 87/87 | yes | done 6 62 292 0 14688 0 {'checked': 87, 'matched': 87} |
| P3_P5 | run_family.py valley: records equal the recorded ones | yes | 168 records on the server (168 recorded, compared 168): 168 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py smoke: records equal the recorded ones | yes | 1 records on the server (1 recorded, compared 1): 1 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py bigvalley: records equal the recorded ones | yes | 90 records on the server (90 recorded, compared 90): 90 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py mechA: records equal the recorded ones | yes | 48 records on the server (48 recorded, compared 48): 48 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py three: records equal the recorded ones | yes | 72 records on the server (72 recorded, compared 72): 72 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py fan: records equal the recorded ones | yes | 144 records on the server (144 recorded, compared 144): 144 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py strip: records equal the recorded ones | yes | 72 records on the server (72 recorded, compared 72): 72 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py merge: records equal the recorded ones | yes | 4226 records on the server (5040 recorded, compared 4226): 4201 identical, 25 identical except floating-point values (largest relative deviation 1.9e-16), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py near45: records equal the recorded ones | yes | 68 records on the server (68 recorded, compared 68): 40 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 28 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py mrand: records equal the recorded ones | yes | 1212 records on the server (4464 recorded, compared 1212): 1187 identical, 22 identical except floating-point values (largest relative deviation 3.0e-16), 3 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py rowjam: records equal the recorded ones | yes | 192 records on the server (2100 recorded, compared 192): 157 identical, 2 identical except floating-point values (largest relative deviation 1.2e-11), 33 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py regime: records equal the recorded ones | yes | 414 records on the server (3159 recorded, compared 414): 404 identical, 10 identical except floating-point values (largest relative deviation 3.0e-16), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py tri: records equal the recorded ones | yes | 318 records on the server (1082 recorded, compared 318): 33 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 285 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | run_family.py mrand 20 9000: records equal res_mrand_9000.jsonl | yes | 53 records on the server (47 recorded, compared 47): 47 identical, 0 identical except floating-point values (largest relative deviation 0.0e+00), 0 with other exact square coordinates (generated from floating-point values); violations on the server 0 |
| P3_P5 | summarize.py: out/summary.json equals the record | yes | byte-identical up to line endings |
| P3_P5 | degen_results.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | death_results.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | mutation_results.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | cross_check_results.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | boundary_exact.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | geom_starmin_0.3.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | geom_star_0.3.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | geom_starmin_0.9.json equals the record | yes | identical (timing fields dropped) |
| P3_P5 | geom_search.py starmin theta=0.01 (not re-run before) equals geom_starmin_0.01.json | yes | identical (timing fields dropped) |
| P3_P5 | geom_search.py starmin theta=1.4 (not re-run before) equals geom_starmin_1.4.json | yes | identical (timing fields dropped) |
| P3_P5 | geom_search.py star theta=1e-6 (not re-run before) equals geom_star_1e-6.json | yes | identical (timing fields dropped) |
| P3_P5 | geom_search.py star theta=1e-3 (not re-run before) equals geom_star_1e-3.json | yes | identical (timing fields dropped) |
| P3_P5 | geom_search.py star theta=1.4 (not re-run before) equals geom_star_1.4.json | yes | identical (timing fields dropped) |
| P3_P5 | geom_search.py two (not re-run before) equals geom_two_all.json | yes | 32 entries on the server, 32 recorded; the recorded ones: identical (timing fields dropped) |
| P3_P5 | geom_search.py three (not re-run before) equals geom_three_all.json | yes | 16 entries on the server, 12 recorded; the recorded ones: identical (timing fields dropped) |
| P4 | res_negctl.json equals the record | yes | identical (timing fields dropped) |
| P4 | res_ties.jsonl equals the record | yes | identical (timing fields dropped) |
| P4 | res_wallsweep.json equals the record | yes | identical (timing fields dropped) |
| P4 | quick_dag.jsonl equals the first 6 records of res_dag.jsonl (_defs read as _paper) | yes | identical (timing fields dropped) |
| P4 | quick_manyvar.jsonl equals the first 6 records of res_manyvar.jsonl (_defs read as _paper) | yes | identical (timing fields dropped) |
| P4 | quick_scaled.jsonl: 0 violations | yes | 10 configurations |
| P4 | quick_paper.jsonl: 0 violations | yes | 7 configurations |
| P4 | summarize.py output equals summary_all.txt | yes | identical (the record starts with a byte order mark) |
| P4 | summarize.py output equals summary_round2.txt | yes | identical (the record starts with a byte order mark) |
| P4 | exceptional_types.py 420 (full budget): the same pattern as the recorded run | yes | 4731 configurations (recorded 2938): type D 8209 heights, 0 failures; type 'side' 13262 heights, 4550 failures, least margin -4.0 |
| P4 | run_tests.py unit: 11 checks, all PASS | yes | 11 PASS |
| P6 | p6x_run.py 3 25: the records equal the first 25 of campaign 3 (except L7_tests and the added cat_* counters) | yes | identical (timing fields dropped) |
| P6 | targeted_out.json equals the record | yes | identical except 3 floating-point value(s), relative deviation at most 1.1e-16; first: /mergeT[55]/R1_leftmost/exc: 0.9862415794891394 vs 0.9862415794891393 |
| P6 | negctl_out.json equals results/rerun_2026-10-07/negctl_out.json | yes | identical (timing fields dropped) |
| P6 | p6x_ties.py console output equals ties_console.txt | yes | identical |
| P6 | p6x_summary.py: violations 0 | yes |  |
| P7 | p7x_run.py blocked2 505: the records equal the first records of campaign D (campD.jsonl, not included) | yes | 45 records; identical (timing fields dropped) |
| P7 | t7.jsonl equals the record | yes | identical (timing fields dropped) |
| P7 | targeted.jsonl equals the record | yes | identical (timing fields dropped) |
| P7 | targeted2.jsonl equals the record | yes | identical (timing fields dropped) |
| P7 | p7x_example.py output equals example_out.txt except the first line (label) | yes | 7 lines |
| P7 | p7x_summary.py: every viol_* counter is 0 | yes |  |
| assembly | asmx_run.py quickC: the 40 records equal the first 40 of campaign C (campC.jsonl, not included; run name quickC read as campC) | yes | identical except 216 floating-point value(s): relative deviation at most 1e-09 or absolute deviation at most 3.6e-15 (values near 0); first: [8]/variants[0]/floor/step3min: 3.3855053563812456 vs 3.385505356381245 |
| assembly | targeted_T2.json equals the record | yes | identical (timing fields dropped) |
| assembly | replay_campC.json equals the record | yes | identical (timing fields dropped) |
| assembly | targeted_T2b.json equals the record | yes | identical except 1 floating-point value(s), relative deviation at most 2.6e-12; first: /M_measure: 0.24381136863020814 vs 0.24381136862958197 |
| assembly | asmx_summary.py: every VIOL_* counter is 0 | yes |  |
| P8 | full mode: p8x_main_out.json equals the record (labels mapped) | yes | identical (timing fields dropped) |
| P8 | full mode: p8x_scan_out.json equals the record (labels mapped) | yes | identical (timing fields dropped) |
| P8 | full mode: p8x_asym_out.json equals the record (labels mapped) | yes | identical (timing fields dropped) |
| P8 | full mode: every number with 4 or more decimals in p8x_main_out.txt equals the record | yes | 952 numbers |
| P8 | full mode: every number with 4 or more decimals in p8x_scan_out.txt equals the record | yes | 197 numbers |
| P8 | full mode: every number with 4 or more decimals in p8x_asym_out.txt equals the record | yes | 98 numbers |
| P8 | p8x_main.py (full): 623 checks, 0 fails, 0 discrepancies | yes | 623 checks, 0 fails, 0 discrepancies, 11 unsafe roundings (draft values, see the README) |
| P8 | p8x_asym_fix_out.json equals the record | yes | identical (timing fields dropped) |
| P8 | p8x_sens_out.json equals the record | yes | identical (timing fields dropped) |
