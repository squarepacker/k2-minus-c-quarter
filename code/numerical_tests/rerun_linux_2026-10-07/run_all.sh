#!/bin/bash
# Re-run of this repository on a clean Ubuntu 24.04 server (2 vCPU); the repository is expected as ~/repo.zip.
# Two lanes in parallel, one process each. Everything is written under ~/results; ~/results/DONE marks the end.
set -u
R=~/results; mkdir -p "$R"
exec > >(tee -a "$R/main.log") 2>&1
echo "start $(date -u +%FT%TZ)"
sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq python3-venv unzip time > /dev/null
python3 -m venv ~/venv
. ~/venv/bin/activate
pip install -q --upgrade pip
pip install -q mpmath numpy scipy
{ python --version; pip freeze; uname -a; nproc; free -m; } > "$R/env.txt" 2>&1
mkdir -p ~/A ~/B
unzip -q ~/repo.zip -d ~/A
unzip -q ~/repo.zip -d ~/B
( cd ~/A && sha256sum -c SHA256SUMS > "$R/sha256_check.txt" 2>&1; echo "sha256sum exit $?" >> "$R/sha256_check.txt" )
for L in A B; do mkdir -p ~/$L/code/numerical_tests/P3_P5/results/original/runs; unzip -q ~/$L/code/numerical_tests/P3_P5/results/original/runs.zip -d ~/$L/code/numerical_tests/P3_P5/results/original/runs; done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0

run() {  # run LANE DIR LABEL CMD...
  local lane=$1 dir=$2 label=$3; shift 3
  echo "[$(date -u +%T)] start $label: $*" >> "$R/${lane}_index.txt"
  ( cd "$dir" && /usr/bin/time -f "ELAPSED %e s  MAXRSS %M KB  EXIT %x" -o "$R/${lane}_${label}.time" "$@" > "$R/${lane}_${label}.log" 2>&1 )
  echo "[$(date -u +%T)] done  $label: $(tail -1 "$R/${lane}_${label}.time")" >> "$R/${lane}_index.txt"
}

laneA() {
  local C=~/A/code T=~/A/code/numerical_tests
  run A $C constants python k14_constants_check.py
  run A $T/P1_P2 p12_thresholds python thresholds.py
  run A $T/P1_P2 p12_reflection python reflection.py
  run A $T/P1_P2 p12_merge python merge_points.py
  run A $T/P1_P2 p12_witnesses python witnesses.py
  run A $T/P1_P2 p12_aggregate python aggregate.py results/original
  run A $T/P1_P2 p12_smoke python run_families.py 60 smoke quick
  run A $T/P1_P2 p12_tvp python run_families.py 120 quick_tvp quick tower,valley,pile
  run A $T/P1_P2 p12_fuzz python run_families.py 60 fuzz quick fuzz
  run A $T/P1_P2 p12_conseq python consequence.py 90 conseq_quick quick
  run A $T/P1_P2 p12_side_full python side_contacts.py
  run A $T/P3_P5 p35_summarize python summarize.py results/original/runs
  run A $T/P3_P5 p35_degenerate python degenerate.py
  run A $T/P3_P5 p35_deaths python deaths.py
  run A $T/P3_P5 p35_mutations python mutations.py
  run A $T/P3_P5 p35_cross python cross_check.py
  run A $T/P3_P5 p35_boundary python boundary_exact.py
  run A $T/P3_P5 p35_starmin03 python geom_search.py starmin 50 0.3
  run A $T/P3_P5 p35_star03 python geom_search.py star 50 0.3
  run A $T/P3_P5 p35_starmin09 python geom_search.py starmin 50 0.9
  for f in valley smoke bigvalley mechA three fan strip; do run A $T/P3_P5 p35_fam_$f python run_family.py $f 60; done
  run A $T/P3_P5 p35_merge python run_family.py merge 50
  run A $T/P3_P5 p35_near45 python run_family.py near45 50
  run A $T/P3_P5 p35_mrand python run_family.py mrand 60
  run A $T/P3_P5 p35_mrand9000 python run_family.py mrand 20 9000
  run A $T/P3_P5 p35_rowjam python run_family.py rowjam 50
  run A $T/P3_P5 p35_regime python run_family.py regime 40
  run A $T/P3_P5 p35_tri python run_family.py tri 50
  run A $T/P4 p4_unit python run_tests.py unit
  run A $T/P4 p4_summ_all python summarize.py results/original/res_scaled.jsonl results/original/res_merge.jsonl results/original/res_defs.jsonl results/original/res_dag.jsonl results/original/res_manyvar.jsonl results/original/res_ties.jsonl
  run A $T/P4 p4_summ_r2 python summarize.py results/original/res_dag.jsonl results/original/res_manyvar.jsonl results/original/res_ties.jsonl
  run A $T/P4 p4_negctl python negative_control.py
  run A $T/P4 p4_ties python ties.py
  run A $T/P4 p4_wall python wall_sweep.py
  run A $T/P4 p4_dag python run_tests.py dag,near45,dag 6 4004 120 out/quick_dag.jsonl 60 manyvar
  run A $T/P4 p4_manyvar python run_tests.py valley,deaths,floor,rows,zigzag,floor_paper 6 5005 120 out/quick_manyvar.jsonl 60 manyvar
  run A $T/P4 p4_scaled python run_tests.py valley,zigzag,rows,near45,deaths,floor,jam,tower,wall,wallmax 10 1001 100 out/quick_scaled.jsonl 30
  run A $T/P4 p4_paper python run_tests.py valley_paper,zigzag_paper,rows_paper,deaths_paper,floor_paper,jam_paper,tower_paper 7 2002 100 out/quick_paper.jsonl 30
  run A $T/P4 p4_exceptional_full python exceptional_types.py 420
  run A $T/P6 p6_run python p6x_run.py 3 25 120 quick_scaled_3.jsonl scaled
  run A $T/P6 p6_summary python p6x_summary.py quick_scaled_3.jsonl
  run A $T/P6 p6_targeted python p6x_targeted.py
  run A $T/P6 p6_ties python p6x_ties.py
  run A $T/P6 p6_negctl python p6x_negctl.py 160 40
  run A $T/P7 p7_run python p7x_run.py blocked2 505 90 quick_D.jsonl
  run A $T/P7 p7_summary python p7x_summary.py quick_D.jsonl
  run A $T/P7 p7_t7 python p7x_t7.py
  run A $T/P7 p7_targeted2 python p7x_targeted2.py
  run A $T/P7 p7_targeted python p7x_targeted.py
  run A $T/P7 p7_example python p7x_example.py
  run A $T/P8 p8_main_quick python p8x_main.py --quick
  run A $T/P8 p8_scan_quick python p8x_scan.py --quick
  run A $T/P8 p8_asym_quick python p8x_asym.py --quick
  run A $T/P8 p8_asym_fix python p8x_asym_fix.py
  run A $T/P8 p8_sens python p8x_sens.py
  run A $T/assembly asm_run python asmx_run.py quickC 13 40 150 valley,hill,columns,floortouch,rows,lshape 20,24,28
  run A $T/assembly asm_summary python asmx_summary.py quickC.jsonl
  run A $T/assembly asm_T2 python asmx_targeted.py T2
  run A $T/assembly asm_t2b python asmx_t2b.py
  run A $T/assembly asm_replay python asmx_replay.py campC 13 valley,hill,columns,floortouch,rows,lshape 20,24,28 results/original/flagged.json
  run A $T/P3_P5 p35_starmin_rest python geom_search.py starmin 900 0.01,1.4
  run A $T/P3_P5 p35_star_rest python geom_search.py star 900 1e-6,1e-3,1.4
  echo "[$(date -u +%T)] lane A finished" >> "$R/A_index.txt"
}

laneB() {
  local T=~/B/code/numerical_tests
  run B $T/P8 p8_main_full python p8x_main.py
  run B $T/P8 p8_scan_full python p8x_scan.py
  run B $T/P8 p8_asym_full python p8x_asym.py
  run B $T/P3_P5 p35_geom_two python geom_search.py two 1500
  run B $T/P3_P5 p35_geom_three python geom_search.py three 2700
  echo "[$(date -u +%T)] lane B finished" >> "$R/B_index.txt"
}

laneA & PA=$!
laneB & PB=$!
wait $PA $PB
cd ~ && tar czf "$R/outputs.tgz" $(find A B -type d -name out) 2>/dev/null
echo "done $(date -u +%FT%TZ)" > "$R/DONE"
echo "done $(date -u +%FT%TZ)"
