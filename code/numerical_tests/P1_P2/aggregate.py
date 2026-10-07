# -*- coding: utf-8 -*-
"""Aggregate the run summaries (*_summary.json of run_families.py and consequence.py) into
one compact file.

Usage:  python aggregate.py [folder]      (default folder: out)
  writes out/aggregate_out.json.  With the folder results/original it recomputes the
  recorded aggregate results/original/aggregate_out.json from the recorded summaries
  (the per-family maxima are carried over for every key that starts with "worst_",
  whatever its name).
"""
import json, glob, os, sys
from common import out_path

SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

runs = {}
for path in sorted(glob.glob(os.path.join(SRC, "*_summary.json"))):
    fn = os.path.basename(path)
    try:
        with open(path, encoding="utf-8") as f:
            runs[fn] = json.load(f)
    except Exception as e:
        runs[fn] = {"error": str(e)}

agg = dict(trace_runs={}, conseq_runs={})
tot = dict(cfg=0, traces=0, contacts=0, bot=0, vbot=0, fails=0, anom=0, aux_viol=0, aux_false=0)
worst = {}
fam_all = {}
kinds = {}
terms = {}
flags = {}
for fn, s in runs.items():
    if "n_contacts" in s:
        agg["trace_runs"][fn] = dict(cfg=s["n_cfg"], traces=s["n_traces"], contacts=s["n_contacts"],
                                     bot=s["n_bot"], vbot=s["n_vbot"], fails=s["n_fail"], anom=s["n_anom"],
                                     elapsed=s.get("elapsed"))
        tot["cfg"] += s["n_cfg"]
        tot["traces"] += s["n_traces"]
        tot["contacts"] += s["n_contacts"]
        tot["bot"] += s["n_bot"]
        tot["vbot"] += s["n_vbot"]
        tot["fails"] += s["n_fail"]
        tot["anom"] += s["n_anom"]
        tot["aux_viol"] += s.get("n_aux_viol", 0)
        tot["aux_false"] += s.get("aux_ineq_false_contacts", 0)
        for k_, v in s.get("worst", {}).items():
            if k_ not in worst or v["val"] > worst[k_]["val"]:
                worst[k_] = dict(v, run=fn)
        for fam, st in s.get("families", {}).items():
            a = fam_all.setdefault(fam, dict(cfg=0, traces=0, contacts=0, bot=0))
            for kk in ("cfg", "traces", "contacts", "bot"):
                a[kk] += st[kk]
            for kk, vv in st.items():
                if kk.startswith("worst_"):
                    a[kk] = max(a.get(kk, 0), vv)
        for dct, src in ((kinds, "kinds"), (terms, "terms"), (flags, "flags")):
            for kk, vv in s.get(src, {}).items():
                dct[kk] = dct.get(kk, 0) + vv
        if s["n_fail"] or s["n_anom"]:
            agg.setdefault("fail_examples", []).append(dict(run=fn, fails=s["fails"][:5], anom=s["anomalies"][:5]))
        if s.get("aux_viol_examples"):
            agg.setdefault("aux_viol_examples", []).append(dict(run=fn, ex=s["aux_viol_examples"][:2]))
    elif "n_Zs" in s:
        agg["conseq_runs"][fn] = dict(cfg=s["n_cfg"], Zs=s["n_Zs"], P2f=s["n_P2f"], traces=s["n_traces"],
                                      contacts_withZ=s["n_contacts_withZ"], kinds_withZ=s["kinds_withZ"],
                                      viol=len(s["viol"]), back_q=s["back_q"],
                                      back_found=s["back_found_forb"], sanity=s["back_sanity"],
                                      min_slack_P2f=(s["min_slack_P2f"][0] if s["min_slack_P2f"] else None),
                                      min_margin_Zs=(s["min_margin_Zs"][0] if s["min_margin_Zs"] else None),
                                      anomalies=s["anomalies"][:5],
                                      c0s=sorted(set(r_.get("c0") for r_ in s["runs"])),
                                      ss=sorted(set(r_.get("s") for r_ in s["runs"])),
                                      skipped=[r_ for r_ in s["runs"] if "skipped" in r_][:6])
agg["totals"] = tot
agg["worst"] = worst
agg["families"] = fam_all
agg["kinds"] = kinds
agg["terms"] = terms
agg["flags"] = flags
with open(out_path("aggregate_out.json"), "w", encoding="utf-8") as f:
    json.dump(agg, f, indent=1, default=str)
print(json.dumps(dict(totals=tot, kinds=kinds, terms=terms, flags=flags,
                      worst={k: (v["val"], v["info"].get("fam"), v["info"].get("par"), v["info"].get("kind"),
                                 v["info"].get("m"), v["info"].get("aZ")) for k, v in worst.items()},
                      conseq=agg["conseq_runs"]), indent=1, default=str))
