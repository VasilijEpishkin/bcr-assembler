"""Сводит оценки всех полных прогонов сборщиков и пилотов в одну папку истории.

Выход: папка evaluation_history
  full_runs_summary.tsv    строки recall / CDR3 каждого полного бенчмарка (прогон, датасет, ветка, схема оценки + исходные колонки)
  full_runs_precision.tsv  precision контигов по образцам каждого полного бенчмарка
  full_runs_structure.tsv  структура контигов по поддержке ридами и пороги поддержки (только новая схема)
  pilots.tsv               пилоты в длинном формате: pilot, condition, metric, value
"""
import glob
import json
import os
from pathlib import Path

import pandas as pd

V = Path("/data/user/epishkin")
R = V / "results"
PL = V / "one-q/29eed1d7-f79d-4362-be81-21e7f7bfb640/pilots"      # пилоты и проверки (Jupyter-каталог третьей таски)
OUT = R / "evaluation_history"
OUT.mkdir(exist_ok=True)
HB = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_umicons_in5ng_rb3x_ss250"
MB = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_ec1x3_in5ng_rb3x_ss250"

HN = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_umicons_in5ng_rb3x_ssnone"
H3 = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_umicons_in100ng_rb3x_ss300_pcr1c30"
HF = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_umicons_in200ng_rb3x_ss200-350h25_pcr1c30_pcr2c5_seq1e-3_al1e-1"
MF = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_ec1x3_in200ng_rb3x_ss200-350h25_pcr1c30_pcr2c5_seq1e-3_al1e-1"
MN = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_ec1x3_in5ng_rb3x_ssnone"
M3 = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_ec1x3_in100ng_rb3x_ss300_pcr1c30"

# (прогон, датасет, ветка, схема оценки, каталог с assemblers_*.tsv)
FULL = [
    ("run1_v2", "PRJEB30386", "random_cut_amp_umimode_in10x_rb3x_ss250 (deleted)", "templates_whole_vj", R / "PRJEB30386/benchmark_v2"),
    ("run1_v2", "ERP003950", "random_cut_amp_in10x_rb3x_ss250 (deleted)", "templates_whole_vj", R / "ERP003950/benchmark_v2"),
    ("run2_base", "ERP003950", MB, "templates_whole_vj", R / "ERP003950/simulated" / MB / "benchmark"),
    ("run2_base", "ERP003950", MB, "unique_vj_trim20", R / "ERP003950/simulated" / MB / "benchmark_unique_trim20"),
    ("run2_base", "PRJEB30386", HB, "unique_vj_trim0", R / "PRJEB30386/simulated" / HB / "benchmark_unique_trim0"),
    ("run3_ssnone", "ERP003950", MN, "unique_vj_trim20", R / "ERP003950/simulated" / MN / "benchmark_unique_trim20"),
    ("run4_ss300", "ERP003950", M3, "unique_vj_trim20", R / "ERP003950/simulated" / M3 / "benchmark_unique_trim20"),
    ("run3_ssnone", "PRJEB30386", HN, "unique_vj_trim0", R / "PRJEB30386/simulated" / HN / "benchmark_unique_trim0"),
    ("run4_ss300", "PRJEB30386", H3, "unique_vj_trim0", R / "PRJEB30386/simulated" / H3 / "benchmark_unique_trim0"),
    ("run5_final", "ERP003950", MF, "unique_vj_trim20", R / "ERP003950/simulated" / MF / "benchmark_unique_trim20"),
    ("run5_final", "PRJEB30386", HF, "unique_vj_trim0", R / "PRJEB30386/simulated" / HF / "benchmark_unique_trim0"),
]
# оценка по внутренней части V…J: без первых 45 и последних 15 нт эталона, т.е. exact full[45:-15]
for run, ds, b in (("run2_base", "ERP003950", MB), ("run3_ssnone", "ERP003950", MN),
                         ("run4_ss300", "ERP003950", M3), ("run2_base", "PRJEB30386", HB),
                         ("run3_ssnone", "PRJEB30386", HN), ("run4_ss300", "PRJEB30386", H3),
                   ("run5_final", "ERP003950", MF), ("run5_final", "PRJEB30386", HF)):
    FULL.append((run, ds, b, "unique_vj_trim45_15", R / ds / "simulated" / b / "benchmark_unique_trim45_15"))


def tagged(path, run, ds, branch, scoring, **extra):
    df = pd.read_csv(path, sep="\t")
    for k, v in reversed([("run", run), ("dataset", ds), ("branch", branch), ("scoring", scoring)]):
        df.insert(0, k, v)
    for k, v in extra.items():
        df[k] = v
    return df


parts = {"summary": [], "precision": [], "structure": []}
for run, ds, branch, scoring, d in FULL:
    if not (d / "assemblers_summary.tsv").exists():
        print("нет результатов:", d)
        continue
    parts["summary"].append(tagged(d / "assemblers_summary.tsv", run, ds, branch, scoring))
    if (d / "assemblers_precision.tsv").exists():
        parts["precision"].append(tagged(d / "assemblers_precision.tsv", run, ds, branch, scoring))
    for name in ("assemblers_structure_by_support.tsv", "assemblers_support_thresholds.tsv"):
        if (d / name).exists() and (d / name).stat().st_size > 1:   # пустые таблицы (сборщика ещё нет) пропускаются
            parts["structure"].append(tagged(d / name, run, ds, branch, scoring, table=name[11:-4]))
for k, dfs in parts.items():
    pd.concat(dfs, ignore_index=True).to_csv(OUT / f"full_runs_{k}.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- пилоты (длинный формат)
rows = []


def add(pilot, condition, metrics, source):
    for m, v in metrics.items():
        if isinstance(v, (int, float, str, bool)) and v is not None:
            rows.append({"pilot": pilot, "condition": condition, "metric": m, "value": v, "source": str(source)})


p = PL / "trust4_param_pilot/decision.json"                       # параметры запуска TRUST4 (мышь ERR346600, база)
if p.exists():
    d = json.load(open(p))
    for cond, m in d["metrics"].items():
        add("trust4_params_ERR346600", cond, m, p)
    add("trust4_params_ERR346600", "decision", {"chosen": d["chosen"] or "default", "rule": d["rule"]}, p)

p = PL / "rnaspades_diag/summary.tsv"                              # режимы rnaSPAdes
if p.exists():
    for r in pd.read_csv(p, sep="\t").to_dict("records"):
        add("rnaspades_modes", r.pop("run"), r, p)

for s in sorted(glob.glob(str(PL / "vdjer_diag*/**/summary.json"), recursive=True)):   # оценки V'DJer
    j = json.load(open(s))
    rec = next((x for x in j.get("chain_recall", j.get("igl_recall", [])) if x["reconstructability"] == "ALL"), {})
    add("vdjer_benchmark", os.path.relpath(s, V), {**{f"recall_{k}": v for k, v in rec.items() if k != "reconstructability"},
        "contigs": j.get("legacy_precision", {}).get("contigs"), "contig_length_max": j.get("contig_length_max"),
        "cdr3_found": j["cdr3_exact_either_strand"]["pilot_present"]["found_exact_either_strand"],
        "cdr3_present": j["cdr3_exact_either_strand"]["pilot_present"]["clonotypes"]}, s)
for s in sorted(glob.glob(str(PL / "vdjer_diag*/auto_decision*.json"))):                # автопилот V'DJer
    j = json.load(open(s))
    for m in j["variants"]:
        add("vdjer_autopilot", f"{j['species']}:{Path(s).stem}:{m['variant']}",
            {"contigs": m["contigs"], **{f"param_{k}": v for k, v in m["params"].items()}}, s)
    add("vdjer_autopilot", f"{j['species']}:{Path(s).stem}:decision", {"approved": j["approved"]}, s)
for s in sorted(glob.glob(str(PL / "vdjer_diag*/**/exit.json"), recursive=True)):        # все запуски V'DJer
    j = json.load(open(s))
    fa = Path(s).parent / "vdj_contigs.fa"
    n = sum(1 for l in open(fa) if l.startswith(">")) if fa.exists() else 0
    add("vdjer_runs", os.path.relpath(Path(s).parent, V), {"returncode": j.get("returncode"), "contigs": n,
        "args": str(j.get("args") or j.get("command"))}, s)

p = PL / "collapse_calibration/report_human_umi_20260929.json"     # правило схлопывания против эталона по UMI человека
if p.exists():
    for run, ratios in json.load(open(p))["runs"].items():
        for ratio, m in ratios.items():
            add("collapse_rule_on_human_umi", f"{run}:ratio{ratio}", m, p)
for s in sorted(glob.glob(str(PL / "igrec_diag/fair/*.json"))):     # IgReC против правила, один оценщик
    add("igrec_vs_rule_trim30", Path(s).stem, json.load(open(s)), s)

for s in sorted(glob.glob(str(PL / "trust4_inexact/ERR*.json"))):   # разбор неточных сборок TRUST4 (вся V..J)
    j = json.load(open(s))
    add("trust4_inexact_whole_vj", Path(s).stem, {**{f"templates_{k}": v for k, v in j["templates"].items()},
        **{f"sites_{k}": v for k, v in j["own_contig_difference_sites"].items()},
        "chimera_candidates": j["chimera_candidates_5p_3p_exact_to_different_templates"], "contigs": j["contigs"]}, s)

for name, path in (("spades_modes_ERR346600", PL / "spades_modes_pilot/summary.tsv"),          # режимы SPAdes (мышь, база)
                   ):
    if path.exists():
        for r in pd.read_csv(path, sep="\t").to_dict("records"):
            add(name, r.pop("variant"), r, path)
p = PL / "bench_audit/out/audit_tables.tsv"                         # независимая проверка сводных таблиц бенчмарка
if p.exists():
    for r in pd.read_csv(p, sep="\t").to_dict("records"):
        add("benchmark_audit_tables", f"{r['branch']}:{r['dataset']}:{r['tool']}:{r['check']}",
            {"recomputed": str(r["recomputed"]), "reported": str(r["reported"])}, p)
p = PL / "bench_audit/out/audit.tsv"                                # независимая проверка построчных классов бенчмарка
if p.exists():
    for r in pd.read_csv(p, sep="\t").to_dict("records"):
        cond = f"{r.pop('branch')}:{r.pop('dataset')}:{r.pop('tool')}:{r.pop('sample')}"
        add("benchmark_audit", cond, r, p)

# сравнение запусков: строка = датасет × ветка × сборщик, столбцы = схемы оценки (доля эталонных V…J, собранных точно / с CDR3)
allrows = pd.concat(parts["summary"], ignore_index=True)
allrows = allrows[(allrows.group_by == "reconstructability") & (allrows.group == "ALL") & allrows.scoring.str.startswith("unique_vj")]
cmp_ = allrows.pivot_table(index=["dataset", "branch", "tool"], columns="scoring", values=["exact_full_frac", "cdr3_in_contig_frac"])
cmp_.columns = [f"{m}__{sc}" for m, sc in cmp_.columns]
cmp_.reset_index().to_csv(OUT / "full_runs_comparison.tsv", sep="\t", index=False)

rows = [r for r in rows if "tails" not in r["condition"] + r["source"]]   # в историю входят только основные ветки
pd.DataFrame(rows).to_csv(OUT / "pilots.tsv", sep="\t", index=False)
print({k: sum(len(x) for x in v) for k, v in parts.items()}, "pilot rows:", len(rows),
      collections_count := pd.DataFrame(rows).groupby("pilot").size().to_dict())
