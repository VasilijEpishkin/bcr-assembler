"""Сводит оценки всех полных прогонов сборщиков и пилотов в одну папку истории.

Выход: папка docs/evaluation_history корня проекта
  full_runs_summary.tsv    строки recall / CDR3 каждого полного бенчмарка (прогон, датасет, ветка, схема оценки + исходные колонки)
  full_runs_precision.tsv  precision контигов по образцам каждого полного бенчмарка
  full_runs_structure.tsv  структура контигов по поддержке ридами и пороги поддержки (только новая схема)
  pilots.tsv               пилоты в длинном формате: pilot, condition, metric, value
  assembler_comparison_report.md  разделы 1 (человек) и 2 (мышь) отчёта: ветки 200 нг и текущая,
                           эталоны trim0 / trim45_15 (мышь ещё trim20), знаменатели «проверенные» и «все» V…J
"""
import csv
import glob
import json
import os
from pathlib import Path

import pandas as pd

V = Path(os.environ.get("BCR_VOLUME", "/data/user/epishkin/bcr-assembler"))
R = V / "results"
PL = Path(os.environ.get("BCR_PILOTS", "/data/user/epishkin/one-q/29eed1d7-f79d-4362-be81-21e7f7bfb640/pilots"))   # пилоты и проверки (Jupyter-каталог третьей таски)
OUT = V / "docs/evaluation_history"
OUT.mkdir(parents=True, exist_ok=True)
P = "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_"
T = "_rb3x_ss200-350h25_pcr1c30_pcr2c5_seq1e-3_al1e-1"
HF, MF = P + "umicons_in200ng" + T, P + "ec1x3_in200ng" + T   # 200 нг
HC, MC = P + "umicons" + T, P + "ec1x3" + T                                                              # текущая ветка

# (прогон, датасет, ветка, схема оценки, каталог с assemblers_*.tsv): только два последних сценария
FULL = [(run, ds, b, f"unique_vj_{tr}", R / ds / "simulated" / b / f"benchmark_unique_{tr}")
        for run, hb, mb in (("200ng", HF, MF), ("current", HC, MC))
        for ds, b, trims in (("ERP003950", mb, ("trim0", "trim20", "trim45_15")), ("PRJEB30386", hb, ("trim0", "trim45_15")))
        for tr in trims]


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
        if (d / name).exists() and (d / name).stat().st_size > 1:   # пустые таблицы (выходов сборщика нет) пропускаются
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
print({k: sum(len(x) for x in v) for k, v in parts.items()}, "pilot rows:", len(rows))

# ---------------------------------------------------------------- отчёт run5: разделы «Человек» и «Мышь»
TOOLS = [("trust4", "TRUST4"), ("mixcr", "MiXCR"), ("gather", "GATHeR"), ("rnaspades", "rnaSPAdes"),
         ("spades_meta", "metaSPAdes"), ("spades_rnaviral", "rnaviralSPAdes"), ("trinity", "Trinity")]
T0 = ("trim0", "без обрезки")
T20 = ("trim20", "обрезка 20 нт слева: специфичная часть FR1-праймера (20 нт у 18 из 19 праймеров Greiff 2014, у fw16 — 22 нт)")
T4515 = ("trim45_15", "обрезка 45 нт слева и 15 нт справа")
RECALL = ["exact_full", "one_diff_full", "near_full", "partial", "missing"]
PRECISION = ["exact_full", "near_full", "correct_partial", "unverified_match", "incorrect", "unmapped"]
STRUCTURE = ["exact_verified", "exact_unverified", "exact_fragment", "mosaic", "unverified_variant", "novel_error",
             "chimera", "other", "unmapped"]
RECON = ["full_coverage_read_overlap", "full_coverage_mate_linked_only", "full_coverage_unlinked",
         "partial_read_coverage", "no_reads"]
# (номер раздела, вид, датасет, локусы, эталоны, [(ветка, каталог)])
ORGS = [
    ("1", "Человек", "PRJEB30386", ["IGH", "IGK", "IGL"], [T0, T4515], [("200 нг", HF), ("текущая ветка", HC)]),
    ("2", "Мышь", "ERP003950", ["IGH"], [T0, T20, T4515], [("200 нг", MF), ("текущая ветка", MC)]),
]
NAMES = ("summary", "summary_all_templates", "precision", "structure_by_support", "support_thresholds")


def n(x):
    return f"{int(x):,}".replace(",", " ")


def cell(frac, total):
    return f"{n(round(float(frac) * total))} ({float(frac) * 100:.1f}%)"


def read(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if path.exists() else None


def load(dataset, branch, trim):
    d = R / dataset / "simulated" / branch / f"benchmark_unique_{trim}"
    return {k: read(d / f"assemblers_{k}.tsv") for k in NAMES}


def pick(rows, tool, by, group):
    for r in rows or []:
        if r["tool"] == tool and r["group_by"] == by and r["group"] == group:
            return r


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "| " + " | ".join(["---"] * len(head)) + " |"]
    return "\n".join(out + ["| " + " | ".join(r) + " |" for r in rows]) + "\n"


def tools_in(data):
    have = {r["tool"] for d in data.values() for k in ("summary", "precision") for r in (d[k] or [])}
    return [(t, name) for t, name in TOOLS if t in have] or TOOLS[:2]


def first(rows):
    return next((r for r in rows if r), None)


def section_recall(data, tools, trims):
    parts = []
    for tk, desc in trims:
        d = data[tk]
        if not d["summary"]:
            parts.append(f"**{tk}** ({desc}): не посчитано.\n")
            continue
        ver = {t: pick(d["summary"], t, "reconstructability", "ALL") for t, _ in tools}
        al = {t: pick(d["summary_all_templates"], t, "reconstructability", "ALL_TEMPLATES") for t, _ in tools}
        rv, ra = first(ver.values()), first(al.values())
        head = ["класс"]
        for _, nm in tools:
            head += [f"{nm}, проверенные ({n(rv['templates'])})" if rv else f"{nm}, проверенные",
                     f"{nm}, все шаблоны ({n(ra['templates'])})" if ra else f"{nm}, все шаблоны"]
        rows = [[f"`{c}`"] + [cell(r[f"{c}_frac"], int(r["templates"])) if r else "—"
                             for t, _ in tools for r in (ver[t], al[t])] for c in RECALL]
        parts.append(f"**{tk}** ({desc})\n\n" + table(head, rows))
    return "\n".join(parts)


def section_slices(data, tools, loci, trims):
    parts = []
    for tk, desc in trims:
        d = data[tk]
        if not d["summary"]:
            continue
        rows = []
        for by, groups in (("copies", ["1", "2-4", "5+"]), ("locus", loci), ("reconstructability", RECON)):
            for g in groups:
                vs = {t: pick(d["summary"], t, by, g) for t, _ in tools}
                al = {t: pick(d["summary_all_templates"], t, by, g) for t, _ in tools}
                rv, ra = first(vs.values()), first(al.values())
                if not rv and not ra:
                    continue
                rows.append([by, g, n(rv["templates"]) if rv else "—"]
                            + [cell(vs[t]["exact_full_frac"], int(vs[t]["templates"])) if vs[t] else "—" for t, _ in tools]
                            + [n(ra["templates"]) if ra else "—"]
                            + [cell(al[t]["exact_full_frac"], int(al[t]["templates"])) if al[t] else "—" for t, _ in tools])
        head = ["разрез", "подгруппа", "V…J проверенные"] + [nm for _, nm in tools] \
            + ["V…J все"] + [f"{nm}, все" for _, nm in tools]
        parts.append(f"**{tk}** ({desc}), exact_full\n\n" + table(head, rows))
    return "\n".join(parts)


def section_cdr3(data, tools, loci, trims):
    head = ["клонотипы"] + [f"{tk}: {nm}" for tk, _ in trims for _, nm in tools]
    rows = []
    for label, grp in [("все локусы", "ALL")] + ([(l, l) for l in loci] if len(loci) > 1 else []):
        row = [label]
        for tk, _ in trims:
            for t, _ in tools:
                r = pick(data[tk]["summary"], t, "cdr3_clonotypes", f"{grp}:sequenced")
                row.append(f"{cell(r['cdr3_in_contig_frac'], int(r['templates']))} из {n(r['templates'])}" if r else "—")
        rows.append(row)
    row = ["контроль (не прочитаны)"]
    for tk, _ in trims:
        for t, _ in tools:
            r = pick(data[tk]["summary"], t, "cdr3_clonotypes", "ALL:no_reads_control")
            row.append(f"{round(float(r['cdr3_in_contig_frac']) * int(r['templates']))} / {r['templates']}" if r else "—")
    return table(head, rows + [row])


def agg_precision(rows, tool, cols):
    rs = [r for r in rows or [] if r["tool"] == tool]
    tot = sum(int(r["contigs"]) for r in rs)
    return (tot, {c: sum(float(r[c]) * int(r["contigs"]) for r in rs) / tot for c in cols}) if tot else None


def section_precision(data, tools, trims):
    parts = []
    for tk, desc in trims:
        d = data[tk]
        if not d["precision"]:
            continue
        a = {t: agg_precision(d["precision"], t, [f"{c}_frac" for c in PRECISION]) for t, _ in tools}
        s = {t: agg_precision(d["precision"], t, [f"structure_{c}_frac" for c in STRUCTURE]) for t, _ in tools}
        head = ["класс"] + [nm for _, nm in tools]
        rows = [["контигов"] + [n(a[t][0]) if a[t] else "—" for t, _ in tools]]
        rows += [[f"`{c}`"] + [cell(a[t][1][f"{c}_frac"], a[t][0]) if a[t] else "—" for t, _ in tools] for c in PRECISION]
        parts.append(f"**{tk}** ({desc}), точность контигов\n\n" + table(head, rows))
        rows = [[f"`{c}`"] + [cell(s[t][1][f"structure_{c}_frac"], s[t][0]) if s[t] else "—" for t, _ in tools]
                for c in STRUCTURE]
        parts.append(f"**{tk}** ({desc}), структура контигов\n\n" + table(head, rows))
    return "\n".join(parts)


def section_support(data, trims):
    parts = []
    for tk, desc in trims:
        sup, thr = data[tk]["structure_by_support"], data[tk]["support_thresholds"]
        for t, nm in TOOLS:
            s = [r for r in sup or [] if r["tool"] == t]
            if not s:
                continue
            kind = "число ридов" if s[0]["support_type"] == "reads" else "k-мерное покрытие"
            rows = [[r["support_bin"], n(r["contigs"])] + [cell(r[c], int(r["contigs"])) for c in
                    ("novel_error", "mosaic", "exact_fragment", "exact_verified")] for r in s]
            parts.append(f"**{nm}, {tk}**: структура контигов по поддержке ({kind}; знаменатель — контиги в строке)\n\n"
                         + table(["поддержка", "контигов", "novel_error", "mosaic", "exact_fragment", "exact_verified"], rows))
            rows = [[r["min_support"], n(r["contigs_kept"]), f"{float(r['verified_exact_recall']) * 100:.1f}%",
                     f"{float(r['novel_error']) * 100:.1f}%"] for r in thr or [] if r["tool"] == t]
            parts.append(f"**{nm}, {tk}**: отсечение по порогу поддержки\n\n"
                         + table(["порог ≥", "осталось контигов", "exact recall (проверенные)", "novel_error среди оставшихся"], rows))
    return "\n".join(parts) if parts else "Данные поддержки не готовы.\n"


def build():
    out = []
    for num, org, dataset, loci, trims, scenarios in ORGS:
        out.append(f"## {num}. {org}\n")
        out.append("Сравнение сборщиков внутри вида. Для каждой ветки несколько вариантов эталона ("
                   + ", ".join(t for t, _ in trims) + ") и для каждого два знаменателя: **проверенные** V…J (основной) "
                   "и **все** уникальные V…J, включая непроверенные. «—» — для этой пары нет результатов benchmark. "
                   "В таблицах только сборщики с посчитанным benchmark.\n")
        out.append("\n".join(f"- **{t}**: {d}" for t, d in trims) + "\n")
        for si, (label, branch) in enumerate(scenarios, 1):
            data = {tk: load(dataset, branch, tk) for tk, _ in trims}
            tools = tools_in(data)
            out.append(f"### {num}.{si}. {org}, {label}\n")
            out.append("#### Полнота: какая доля эталонных V…J собрана\n\n" + section_recall(data, tools, trims))
            out.append("#### Разрезы\n\n" + section_slices(data, tools, loci, trims))
            out.append("#### CDR3\n\nКлонотип — (локус, junction) проверенных V…J; найден, если junction целиком входит в контиг. "
                       "Контроль — клонотипы, у которых все шаблоны в классе no_reads.\n\n" + section_cdr3(data, tools, loci, trims))
            out.append("#### Точность и структура контигов\n\nЗнаменатель — все контиги сборщика"
                       + ("; у мыши агрегировано по шести образцам" if org == "Мышь" else "") + ".\n\n"
                       + section_precision(data, tools, trims))
            out.append("#### Поддержка и порог\n\n" + section_support(data, trims))
    return "\n".join(out)


REPORT = OUT / "assembler_comparison_report.md"
if REPORT.exists():     # заменяются разделы 1 и 2; вводный и справочные разделы остаются как есть
    lines = REPORT.read_text().split("\n")
    i1 = next(i for i, l in enumerate(lines) if l.startswith("## 1."))
    i3 = next(i for i, l in enumerate(lines) if l.startswith("## 3. Сборщики"))
    while lines[i3 - 1].strip() in ("", "---"):
        i3 -= 1
    REPORT.write_text("\n".join(lines[:i1] + build().split("\n") + ["", "---", ""] + lines[i3:]))
