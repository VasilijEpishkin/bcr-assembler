"""Калибровка правила схлопывания ошибок мыши на ридах человека, чья истинная молекула известна по UMI.

Риды человека несут UMI длиной 20 нт, поэтому у семьи UMI с >= MIN_UMI_READS ридов консенсус V..J — это
«истинная» молекула. Правило мыши (`collapse_errors`, дословно из ноутбука симуляции мыши) применяется
к тем же ридам, НЕ глядя на UMI, а его решения сверяются с консенсусом UMI:

  * absorbed_true_variants  - разные консенсусы UMI, поглощённые правилом в другую последовательность (ложные слияния)
  * error_reads_*           - риды, отличающиеся от консенсуса своей семьи (ошибки секвенирования/ПЦР): поглощены
                              в свой консенсус (верно), в другую последовательность (неверно) или оставлены
  * residual_false_variants - последовательности после схлопывания, не являющиеся консенсусом ни одной семьи UMI

Считается по прогонам (каждый прогон — одна библиотека, как образец мыши). Только читает данные; пишет один JSON-отчёт.
"""
import argparse
import ast
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

MIN_UMI_READS = 2


def barcode(sequence_id):
    m = re.search(r"(?:^|\|)BARCODE=([^|\s]+)", sequence_id or "")
    return m.group(1) if m else None


def consensus(seqs):
    """Поколоночный консенсус ридов одной длины; None при равенстве (как в ноутбуке симуляции)."""
    if len(set(seqs)) == 1:
        return seqs[0]
    out = []
    for col in zip(*seqs):
        top = Counter(col).most_common(2)
        if len(top) > 1 and top[0][1] == top[1][1]:
            return None
        out.append(top[0][0])
    return "".join(out)


def load_collapse(notebook, ratio):
    nb = json.loads(Path(notebook).read_text())
    source = next("".join(c["source"]) for c in nb["cells"] if "def collapse_errors(" in "".join(c["source"]))
    tree = ast.parse(source)
    funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in {"collapse_errors", "_one_mismatch"}]
    assert len(funcs) == 2, "collapse_errors/_one_mismatch not found"
    env = {"defaultdict": defaultdict, "ERROR_COLLAPSE_MAX_MISMATCHES": 1, "ERROR_COLLAPSE_MIN_RATIO": ratio}
    exec(compile(ast.Module(body=funcs, type_ignores=[]), str(notebook), "exec"), env)
    return env["collapse_errors"]


def read_run(airr, with_ids=False):
    """-> список (umi, locus, vj_sequence) для ридов с корректным UMI длиной 20 нт и без N; with_ids добавляет sequence_id в начало."""
    reads = []
    with open(airr, newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            full = (r.get("sequence") or "").upper()
            if "N" in full or not r.get("v_call") or not r.get("j_call"):
                continue
            bc = barcode(r.get("sequence_id"))
            if not bc or len(bc) != 20 or set(bc.upper()) - set("ACGT"):
                continue
            seq = full[int(r["v_sequence_start"]) - 1:int(r["j_sequence_end"])]
            reads.append((r["sequence_id"], bc, r["locus"], seq) if with_ids else (bc, r["locus"], seq))
    return reads


def evaluate(reads, collapse):
    fam = defaultdict(list)
    for umi, locus, seq in reads:
        fam[umi].append(seq)
    umi_cons = {}                                  # UMI -> консенсус ридов модальной длины (семьи с >= MIN_UMI_READS ридов)
    for umi, seqs in fam.items():
        length = Counter(map(len, seqs)).most_common(1)[0][0]
        group = [s for s in seqs if len(s) == length]
        if len(group) >= MIN_UMI_READS:
            c = consensus(group)
            if c is not None:
                umi_cons[umi] = c
    counts = Counter(seq for _, _, seq in reads)   # слепой подсчёт: правило видит риды, а не UMI
    parent_of = collapse(counts)
    cons_umis = Counter(umi_cons.values())         # консенсусная последовательность -> число независимых семей UMI
    res = {"reads": len(reads), "distinct_sequences": len(counts), "umi_families": len(fam),
           "verified_families": len(umi_cons), "distinct_true_sequences": len(cons_umis),
           "absorbed_variants": len(parent_of)}
    absorbed_true = [c for c in cons_umis if c in parent_of]
    res["absorbed_true_variants"] = len(absorbed_true)
    res["absorbed_true_variants_ge2_umis"] = sum(cons_umis[c] >= 2 for c in absorbed_true)
    res["absorbed_true_variants_umis_total"] = sum(cons_umis[c] for c in absorbed_true)
    res["true_variants_ge2_umis"] = sum(n >= 2 for n in cons_umis.values())
    # риды проверенных семей: ошибочен ли рид и что с ним сделало правило
    tally = Counter()
    for umi, locus, seq in reads:
        c = umi_cons.get(umi)
        if c is None:
            continue
        if seq == c:
            tally["reads_correct"] += 1
            tally["correct_reads_lost_by_absorption"] += int(seq in parent_of)
            continue
        tally["error_reads"] += 1
        target = parent_of.get(seq)
        if target is None:
            tally["error_reads_kept_as_variant"] += 1
        elif target == c:
            tally["error_reads_absorbed_into_own_consensus"] += 1
        else:
            tally["error_reads_absorbed_into_other_sequence"] += 1
    res.update(tally)
    retained = [s for s in counts if s not in parent_of]
    res["retained_sequences"] = len(retained)
    res["residual_false_variants"] = sum(s not in cons_umis for s in retained)
    res["retained_true_sequences"] = sum(s in cons_umis for s in retained)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--airr-dir", type=Path, required=True, help="results/PRJEB30386/post_annotation_filtered/airr_pass")
    ap.add_argument("--notebook", type=Path, required=True, help="ноутбук симуляции мыши с функцией collapse_errors")
    ap.add_argument("--runs", nargs="+", default=["ERR3004229", "ERR3004230", "ERR3004231", "ERR3004232"])
    ap.add_argument("--ratios", nargs="+", type=int, default=[3, 5, 10])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    assert not a.out.exists(), "choose a new report path; do not overwrite reports"
    report = {"airr_dir": str(a.airr_dir), "notebook": str(a.notebook), "min_umi_reads": MIN_UMI_READS, "runs": {}}
    for run in a.runs:
        reads = read_run(a.airr_dir / f"{run}.airr.tsv")
        report["runs"][run] = {}
        for ratio in a.ratios:
            r = evaluate(reads, load_collapse(a.notebook, ratio))
            report["runs"][run][str(ratio)] = r
            print(run, ratio, json.dumps(r), flush=True)
    a.out.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
