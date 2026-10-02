"""Почему сборки TRUST4 близки к эталону, но не точны: разбор одного образца готового бенчмарка.

Для каждого проверенного шаблона, чей лучший контиг не совпадает с ним точно по всей длине:
  * где отличия (V-часть / CDR3 / J-часть, непокрытый 5'- или 3'-конец) и какого они типа
    (замена, вставка, делеция) — по minimap2 --cs эталонной V..J против контигов;
  * чей это контиг: лучший контиг может быть точной сборкой ДРУГОГО (родственного) шаблона,
    то есть сам шаблон не собран, а собран его клональный родственник («затенение»);
  * кандидаты в химеры: контиги без точного совпадения с шаблоном, у которых 5'- и 3'-части
    точно совпадают с разными шаблонами.
Только читает данные; пишет одну JSON-сводку.
"""
import argparse
import collections
import json
import re
import subprocess
from pathlib import Path

import pandas as pd


def fasta(path):
    name, out = None, {}
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v).upper() for k, v in out.items()}


def paf_cs(ref, query, threads):
    """Лучшее выравнивание каждого запроса со строкой cs; пресет minimap2 для коротких ридов, как в бенчмарке."""
    cmd = ["minimap2", "-c", "--cs", "-x", "sr", "-N", "20", "-p", "0.5", "--secondary=yes", "-t", str(threads), str(ref), str(query)]
    best = {}
    for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        f = line.split("\t")
        nm = next(int(x[5:]) for x in f[12:] if x.startswith("NM:i:"))
        cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
        h = dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], t=f[5], tlen=int(f[6]),
                 ts=int(f[7]), te=int(f[8]), match=int(f[9]), block=int(f[10]), nm=nm, cs=cs)
        key = (h["match"] - nm, h["qe"] - h["qs"])
        if f[0] not in best or key > best[f[0]][0]:
            best[f[0]] = (key, h)
    return {q: v[1] for q, v in best.items()}


def diffs(h):
    """Координаты в запросе (эталоне) и тип каждого отличия. cs идёт по референсу:
    '*xy' замена, '+seq' основания эталона, которых нет в контиге, '-seq' лишние основания контига."""
    minus = h["strand"] == "-"
    pos, out = (h["qlen"] - h["qe"]) if minus else h["qs"], []
    for op, val in re.findall(r"([:*+\-])([0-9]+|[a-z]+)", h["cs"]):
        if op == ":":
            pos += int(val)
        elif op == "*":
            out.append((pos, "sub")); pos += 1
        elif op == "+":
            out.append((pos, "del_in_contig")); pos += len(val)
        elif op == "-":
            out.append((pos, "ins_in_contig"))
    return [(h["qlen"] - 1 - p, t) for p, t in out] if minus else out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench-dir", type=Path, required=True, help="<ветка>/benchmark/trust4/<образец> (содержит work/truth_vj.fasta, work/contigs.fasta)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--end", type=int, default=5, help="отличия не дальше стольких нт от конца V..J считаются краевыми")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    work = a.bench_dir / "work"
    truth, contigs = fasta(work / "truth_vj.fasta"), fasta(work / "contigs.fasta")
    rows = pd.read_csv(a.bench_dir / "per_template.tsv", sep="\t", dtype=str)
    rows = rows[rows.verified == "True"]
    junction = dict(zip(rows.template_id, rows.junction.fillna("")))

    t2c = paf_cs(work / "contigs.fasta", work / "truth_vj.fasta", a.threads)   # truth -> contigs
    c2t = paf_cs(work / "truth_vj.fasta", work / "contigs.fasta", a.threads)   # contigs -> truth
    exact_contig_of = {c: h["t"] for c, h in c2t.items() if h["nm"] == 0 and h["te"] - h["ts"] == h["tlen"]}

    summary, where, ndiff = collections.Counter(), collections.Counter(), collections.Counter()
    for tid in rows.template_id:
        h = t2c.get(tid)
        if h is None:
            summary["missing"] += 1
            continue
        if h["qs"] == 0 and h["qe"] == h["qlen"] and h["nm"] == 0:
            summary["exact"] += 1
            continue
        summary["inexact"] += 1
        owner = exact_contig_of.get(h["t"])
        if owner and owner != tid:
            summary["inexact_best_contig_is_exact_copy_of_other_template"] += 1
            continue
        summary["inexact_own_contig"] += 1
        L = h["qlen"]
        if h["qs"] > 0:
            where["5prime_not_covered"] += 1
        if h["qe"] < L:
            where["3prime_not_covered"] += 1
        j = junction.get(tid, "")
        js = truth[tid].find(j) if j else -1
        d = diffs(h)
        ndiff[min(len(d), 5)] += 1
        for p, kind in d:
            if p < a.end or p >= L - a.end:
                region = "end"
            elif js >= 0 and js <= p < js + len(j):
                region = "CDR3"
            elif js >= 0 and p >= js + len(j):
                region = "J"
            else:
                region = "V"
            where[f"{region}:{kind}"] += 1

    # кандидаты в химеры: контиг без точного шаблона, чьи 5'- и 3'-части по 60% точно совпадают с разными шаблонами
    parts = a.out.with_suffix(".chimera_parts.fasta")
    with open(parts, "w") as out:
        for c, s in contigs.items():
            if c in exact_contig_of or len(s) < 200:
                continue
            k = int(len(s) * 0.6)
            out.write(f">{c}|5\n{s[:k]}\n>{c}|3\n{s[-k:]}\n")
    ph = paf_cs(work / "truth_vj.fasta", parts, a.threads)
    exact_part = {q: h["t"] for q, h in ph.items() if h["nm"] == 0 and h["qe"] - h["qs"] == h["qlen"]}
    chim = [c for c in contigs if exact_part.get(f"{c}|5") and exact_part.get(f"{c}|3")
            and exact_part[f"{c}|5"] != exact_part[f"{c}|3"]]
    res = {"bench_dir": str(a.bench_dir), "verified_templates": len(rows), "templates": dict(summary),
           "own_contig_difference_sites": dict(where), "own_contig_differences_per_template(5=5+)": dict(sorted(ndiff.items())),
           "contigs": len(contigs), "contigs_exact_copy_of_a_template": len(exact_contig_of),
           "contigs_checked_for_chimera": sum(1 for c, s in contigs.items() if c not in exact_contig_of and len(s) >= 200),
           "chimera_candidates_5p_3p_exact_to_different_templates": len(chim), "chimera_examples": chim[:10]}
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
