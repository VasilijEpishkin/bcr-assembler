"""Оценка репертуара IgReC по эталону из консенсусов UMI человека (тот же эталон, что в calibrate_collapse_on_human_umi.py).

IgReC запускается вслепую (без UMI) на post_annotation_filtered/fastq/<прогон>_filtered.fastq.gz. Эталон — консенсус V..J
семей UMI с >= 2 ридами. Кластер «несёт» истинную последовательность, если её V..J без TRIM нт с каждого конца целиком
входит в консенсус кластера: IgReC обрезает концы иначе, чем интервал V..J в AIRR (на ERR3004230 целиком находится
54% истинных последовательностей, с обрезкой 30 нт — 93%), поэтому сравнение целых последовательностей было бы нечестным.

Метрики совпадают с отчётом правила схлопывания, чтобы их можно было сравнить построчно.
"""
import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from calibrate_collapse_on_human_umi import MIN_UMI_READS, consensus, read_run

K = 25
TRIM = 30


def read_fasta(path):
    name, out = None, {}
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            out[name] = []
        elif name:
            out[name].append(line.upper())
    return {n: "".join(s) for n, s in out.items()}


def cluster_id(header):
    m = re.search(r"cluster___(\d+)___size___(\d+)", header)
    return (m.group(1), int(m.group(2))) if m else (header.split()[0], 0)


def evaluate(reads, repertoire_fa, rcm, trim=TRIM):
    core = (lambda t: t[trim:len(t) - trim]) if trim else (lambda t: t)
    fam, read_seq, read_umi = defaultdict(list), {}, {}
    for rid, umi, locus, seq in reads:
        fam[umi].append(seq)
        read_seq[rid], read_umi[rid] = seq, umi
    umi_cons = {}
    for umi, seqs in fam.items():
        length = Counter(map(len, seqs)).most_common(1)[0][0]
        group = [s for s in seqs if len(s) == length]
        if len(group) >= MIN_UMI_READS:
            c = consensus(group)
            if c is not None:
                umi_cons[umi] = c
    cons_umis = Counter(umi_cons.values())

    clusters = {}                                      # id кластера -> (консенсус, размер)
    for h, s in read_fasta(repertoire_fa).items():
        cid, size = cluster_id(h)
        clusters[cid] = (s, size)
    by_kmer = defaultdict(set)                         # 25-мер из середины каждой истинной последовательности
    for t in cons_umis:
        mid = len(t) // 2
        by_kmer[t[mid:mid + K]].add(t)                 # средний 25-мер лежит внутри обрезанной сердцевины
    carries = {}                                       # id кластера -> множество истинных последовательностей в нём
    for cid, (s, _) in clusters.items():
        found = set()
        for i in range(len(s) - K + 1):
            for t in by_kmer.get(s[i:i + K], ()):
                if core(t) in s:
                    found.add(t)
        carries[cid] = found
    recovered = set().union(*carries.values()) if carries else set()
    recovered5 = set().union(*(carries[cid] for cid, (_, sz) in clusters.items() if sz >= 5))  # «большой» репертуар IgReC

    read_cluster = {}
    for line in open(rcm):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2 and parts[1] != "":
            read_cluster[parts[0].split()[0]] = parts[1]

    res = {"reads": len(reads), "verified_families": len(umi_cons), "distinct_true_sequences": len(cons_umis),
           "true_variants_ge2_umis": sum(n >= 2 for n in cons_umis.values()),
           "clusters": len(clusters), "clusters_size_ge5": sum(sz >= 5 for _, sz in clusters.values()),
           "reads_in_rcm": len(read_cluster), "trim_nt": trim}
    lost = [t for t in cons_umis if t not in recovered]
    res["lost_true_variants"] = len(lost)
    res["lost_true_variants_ge2_umis"] = sum(cons_umis[t] >= 2 for t in lost)
    lost5 = [t for t in cons_umis if t not in recovered5]
    res["lost_true_variants_if_size_ge5"] = len(lost5)
    res["lost_true_variants_ge2_umis_if_size_ge5"] = sum(cons_umis[t] >= 2 for t in lost5)
    res["false_clusters"] = sum(1 for cid in clusters if not carries[cid])
    res["false_clusters_size_ge5"] = sum(1 for cid, (_, sz) in clusters.items() if sz >= 5 and not carries[cid])
    tally = Counter()
    for rid, seq in read_seq.items():
        c = umi_cons.get(read_umi[rid])
        if c is None:
            continue
        cid = read_cluster.get(rid)
        cons = clusters[cid][0] if cid in clusters else None
        kind = "correct_reads" if seq == c else "error_reads"
        tally[kind] += 1
        if cons is None:
            tally[kind + "_not_in_repertoire"] += 1
        elif core(c) in cons:
            tally[kind + "_in_own_consensus_cluster"] += 1
        elif core(seq) in cons:
            tally[kind + "_in_cluster_of_own_erroneous_sequence"] += 1
        else:
            tally[kind + "_in_other_cluster"] += 1
    res.update(tally)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--airr", type=Path, required=True)
    ap.add_argument("--igrec-out", type=Path, required=True, help="каталог выхода IgReC")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--trim", type=int, default=TRIM, help="сколько нт с каждого конца истинной V..J не сравнивается (0 — вся)")
    a = ap.parse_args()
    assert not a.out.exists(), "choose a new report path"
    reads = read_run(a.airr, with_ids=True)
    res = evaluate(reads, a.igrec_out / "final_repertoire.fa", a.igrec_out / "final_repertoire.rcm", a.trim)
    res["airr"], res["igrec_out"] = str(a.airr), str(a.igrec_out)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res))


if __name__ == "__main__":
    main()
