"""Оценка пилота V'DJer на IGL человека без изменения основных результатов.

Recall считается относительно ридов, реально попавших в пилот; используется
бенчмарк проекта с дополнительной проверкой полных контигов.
"""
import collections
import hashlib
import json
import os
from pathlib import Path


def main():
    root = Path(os.environ.get("BCR_VOLUME", "/data/user/epishkin"))
    nbroot = root / "one-q/bbc68913-4463-433d-b647-c3b8a76555ba"
    dataset = os.environ.get("BCR_DATASET", "PRJEB30386")
    sample = os.environ.get("BCR_SAMPLE", "PRJEB30386_all_chains")
    chain = os.environ.get("BCR_CHAIN", "IGL")
    branch = root / "results" / dataset / "simulated" / os.environ.get("BCR_BRANCH", "insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_umicons_in5ng_rb3x_ss250")
    pilot = Path(os.environ.get("BCR_PILOT_ROOT", str(root / "vdjer_diag")))
    output = Path(os.environ.get("BCR_PILOT_OUTPUT", str(pilot / "benchmark_r4_rf0_audit_20260929")))
    output.mkdir(exist_ok=True)
    report_path = output / "summary.json"
    if report_path.exists():
        print(report_path.read_text())
        return
    notebook = nbroot / dataset / ("benchmark_assemblers_human.ipynb" if dataset == "PRJEB30386" else "benchmark_assemblers_mouse.ipynb")
    cells = json.loads(notebook.read_text())["cells"]
    source = next("".join(c["source"]) for c in cells if c["cell_type"] == "code" and "def load_truth(" in "".join(c["source"]))
    namespace = {}
    exec(compile(source, str(notebook), "exec"), namespace)
    truth = namespace["load_truth"](branch, sample)
    represented = collections.Counter()
    n_pairs = 0
    with open(pilot / "sub_R1.fastq") as r1, open(pilot / "sub_R2.fastq") as r2:
        while True:
            a, b = r1.readline(), r2.readline()
            if not a or not b:
                assert not a and not b
                break
            aid, bid = a.split()[0].removesuffix("/1"), b.split()[0].removesuffix("/2")
            assert aid == bid and aid.startswith("@")
            tid = aid[1:].rsplit("_rc_", 1)[0]
            assert tid in truth, tid
            represented[tid] += 1
            for _ in range(3):
                assert r1.readline() and r2.readline()
            n_pairs += 1
    for tid, item in truth.items():
        item["class"] = "pilot_present" if tid in represented else "no_reads"
    contigs_path = Path(os.environ.get("BCR_PILOT_CONTIGS", str(pilot / "r4_rf0/IGL/vdj_contigs.fa")))
    contigs = dict(namespace["iter_fasta"](contigs_path))
    assert contigs, "No diagnostic contigs"
    rows, crow = namespace["benchmark"](truth, contigs, output / "work", 2)
    namespace["write_tsv"](output / "per_template.tsv", rows)
    namespace["write_tsv"](output / "per_contig.tsv", crow)
    igl = [r for r in rows if r["locus"] == chain]
    summary, precision = namespace["summarize"](igl, crow)
    best, _ = namespace["best_by_query"](output / "work/contigs_to_truth.paf")
    strict = collections.Counter()
    for cid in contigs:
        hit = best.get(cid)
        if hit is None:
            strict["unmapped"] += 1
            continue
        qcov = (hit["qe"] - hit["qs"]) / hit["qlen"]
        identity = hit["match"] / hit["block"]
        key = "query_95pct_identity_98pct" if qcov >= .95 and identity >= .98 else "fails_query_95pct_identity_98pct"
        strict[key] += 1
        strict["best_template_" + truth[hit["t"]]["locus"]] += 1
    rc = str.maketrans("ACGT", "TGCA")
    sequences = list(contigs.values())
    sequences += [s.translate(rc)[::-1] for s in sequences]
    junctions = {}
    for row in igl:
        if row["verified"] and row["junction"]:
            j = row["junction"]
            junctions[j] = junctions.get(j, False) or row["reconstructability"] == "pilot_present"
    cdr3 = {}
    for state, present in (("pilot_present", True), ("no_reads_control", False)):
        js = [j for j, p in junctions.items() if p == present]
        cdr3[state] = {"clonotypes": len(js), "found_exact_either_strand": sum(any(j in s for s in sequences) for j in js)}
    report = {
        "pilot": f"{dataset} {sample} {chain}; not a full-run result",
        "branch": str(branch), "pairs": n_pairs,
        "represented_templates_by_locus": dict(collections.Counter(truth[t]["locus"] for t in represented)),
        "represented_pairs_by_locus": dict(sum((collections.Counter({truth[t]["locus"]: n}) for t, n in represented.items()), collections.Counter())),
        "benchmark_source_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "contigs_sha256": hashlib.sha256(contigs_path.read_bytes()).hexdigest(),
        "contig_length_min": min(map(len, contigs.values())),
        "contig_length_max": max(map(len, contigs.values())),
        "chain_recall": summary, "legacy_precision": precision,
        "native_header_junctions_contained": sum(cid.rsplit("_", 1)[-1] in seq for cid, seq in contigs.items()),
        "full_contig_alignment_check": dict(strict), "cdr3_exact_either_strand": cdr3,
        "limitations": ["Pilot-present does not imply complete V-J coverage.",
                       "Legacy precision can call a contig correct from a local alignment; inspect full-contig check.",
                       "Short reads may not phase distant SHM variants; even an exact CDR3 is not full V-J recovery."]}
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
