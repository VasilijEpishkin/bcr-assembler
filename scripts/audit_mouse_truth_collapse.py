"""Проверка чувствительности правила схлопывания эталона мыши (только чтение).

Использует ту же реализацию, что ноутбук симуляции. Показывает, что удаляется,
но не может определить, ошибки ли это секвенирования или SHM.
"""
import argparse
import ast
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--airr", type=Path, required=True)
    parser.add_argument("--notebook", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "Choose a new report path; do not overwrite audits"
    counts, junctions = Counter(), {}
    total = excluded = 0
    opener = gzip.open if str(args.airr).endswith(".gz") else open
    with opener(args.airr, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            total += 1
            full = row["sequence"].upper()
            if "N" in full:
                excluded += 1
                continue
            start, end = int(row["v_sequence_start"]) - 1, int(row["j_sequence_end"])
            seq = full[start:end]
            counts[seq] += 1
            if row.get("cdr3_start") and row.get("junction"):
                jo, length = int(row["cdr3_start"]) - 4 - start, len(row["junction"])
                if 0 <= jo and jo + length <= len(seq):
                    junctions.setdefault(seq, Counter())[(jo, length)] += 1
    print(f"Loaded {total} rows; {len(counts)} exact variants", flush=True)
    nb = json.loads(args.notebook.read_text())
    source = next("".join(c["source"]) for c in nb["cells"] if "def collapse_errors(" in "".join(c["source"]))
    tree = ast.parse(source)
    funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in {"collapse_errors", "_one_mismatch"}]
    assert len(funcs) == 2
    env = {"defaultdict": defaultdict, "ERROR_COLLAPSE_MAX_MISMATCHES": 1}
    exec(compile(ast.Module(body=funcs, type_ignores=[]), str(args.notebook), "exec"), env)
    junction = {}
    for seq, spans in junctions.items():
        jo, length = spans.most_common(1)[0][0]
        junction[seq] = seq[jo:jo + length]
    stats = []
    for ratio in (3, 5, 10):
        env["ERROR_COLLAPSE_MIN_RATIO"] = ratio
        parent_of = env["collapse_errors"](counts)
        absorbed_reads = sum(counts[c] for c in parent_of)
        changed = [c for c, p in parent_of.items() if c in junction and p in junction and junction[c] != junction[p]]
        original_js = set(junction.values())
        retained_js = {j for s, j in junction.items() if s not in parent_of}
        examples = sorted(parent_of, key=lambda c: counts[c], reverse=True)[:5]
        result = {
            "ratio": ratio, "absorbed_variants": len(parent_of), "absorbed_reads": absorbed_reads,
            "retained_variants": len(counts) - len(parent_of),
            "absorbed_variants_read_support": {str(n): sum(counts[c] >= n for c in parent_of) for n in (2, 3, 5, 10, 20)},
            "absorbed_with_changed_junction": len(changed),
            "unique_junctions_before": len(original_js), "unique_junctions_after": len(retained_js),
            "unique_junctions_lost": len(original_js - retained_js),
            "highest_support_absorbed_examples": [
                {"child_reads": counts[c], "parent_reads": counts[parent_of[c]],
                 "child_sha256": hashlib.sha256(c.encode()).hexdigest(),
                 "parent_sha256": hashlib.sha256(parent_of[c].encode()).hexdigest(),
                 "changed_position_vj_0based": next(i for i, (a, b) in enumerate(zip(c, parent_of[c])) if a != b),
                 "junction_changed": junction.get(c) != junction.get(parent_of[c])} for c in examples]}
        stats.append(result)
        print(json.dumps(result), flush=True)
    report = {"source": str(args.airr), "source_rows": total, "excluded_N": excluded,
              "eligible_reads": sum(counts.values()), "exact_variants": len(counts),
              "exact_variants_at_least_two_reads": sum(n >= 2 for n in counts.values()),
              "collapse_code_sha256": hashlib.sha256(ast.unparse(ast.Module(body=funcs, type_ignores=[])).encode()).hexdigest(),
              "sensitivity": stats,
              "interpretation": "Counts describe algorithmic removals, not proven errors. Reads are not independent molecules without UMI. No base-quality or replicate test was performed."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
