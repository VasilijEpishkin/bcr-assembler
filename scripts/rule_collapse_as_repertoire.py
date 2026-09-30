"""Write the mouse collapse rule's output in IgReC repertoire format, so both are scored by evaluate_igrec_on_human_umi.py.

Each retained sequence becomes a cluster (size = its reads + reads of variants absorbed into it); every read maps to
the cluster of its own sequence or of the parent that absorbed it. Human reads, UMI ignored (blind, like IgReC).
"""
import argparse
from collections import Counter
from pathlib import Path

from calibrate_collapse_on_human_umi import load_collapse, read_run


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--airr", type=Path, required=True)
    ap.add_argument("--notebook", type=Path, required=True, help="mouse simulation notebook with collapse_errors")
    ap.add_argument("--ratio", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True, help="directory for final_repertoire.fa / .rcm")
    a = ap.parse_args()
    reads = read_run(a.airr, with_ids=True)
    counts = Counter(seq for _, _, _, seq in reads)
    parent_of = load_collapse(a.notebook, a.ratio)(counts)
    size = Counter()
    for seq, n in counts.items():
        size[parent_of.get(seq, seq)] += n
    cid = {seq: i for i, seq in enumerate(size)}
    a.out.mkdir(parents=True, exist_ok=False)
    with open(a.out / "final_repertoire.fa", "w") as fa:
        for seq, i in cid.items():
            fa.write(f">cluster___{i}___size___{size[seq]}\n{seq}\n")
    with open(a.out / "final_repertoire.rcm", "w") as rcm:
        for rid, _, _, seq in reads:
            rcm.write(f"{rid}\t{cid[parent_of.get(seq, seq)]}\n")
    print(a.out, "clusters", len(cid), "reads", len(reads))


if __name__ == "__main__":
    main()
