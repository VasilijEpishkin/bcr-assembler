#!/usr/bin/env python3
"""Собирает референсный набор V'DJer для вида/цепи без готового набора (например, IGH мыши).

Повторяет структуру vdjer_human_references/<цепь>/:
  ig_vdj.fa      аллели V/D/J (V'DJer читает только последовательности, заголовки справочные)
  v_anchors.txt  16-меры, по одному на аллель V: --v-rule cys (по умолчанию; КОНЧАЕТСЯ за 7 нт до консервативного кодона Cys) или fixed16 (seq[-32:-16])
  j_anchors.txt  16-меры, по одному на аллель J: --j-rule fixed16 (по умолчанию, официальное seq[16:32]) или trp (НАЧИНАЕТСЯ за 2 нт до TGG)
  v_index/j_index  все 16-меры на расстоянии Хэмминга <= --max-dist от якоря, "код<TAB>расстояние"
                   (код: 2 бита на основание, A=0 T=1 C=2 G=3, первое основание старшее; см. seq_to_kmer.c)
  v_region.fa    геномная область V в ориентации ГЕНА (обратный комплемент прямого отрезка для локусов на обратной цепи)

Положение якорей сверено с набором человека: официальные якоря стоят на ФИКСИРОВАННОМ смещении (V: seq[-32:-16] для всех
142 генов V человека; J: seq[16:32]); правило от Cys воспроизводит его для аллелей, оканчивающихся на TGT GCn AGA.
--max-dist по умолчанию 4 (в V'DJer --am 4); набор человека хранит и расстояние 5, которое при --am 4 не используется.
"""
import argparse, collections, gzip, hashlib, itertools, json, re, sys
from pathlib import Path

import numpy as np

CODE = {"A": 0, "T": 1, "C": 2, "G": 3}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta(path):
    name, chunks, out = None, [], {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as h:
        for line in h:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        out[name] = "".join(chunks).upper()
    return out


def v_anchor(seq):
    """16-мер, кончающийся за 7 нт до последнего кодона Cys (TGT/TGC) в рамке среди последних 39 нт; None, если его нет."""
    for i in range(len(seq) - 3, max(len(seq) - 40, 0), -1):
        if seq[i:i + 3] in ("TGT", "TGC") and (len(seq) - i) % 3 == 0:
            end = i - 7
            return seq[end - 16:end] if end - 16 >= 0 else None
    return None


def j_anchor(seq):
    """16-мер, начинающийся за 2 нт до консервативного кодона Trp (TGG мотива W-G-x-G: TGGG[GCT]); None, если его нет."""
    m = re.search("TGGG[GCT]", seq)
    if not m or m.start() < 2 or m.start() + 14 > len(seq):
        return None
    return seq[m.start() - 2:m.start() + 14]


def v_anchor_fixed16(seq):
    """Официальное правило V'DJer: 16-мер, кончающийся за 16 нт до 3'-конца аллеля V (все 142 гена V человека)."""
    return seq[-32:-16] if len(seq) >= 32 else None


def j_anchor_fixed16(seq):
    """Официальное правило V'DJer: 16-мер, начинающийся через 16 нт от 5'-конца аллеля J."""
    return seq[16:32] if len(seq) >= 32 else None


def to_code(kmer):
    v = 0
    for b in kmer:
        v = (v << 2) | CODE[b]
    return v


def neighbour_masks(max_dist, k=16):
    """XOR-маски, меняющие ровно d оснований (d <= max_dist); XOR с ненулевым 2-битным значением всегда меняет основание."""
    masks, dists = [np.zeros(1, dtype=np.uint64)], [np.zeros(1, dtype=np.uint8)]
    for d in range(1, max_dist + 1):
        cur = []
        for pos in itertools.combinations(range(k), d):
            shifts = [2 * (k - 1 - p) for p in pos]
            for subs in itertools.product((1, 2, 3), repeat=d):
                cur.append(sum(s << sh for s, sh in zip(subs, shifts)))
        masks.append(np.array(cur, dtype=np.uint64))
        dists.append(np.full(len(cur), d, dtype=np.uint8))
    return np.concatenate(masks), np.concatenate(dists)


def write_index(anchors, path, max_dist):
    masks, mdist = neighbour_masks(max_dist)
    codes = np.concatenate([np.uint64(to_code(a)) ^ masks for a in anchors])
    dists = np.tile(mdist, len(anchors))
    order = np.lexsort((dists, codes))                       # по коду, затем по расстоянию
    codes, dists = codes[order], dists[order]
    first = np.ones(len(codes), dtype=bool)
    first[1:] = codes[1:] != codes[:-1]                      # для каждого кода оставить минимальное расстояние
    codes, dists = codes[first], dists[first]
    with open(path, "w") as out:
        out.write("".join(f"{c}\t{d}\n" for c, d in zip(codes.tolist(), dists.tolist())))
    return len(codes)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--imgt", required=True, help="FASTA аллелей IMGT с именами вида IGHV1-1*01 (например, mouse_IMGT+C.fa из TRUST4)")
    ap.add_argument("--chain", required=True, help="IGH, IGK или IGL")
    ap.add_argument("--genome", required=True, help="FASTA(.gz) хромосомы с локусом")
    ap.add_argument("--contig", required=True, help="имя контига в BAM, например chr12")
    ap.add_argument("--v-region", required=True, help="НАЧАЛО-КОНЕЦ геномной области V (с 1, включительно)")
    ap.add_argument("--reverse-strand", action="store_true", help="локус на обратной цепи: v_region.fa пишется обратным комплементом")
    ap.add_argument("--v-rule", choices=["cys", "fixed16"], default="cys",
                    help="cys: 16-мер, кончающийся за 7 нт до консервативного Cys (совпадает с fixed16 для аллелей на TGT GCn AGA); "
                         "fixed16: official V'DJer rule seq[-32:-16]")
    ap.add_argument("--j-rule", choices=["fixed16", "trp"], default="fixed16",
                    help="fixed16: официальное правило V'DJer seq[16:32]; trp: 16-мер, начинающийся за 2 нт до консервативного TGG")
    ap.add_argument("--max-dist", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    seqs = {n: s.replace(".", "") for n, s in read_fasta(a.imgt).items()}   # выравнивания IMGT с пропусками содержат символы "."
    pref = a.chain.upper()
    clean = re.compile("^[ACGT]+$")      # V'DJer завершается на любом основании, кроме A/C/G/T
    skipped = sorted(n for n, s in seqs.items() if n.startswith(pref) and not clean.match(s))
    V = {n: s for n, s in seqs.items() if n.startswith(pref + "V") and clean.match(s)}
    D = {n: s for n, s in seqs.items() if n.startswith(pref + "D") and clean.match(s)}
    J = {n: s for n, s in seqs.items() if n.startswith(pref + "J") and clean.match(s)}
    v_fn = v_anchor if a.v_rule == "cys" else v_anchor_fixed16
    j_fn = j_anchor_fixed16 if a.j_rule == "fixed16" else j_anchor
    v_all = {n: v_fn(s) for n, s in V.items()}
    j_all = {n: j_fn(s) for n, s in J.items()}
    v_anchors = sorted({x for x in v_all.values() if x})
    j_anchors = sorted({x for x in j_all.values() if x})
    if not v_anchors or not j_anchors:
        sys.exit("no anchors found")

    with open(out / "ig_vdj.fa", "w") as f:
        for group in (V, D, J):
            for n, s in group.items():
                f.write(f">{n}|IMGT|{pref}\n{s}\n")
    (out / "v_anchors.txt").write_text("".join(x + "\n" for x in v_anchors))
    (out / "j_anchors.txt").write_text("".join(x + "\n" for x in j_anchors))
    nv = write_index(v_anchors, out / "v_index", a.max_dist)
    nj = write_index(j_anchors, out / "j_index", a.max_dist)

    start, end = map(int, a.v_region.split("-"))
    chrom = read_fasta(a.genome)
    contig_seq = next(iter(chrom.values())) if a.contig not in chrom else chrom[a.contig]
    region = contig_seq[start - 1:end]
    if a.reverse_strand:
        region = region.translate(COMP)[::-1]
    (out / "v_region.fa").write_text(f">{a.contig}:{start}-{end}\n{region}\n")

    fw = sum(x in region for x in v_anchors)
    prov = {
        "chain": pref, "imgt": str(a.imgt), "imgt_sha256": sha256(a.imgt), "genome": str(a.genome), "genome_sha256": sha256(a.genome),
        "v_rule": a.v_rule, "j_rule": a.j_rule, "contig": a.contig, "v_region": a.v_region, "reverse_strand": a.reverse_strand, "max_dist": a.max_dist,
        "alleles": {"V": len(V), "D": len(D), "J": len(J)}, "alleles_skipped_non_ACGT": skipped,
        "v_alleles_without_anchor": sorted(n for n, x in v_all.items() if not x),
        "j_alleles_without_anchor": sorted(n for n, x in j_all.items() if not x),
        "v_anchors": len(v_anchors), "j_anchors": len(j_anchors), "v_index_lines": nv, "j_index_lines": nj,
        "v_region_length": len(region), "v_region_non_ACGT": sum(1 for c in region if c not in "ACGT"), "v_anchors_found_in_v_region_gene_orientation": fw,
    }
    (out / "provenance.json").write_text(json.dumps(prov, indent=1) + "\n")
    print(json.dumps({k: v for k, v in prov.items() if not k.endswith("_sha256")}, indent=1))


if __name__ == "__main__":
    main()
