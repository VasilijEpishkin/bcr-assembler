#!/usr/bin/env python3
from pathlib import Path
import argparse
import csv
import gzip
import shutil


def open_text(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def read_fastq_records(path):
    with open_text(path, 'rt') as h:
        while True:
            head = h.readline()
            if not head:
                break
            seq = h.readline(); plus = h.readline(); qual = h.readline()
            if not qual:
                raise RuntimeError(f'Truncated FASTQ record in {path}')
            yield head, seq, plus, qual


def key_from_header(header: str) -> str:
    token = header.strip().split()[0]
    token = token.split('|')[0]
    token = token.removesuffix('/1').removesuffix('/2')
    return token.lstrip('@')


def sync_one(r1_path: Path, r2_path: Path, out1: Path, out2: Path, orphan1: Path, orphan2: Path):
    r2_map = {}
    for rec in read_fastq_records(r2_path):
        r2_map[key_from_header(rec[0])] = rec

    paired = orphan_r1 = 0
    seen = set()
    with open_text(out1, 'wt') as o1, open_text(out2, 'wt') as o2, open_text(orphan1, 'wt') as x1:
        for rec1 in read_fastq_records(r1_path):
            k = key_from_header(rec1[0])
            rec2 = r2_map.get(k)
            if rec2 is None:
                x1.writelines(rec1)
                orphan_r1 += 1
                continue
            o1.writelines(rec1)
            o2.writelines(rec2)
            paired += 1
            seen.add(k)

    orphan_r2 = 0
    with open_text(orphan2, 'wt') as x2:
        for k, rec2 in r2_map.items():
            if k not in seen:
                x2.writelines(rec2)
                orphan_r2 += 1

    return paired, orphan_r1, orphan_r2


def main(volume: str, dataset: str, force: bool = False):
    vol = Path(volume)
    src = vol / 'results' / dataset / 'pr_trimmed' / 'fastq'
    base = vol / 'results' / dataset / 'pr_trimmed_sync'
    out = base / 'fastq'
    orphans = base / 'orphan_fastq'
    logs = base / 'logs'

    if force and base.exists():
        shutil.rmtree(base)
    out.mkdir(parents=True, exist_ok=True)
    orphans.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    samples = sorted({p.name.replace('_1.pr.fastq.gz','').replace('_2.pr.fastq.gz','') for p in src.glob('*.pr.fastq.gz')})
    rows = []
    for sample in samples:
        r1 = src / f'{sample}_1.pr.fastq.gz'
        r2 = src / f'{sample}_2.pr.fastq.gz'
        if not (r1.exists() and r2.exists()):
            continue
        out1 = out / f'{sample}_1.pr.fastq.gz'
        out2 = out / f'{sample}_2.pr.fastq.gz'
        orphan1 = orphans / f'{sample}_1.orphan.fastq.gz'
        orphan2 = orphans / f'{sample}_2.orphan.fastq.gz'
        paired, orphan_r1, orphan_r2 = sync_one(r1, r2, out1, out2, orphan1, orphan2)
        rows.append({
            'sample': sample,
            'paired_records': paired,
            'orphan_r1_records': orphan_r1,
            'orphan_r2_records': orphan_r2,
            'r1_in': str(r1),
            'r2_in': str(r2),
            'r1_out': str(out1),
            'r2_out': str(out2),
        })
        print(f'[sync] {sample}: paired={paired:,} orphan_r1={orphan_r1:,} orphan_r2={orphan_r2:,}')

    summary = logs / 'pair_sync_summary.tsv'
    with open(summary, 'w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0].keys()) if rows else ['sample','paired_records','orphan_r1_records','orphan_r2_records','r1_in','r2_in','r1_out','r2_out'], delimiter='\t')
        w.writeheader()
        if rows:
            w.writerows(rows)
    print(summary)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Synchronize mouse primer-trimmed paired FASTQ by read id')
    ap.add_argument('volume')
    ap.add_argument('dataset')
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    main(args.volume, args.dataset, force=args.force)
