#!/usr/bin/env python3
from pathlib import Path
import argparse
import csv
import gzip

COMP = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')
MIN_OVERLAP = 20
MAX_OVERLAP = 180
MAX_ERROR_RATE = 0.12
SEED_LEN_OPTIONS = (16, 14, 12, 10)


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def open_text(path: Path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def iter_fastq(path: Path):
    with open_text(path, 'rt') as h:
        while True:
            head = h.readline()
            if head == '':
                break
            seq_line = h.readline()
            plus_line = h.readline()
            qual_line = h.readline()
            if seq_line == '' or plus_line == '' or qual_line == '':
                raise RuntimeError(f'Truncated FASTQ in {path}')
            seq = seq_line.rstrip('\n')
            plus = plus_line.rstrip('\n')
            qual = qual_line.rstrip('\n')
            yield head.rstrip('\n'), seq, plus, qual


def mismatch_count(a: str, b: str) -> int:
    mm = 0
    for x, y in zip(a, b):
        if x == 'N' or y == 'N':
            continue
        if x != y:
            mm += 1
    return mm


def choose_consensus_base(b1: str, q1: str, b2: str, q2: str):
    if b1 == b2 or b2 == 'N':
        return b1, q1
    if b1 == 'N':
        return b2, q2
    return (b1, q1) if ord(q1) >= ord(q2) else (b2, q2)


def find_overlap(s1: str, s2rc: str):
    # fast seed-based search near the tail of read1
    for seed_len in SEED_LEN_OPTIONS:
        if len(s2rc) < seed_len:
            continue
        seed = s2rc[:seed_len]
        start = max(0, len(s1) - MAX_OVERLAP)
        pos = s1.rfind(seed, start)
        while pos != -1:
            ov = len(s1) - pos
            if MIN_OVERLAP <= ov <= min(MAX_OVERLAP, len(s1), len(s2rc)):
                mm = mismatch_count(s1[-ov:], s2rc[:ov])
                if mm / ov <= MAX_ERROR_RATE:
                    return ov, mm
            pos = s1.rfind(seed, start, pos)
    # fallback bounded scan, prefer longest acceptable overlap
    best = None
    max_ov = min(MAX_OVERLAP, len(s1), len(s2rc))
    for ov in range(max_ov, MIN_OVERLAP - 1, -1):
        mm = mismatch_count(s1[-ov:], s2rc[:ov])
        er = mm / ov
        if er <= MAX_ERROR_RATE:
            if best is None or ov > best[0] or (ov == best[0] and mm < best[1]):
                best = (ov, mm)
    return best


def merge_pair(seq1: str, qual1: str, seq2: str, qual2: str):
    s2 = rc(seq2)
    q2 = qual2[::-1]
    hit = find_overlap(seq1, s2)
    if hit is None:
        return None
    ov, _mm = hit
    prefix = seq1[:-ov]
    suffix = s2[ov:]
    prefix_q = qual1[:-ov]
    suffix_q = q2[ov:]
    cons_seq = []
    cons_qual = []
    a = seq1[-ov:]
    aq = qual1[-ov:]
    b = s2[:ov]
    bq = q2[:ov]
    for x, qx, y, qy in zip(a, aq, b, bq):
        base, q = choose_consensus_base(x, qx, y, qy)
        cons_seq.append(base)
        cons_qual.append(q)
    return prefix + ''.join(cons_seq) + suffix, prefix_q + ''.join(cons_qual) + suffix_q, ov


def key_from_header(h: str) -> str:
    return h.split()[0].removesuffix('/1').removesuffix('/2').lstrip('@')


def process_sample(sample: str, r1: Path, r2: Path, outdir: Path):
    out = outdir / f'{sample}_assemble-pass.fastq.gz'
    fail1 = outdir / f'{sample}-1_assemble-fail.fastq.gz'
    fail2 = outdir / f'{sample}-2_assemble-fail.fastq.gz'
    n1 = n2 = n_pass = n_fail = 0
    with open_text(out, 'wt') as pass_h, open_text(fail1, 'wt') as f1, open_text(fail2, 'wt') as f2:
        for rec1, rec2 in zip(iter_fastq(r1), iter_fastq(r2)):
            h1, s1, p1, q1 = rec1
            h2, s2, p2, q2 = rec2
            n1 += 1; n2 += 1
            if key_from_header(h1) != key_from_header(h2):
                raise RuntimeError(f'Header mismatch in {sample}: {h1} vs {h2}')
            merged = merge_pair(s1, q1, s2, q2)
            if merged is None:
                n_fail += 2
                f1.write(f'{h1}\n{s1}\n{p1}\n{q1}\n')
                f2.write(f'{h2}\n{s2}\n{p2}\n{q2}\n')
                continue
            mseq, mqual, ov = merged
            pass_h.write(f'{h1}\n{mseq}\n+\n{mqual}\n')
            n_pass += 1
    return {
        'sample': sample,
        'input_r1_records': n1,
        'input_r2_records': n2,
        'assembled_pairs': n_pass,
        'failed_reads': n_fail,
        'total_seen_records': n_pass + n_fail,
        'merge_rate_from_outputs': f'{(n_pass/(n_pass+n_fail) if (n_pass+n_fail) else 0):.6f}',
        'pass_fastq': str(out),
        'fail1_fastq': str(fail1),
        'fail2_fastq': str(fail2),
    }


def main(volume: str, dataset: str):
    vol = Path(volume)
    in_dir = vol / 'results' / dataset / 'pr_trimmed' / 'fastq'
    out_base = vol / 'results' / dataset / 'merged'
    out_fastq = out_base / 'fastq'
    out_qc = out_base / 'qc'
    out_logs = out_base / 'logs'
    out_fastq.mkdir(parents=True, exist_ok=True)
    out_qc.mkdir(parents=True, exist_ok=True)
    out_logs.mkdir(parents=True, exist_ok=True)
    samples = sorted({p.name.replace('_1.pr.fastq.gz','').replace('_2.pr.fastq.gz','') for p in in_dir.glob('*.pr.fastq.gz')})
    rows = []
    for i, sample in enumerate(samples, 1):
        r1 = in_dir / f'{sample}_1.pr.fastq.gz'
        r2 = in_dir / f'{sample}_2.pr.fastq.gz'
        print(f'[{i}/{len(samples)}] {sample}', flush=True)
        row = process_sample(sample, r1, r2, out_fastq)
        row['status'] = 'done'
        row['elapsed_sec'] = ''
        row['stdout'] = ''
        row['stderr'] = ''
        row['verbose_log'] = ''
        rows.append(row)
        qc_path = out_qc / 'assembly_qc.tsv'
        with open(qc_path, 'w', newline='') as h:
            w = csv.DictWriter(h, fieldnames=list(rows[0].keys()), delimiter='\t')
            w.writeheader(); w.writerows(rows)
        print(f"    pass={row['assembled_pairs']:,} fail_records={row['failed_reads']:,}", flush=True)
    print(out_qc / 'assembly_qc.tsv')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('volume')
    ap.add_argument('dataset')
    args = ap.parse_args()
    main(args.volume, args.dataset)
