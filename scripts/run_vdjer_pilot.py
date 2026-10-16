#!/usr/bin/env python3
"""Воспроизводимые пилоты V'DJer на подвыборке симулированных ридов мыши и человека: подвыборка → STAR → V'DJer (rf=1, rf=0).

Готовые BAM и контиги не пересчитываются; BAM готовит `make_vdjer_bam.py`, описание — `docs/vdjer_bam.md`.
"""

import argparse
import gzip
import json
import os
import random
import subprocess
import time
from pathlib import Path


ROOT = Path('/data/user/epishkin')
BRANCHES = {
    'mouse': ('ERP003950', 'insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_ec1x3_in5ng_rb3x_ss250', 'ERR346600', 'IGH', 'mm39/star_chr12', 'vdjer_mouse/igh'),
    'human': ('PRJEB30386', 'insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_umicons_in5ng_rb3x_ss250', 'PRJEB30386_all_chains', 'IGL', 'hg38/star_chr2_chr14_chr22', 'vdjer/igl'),
}


def count_fastq(path):
    with gzip.open(path, 'rt') as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f'Incomplete FASTQ: {path}')
    return lines // 4


def sample_pairs(r1, r2, out1, out2, wanted, seed=42):
    total = count_fastq(r1)
    if wanted > total:
        raise ValueError(f'Asked for {wanted} pairs, only {total} available')
    rng = random.Random(seed)
    remaining = wanted
    with gzip.open(r1, 'rt') as a, gzip.open(r2, 'rt') as b, out1.open('w') as x, out2.open('w') as y:
        for i in range(total):
            left = [a.readline() for _ in range(4)]
            right = [b.readline() for _ in range(4)]
            if not all(left + right) or left[0].split()[0].removesuffix('/1') != right[0].split()[0].removesuffix('/2'):
                raise ValueError(f'FASTQ pair mismatch at record {i}')
            if rng.random() < remaining / (total - i):
                x.writelines(left)
                y.writelines(right)
                remaining -= 1
        if a.readline() or b.readline() or remaining:
            raise ValueError('FASTQ lengths differ or sampling failed')
    return total


def select_template_pairs(r1, r2, out1, out2, template, limit):
    selected = 0
    with gzip.open(r1, 'rt') as a, gzip.open(r2, 'rt') as b, out1.open('w') as x, out2.open('w') as y:
        while selected < limit:
            left = [a.readline() for _ in range(4)]
            right = [b.readline() for _ in range(4)]
            if not left[0] and not right[0]:
                break
            if not all(left + right) or left[0].split()[0].removesuffix('/1') != right[0].split()[0].removesuffix('/2'):
                raise ValueError('FASTQ pair mismatch')
            tid = left[0].split()[0][1:].removesuffix('/1').rsplit('_rc_', 1)[0]
            if tid == template:
                x.writelines(left)
                y.writelines(right)
                selected += 1
    return selected


def run(command, log, *, cwd=None, stdout=None):
    started = time.time()
    with log.open('w') as err:
        if stdout:
            with stdout.open('w') as out:
                result = subprocess.run(command, cwd=cwd, stdout=out, stderr=err)
        else:
            result = subprocess.run(command, cwd=cwd, stdout=err, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'{command[0]} exit {result.returncode}; see {log}')
    return round(time.time() - started, 1)


def fasta_lengths(path):
    lengths = []
    length = None
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if length is not None:
                    lengths.append(length)
                length = 0
            elif length is not None:
                length += len(line.strip())
    if length is not None:
        lengths.append(length)
    return lengths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species', choices=BRANCHES, required=True)
    parser.add_argument('--pairs', type=int, default=300_000)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--rf', type=int, choices=(0, 1), action='append', help='Run only selected coverage mode(s)')
    parser.add_argument('--template', help='Positive control: use reads from one known template only')
    args = parser.parse_args()
    if args.pairs <= 0 or args.threads <= 0:
        parser.error('--pairs and --threads must be positive')

    dataset, branch, sample, chain, index_name, ref_name = BRANCHES[args.species]
    fqdir = ROOT / 'results' / dataset / 'simulated' / branch / '06_fastq_pe150'
    r1, r2 = (fqdir / f'{sample}_R{i}.fastq.gz' for i in (1, 2))
    index = ROOT / 'seq_refs' / index_name
    ref = ROOT / 'seq_refs' / ref_name
    vdjer = ROOT / 'tools/vdjer-amplicon/vdjer'
    out_name = f'{args.species}_{sample}_{args.pairs}' + (f'_{args.template}' if args.template else '')
    out = ROOT / 'vdjer_pilots_20261006' / out_name
    out.mkdir(parents=True, exist_ok=True)
    os.environ['PATH'] = str(ROOT / 'conda/envs/bcr_env/bin') + os.pathsep + os.environ['PATH']
    for path in (r1, r2, index / 'SA', ref / 'v_index', ref / 'j_index', vdjer):
        if not path.exists():
            raise FileNotFoundError(path)
    summary = {'species': args.species, 'sample': sample, 'branch': branch, 'requested_pairs': args.pairs,
               'seed': 42, 'threads': args.threads, 'template_control': args.template, 'output': str(out), 'modes': {}}
    sub1, sub2 = out / 'sub_R1.fastq', out / 'sub_R2.fastq'
    if not (sub1.exists() and sub2.exists()):
        if args.template:
            summary['selected_pairs'] = select_template_pairs(r1, r2, sub1, sub2, args.template, args.pairs)
        else:
            summary['source_pairs'] = sample_pairs(r1, r2, sub1, sub2, args.pairs)
    summary['sampled_pairs'] = sum(1 for _ in sub1.open()) // 4
    if not summary['sampled_pairs'] or (not args.template and summary['sampled_pairs'] != args.pairs):
        raise ValueError('Incomplete pilot subset')
    print(json.dumps({'stage': 'sampled', **summary}, ensure_ascii=False), flush=True)

    bam = out / 'star.sort.bam'
    if not Path(str(bam) + '.bai').exists():
        prefix = str(out / 'star_')
        command = ['STAR', '--runThreadN', str(args.threads), '--genomeDir', str(index),
                   '--readFilesIn', str(sub1), str(sub2), '--outSAMtype', 'BAM', 'SortedByCoordinate',
                   '--outSAMunmapped', 'Within', '--limitBAMsortRAM', '4000000000', '--outFileNamePrefix', prefix]
        summary['star_seconds'] = run(command, out / 'star.log')
        (out / 'star_Aligned.sortedByCoord.out.bam').rename(bam)
        run(['samtools', 'index', str(bam)], out / 'samtools_index.log')
    print(json.dumps({'stage': 'aligned', 'bam': str(bam)}, ensure_ascii=False), flush=True)

    for rf in (args.rf or (1, 0)):
        mode = out / f'rf{rf}'
        mode.mkdir(exist_ok=True)
        contigs = mode / 'vdj_contigs.fa'
        if not contigs.exists():
            command = [str(vdjer), '--in', str(bam), '--t', str(args.threads), '--ins', '250',
                       '--chain', chain]
            if args.species == 'mouse':
                command += ['--vr', 'chr12:113218450-115973713', '--cr', 'chr12:113218450-113490429']
            command += ['--ref-dir', str(ref), '-jext', '60', '--ws', '364' if args.species == 'mouse' else '384',
                        '--e0', '20', '--e1', '363' if args.species == 'mouse' else '383', '--rf', str(rf)]
            summary['modes'][f'rf{rf}'] = {'command': command, 'seconds': run(command, mode / 'vdjer.log', cwd=mode, stdout=mode / 'vdjer.sam')}
        lengths = fasta_lengths(contigs)
        summary['modes'].setdefault(f'rf{rf}', {}).update({'contigs': len(lengths), 'length_min': min(lengths) if lengths else None,
                                                           'length_max': max(lengths) if lengths else None})
        (out / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
        print(json.dumps({'stage': f'rf{rf}', **summary['modes'][f'rf{rf}']}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
