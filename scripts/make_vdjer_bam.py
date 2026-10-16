#!/usr/bin/env python3
"""Готовит BAM для V'DJer из симулированных ридов PE150: подвыборка пар → STAR → сортированный индексированный BAM.

Результат лежит в `<out>/<вид>_<образец>_<пары>/` под теми же именами, которые ждут `run_vdjer_pilot.py`
и `score_vdjer_pilot.py` (`sub_R1.fastq`, `sub_R2.fastq`, `star.sort.bam`, `star.sort.bam.bai`), поэтому пилот
использует готовый BAM и не пересчитывает выравнивание. Описание шагов и состава BAM — `docs/vdjer_bam.md`.

Пример (человек, 300 тыс. пар):
    python scripts/make_vdjer_bam.py --species human --pairs 300000 --threads 8
"""
import argparse
import json
import os
from pathlib import Path

from run_vdjer_pilot import BRANCHES, ROOT, run, sample_pairs

DEFAULT_OUT = ROOT / 'vdjer_pilots_20261006'


def star_stats(log_final):
    """Строки Log.final.out STAR в виде словаря {название: значение}."""
    stats = {}
    for line in Path(log_final).read_text().splitlines():
        if '|' in line:
            key, value = line.split('|', 1)
            stats[key.strip()] = value.strip()
    return stats


def samtools_lines(command, log):
    run(command, log, stdout=log)
    return log.read_text().strip().splitlines()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--species', choices=BRANCHES, required=True)
    parser.add_argument('--pairs', type=int, default=300_000, help='число пар в подвыборке (seed 42)')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT, help='каталог пилотов')
    args = parser.parse_args()
    if args.pairs <= 0 or args.threads <= 0:
        parser.error('--pairs и --threads должны быть положительными')

    dataset, branch, sample, _chain, index_name, _ref = BRANCHES[args.species]
    fqdir = ROOT / 'results' / dataset / 'simulated' / branch / '06_fastq_pe150'
    r1, r2 = (fqdir / f'{sample}_R{i}.fastq.gz' for i in (1, 2))
    index = ROOT / 'seq_refs' / index_name
    for path in (r1, r2, index / 'SA'):
        if not path.exists():
            raise FileNotFoundError(path)
    os.environ['PATH'] = str(ROOT / 'conda/envs/bcr_env/bin') + os.pathsep + os.environ['PATH']

    out = args.out / f'{args.species}_{sample}_{args.pairs}'
    out.mkdir(parents=True, exist_ok=True)
    bam = out / 'star.sort.bam'
    if Path(str(bam) + '.bai').exists():
        print('BAM уже готов:', bam)
        return

    sub1, sub2 = out / 'sub_R1.fastq', out / 'sub_R2.fastq'
    source_pairs = sample_pairs(r1, r2, sub1, sub2, args.pairs)
    star = ['STAR', '--runThreadN', str(args.threads), '--genomeDir', str(index),
            '--readFilesIn', str(sub1), str(sub2), '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--outSAMunmapped', 'Within', '--limitBAMsortRAM', '4000000000', '--outFileNamePrefix', str(out / 'star_')]
    star_seconds = run(star, out / 'star.log')
    (out / 'star_Aligned.sortedByCoord.out.bam').rename(bam)
    run(['samtools', 'index', str(bam)], out / 'samtools_index.log')

    provenance = {
        'species': args.species, 'sample': sample, 'branch': branch,
        'source_fastq': [str(r1), str(r2)], 'source_pairs': source_pairs, 'sampled_pairs': args.pairs, 'seed': 42,
        'genome_index': str(index), 'star_command': star, 'star_seconds': star_seconds,
        'star_stats': star_stats(out / 'star_Log.final.out'),
        'samtools_flagstat': samtools_lines(['samtools', 'flagstat', str(bam)], out / 'flagstat.txt'),
        'samtools_idxstats': samtools_lines(['samtools', 'idxstats', str(bam)], out / 'idxstats.txt'),
    }
    (out / 'bam_provenance.json').write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + '\n')
    print('BAM:', bam)


if __name__ == '__main__':
    main()
