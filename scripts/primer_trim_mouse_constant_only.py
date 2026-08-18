#!/usr/bin/env python3
from pathlib import Path
import argparse
import shutil
import subprocess

THREADS = 1
ERROR_RATE = 0.2
MIN_OVERLAP = 10
FORWARD_ADAPTER = 'CARKGGATRRRCHGATGGGG'
REVERSE_ADAPTER = 'CCCCATCDGYYYATCCMYTG'
CONSTANT_PRIMER_MATE = {
    'ERR346596': '1',
    'ERR346597': '1',
    'ERR346598': '1',
    'ERR346599': '2',
    'ERR346600': '1',
    'ERR346601': '1',
}


def run(cmd, logfile: Path):
    with open(logfile, 'a') as lf:
        res = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if res.returncode != 0:
        raise RuntimeError(f'Command failed: {cmd}')


def main(volume: str, dataset: str, force: bool = False):
    if dataset != 'ERP003950':
        raise ValueError(f'This runner is configured only for ERP003950, got {dataset}')

    vol = Path(volume)
    src_dir = vol / 'results' / dataset / 'trimmed' / 'fastq'
    base = vol / 'results' / dataset / 'pr_trimmed'
    out_dir = base / 'fastq'
    logs_dir = base / 'cutadapt_logs'

    if force and base.exists():
        shutil.rmtree(base)

    out_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    pairs = sorted(set(
        f.name.replace('_1.trim.fastq.gz', '').replace('_2.trim.fastq.gz', '')
        for f in src_dir.glob('*.trim.fastq.gz')
    ))
    print(f'[primer_trim_mouse_constant_only] {dataset}: {len(pairs)} pairs', flush=True)

    for sample in pairs:
        constant_mate = CONSTANT_PRIMER_MATE.get(sample)
        if constant_mate is None:
            raise RuntimeError(f'No constant-primer mate mapping for sample {sample}')
        for mate in ('1', '2'):
            src = src_dir / f'{sample}_{mate}.trim.fastq.gz'
            dest = out_dir / f'{sample}_{mate}.pr.fastq.gz'
            if not src.exists():
                raise FileNotFoundError(src)
            if dest.exists():
                print(f'  [{sample}_{mate}] already done, skip', flush=True)
                continue
            if mate != constant_mate:
                shutil.copy2(src, dest)
                print(f'  [{sample}_{mate}] copied unchanged (non-constant mate)', flush=True)
                continue
            cmd = [
                'cutadapt',
                '-j', str(THREADS),
                '--overlap', str(MIN_OVERLAP),
                '-e', str(ERROR_RATE),
                '--action=trim',
                '-g', FORWARD_ADAPTER,
                '-g', REVERSE_ADAPTER,
                '-o', str(dest),
                str(src),
            ]
            print(f'  [{sample}_{mate}] cutadapt constant-region trim ...', flush=True)
            run(cmd, logs_dir / f'{sample}_{mate}.cutadapt.log')

    n = len(list(out_dir.glob('*.pr.fastq.gz')))
    print(f'[primer_trim_mouse_constant_only] DONE: {n} pr-trimmed files', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Mouse ERP003950 primer trim: cut constant-region primer only with cutadapt')
    ap.add_argument('volume')
    ap.add_argument('dataset')
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    main(args.volume, args.dataset, force=args.force)
