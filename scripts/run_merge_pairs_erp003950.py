#!/usr/bin/env python3
import csv, gzip, shutil, subprocess, time
from pathlib import Path

VOLUME = '/data/user/epishkin'
DATASET = 'ERP003950'
LABEL = 'pr_trimmed'
NPROC = 8
GZIP_OUTPUT = True
FORCE = False
WRITE_FAILED_READS = False
WRITE_VERBOSE_LOG = False


def count_fastq_records(path: Path) -> int:
    if not path.exists():
        return 0
    opener = gzip.open if path.suffix == '.gz' else open
    line_count = 0
    with opener(path, 'rt') as handle:
        for line_count, _ in enumerate(handle, start=1):
            pass
    return line_count // 4


def run_merge_pairs_optimized(volume, dataset, label='pr_trimmed', nproc=8,
                              write_verbose_log=False, write_failed_reads=False,
                              gzip_output=True, force=False):
    vol = Path(volume)
    suffix1 = '_1.pr.fastq.gz' if label == 'pr_trimmed' else '_1.trim.fastq.gz'
    suffix2 = '_2.pr.fastq.gz' if label == 'pr_trimmed' else '_2.trim.fastq.gz'
    src = vol / 'results' / dataset / label / 'fastq'
    out_base = vol / 'results' / dataset / 'merged'
    out_fastq = out_base / 'fastq'
    out_logs = out_base / 'logs'
    out_qc = out_base / 'qc'
    out_fastq.mkdir(parents=True, exist_ok=True)
    out_logs.mkdir(parents=True, exist_ok=True)
    out_qc.mkdir(parents=True, exist_ok=True)

    pairs = []
    for r1 in sorted(src.glob(f'*{suffix1}')):
        sample = r1.name[:-len(suffix1)]
        r2 = src / f'{sample}{suffix2}'
        if r2.exists():
            pairs.append((sample, r1, r2))

    print(f'[merge_pairs_optimized] dataset={dataset} label={label} pairs={len(pairs)} nproc={nproc}')
    rows = []
    for idx, (sample, r1, r2) in enumerate(pairs, start=1):
        pass_fastq = out_fastq / f'{sample}_assemble-pass.fastq.gz'
        fail1_fastq = out_fastq / f'{sample}_assemble-fail-1.fastq.gz'
        fail2_fastq = out_fastq / f'{sample}_assemble-fail-2.fastq.gz'
        stdout_path = out_logs / f'{sample}.stdout.txt'
        stderr_path = out_logs / f'{sample}.stderr.txt'
        log_path = out_logs / f'{sample}.assemble.log'

        if pass_fastq.exists() and not force:
            print(f'  [{idx}/{len(pairs)}] [skip] {sample}: {pass_fastq.name} exists')
            status = 'skipped'
            elapsed = 0.0
        else:
            cmd = [
                'AssemblePairs.py', 'align',
                '-1', str(r1),
                '-2', str(r2),
                '--coord', 'illumina',
                '--rc', 'tail',
                '--outname', sample,
                '--outdir', str(out_fastq),
                '--nproc', str(nproc),
            ]
            if gzip_output:
                cmd.append('--gzip-output')
            if write_failed_reads:
                cmd.append('--failed')
            if write_verbose_log:
                cmd += ['--log', str(log_path)]
            print(f'  [{idx}/{len(pairs)}] [run] {sample}')
            t0 = time.time()
            res = subprocess.run(cmd, capture_output=True, text=True)
            elapsed = time.time() - t0
            stdout_path.write_text(res.stdout)
            stderr_path.write_text(res.stderr)
            if res.returncode != 0:
                raise RuntimeError(f'AssemblePairs failed for {sample}; see {stderr_path}')
            status = 'done'

        n_pass = count_fastq_records(pass_fastq)
        n_fail = count_fastq_records(fail1_fastq) + count_fastq_records(fail2_fastq)
        total = n_pass + n_fail
        rate = n_pass / total if total else 0.0
        rows.append({
            'sample': sample,
            'status': status,
            'elapsed_sec': f'{elapsed:.1f}',
            'assembled_pairs': str(n_pass),
            'failed_reads': str(n_fail),
            'total_seen_records': str(total),
            'merge_rate_from_outputs': f'{rate:.6f}',
            'pass_fastq': str(pass_fastq),
            'fail1_fastq': str(fail1_fastq),
            'fail2_fastq': str(fail2_fastq),
            'stdout': str(stdout_path),
            'stderr': str(stderr_path),
            'verbose_log': str(log_path) if write_verbose_log else '',
        })
        qc_path = out_qc / 'assembly_qc.tsv'
        with open(qc_path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        print(f'    pass={n_pass:,} fail_records={n_fail:,} elapsed={elapsed/60:.1f} min')

    total_pass = sum(int(r['assembled_pairs']) for r in rows)
    total_fail = sum(int(r['failed_reads']) for r in rows)
    print(f"[merge_pairs_optimized] qc: {out_qc / 'assembly_qc.tsv'}")
    print(f'[merge_pairs_optimized] assembled_pairs={total_pass:,}')
    print(f'[merge_pairs_optimized] failed_records={total_fail:,}')

if __name__ == '__main__':
    run_merge_pairs_optimized(VOLUME, DATASET, label=LABEL, nproc=NPROC,
                              write_verbose_log=WRITE_VERBOSE_LOG,
                              write_failed_reads=WRITE_FAILED_READS,
                              gzip_output=GZIP_OUTPUT, force=FORCE)
