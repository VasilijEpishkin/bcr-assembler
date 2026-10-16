#!/usr/bin/env python3
"""Оценка пилота V'DJer: точные полные V…J и junction эталона среди контигов (для шаблонов, чьи риды попали в подвыборку)."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from run_vdjer_pilot import BRANCHES, ROOT


def fasta(path):
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, ''.join(chunks).upper()


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def sequence_index(contigs, k):
    oriented = []
    seeds = defaultdict(set)
    for seq in contigs:
        for strand in (seq, reverse_complement(seq)):
            cid = len(oriented)
            oriented.append(strand)
            for i in range(len(strand) - k + 1):
                seeds[strand[i:i+k]].add(cid)
    return oriented, seeds


def find_contigs(seq, oriented, seeds, k):
    if len(seq) < k:
        return {i // 2 for i, contig in enumerate(oriented) if seq in contig}
    return {i // 2 for i in seeds.get(seq[:k], ()) if seq in oriented[i]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--species', choices=BRANCHES, required=True)
    parser.add_argument('--pairs', type=int, default=300_000)
    parser.add_argument('--rf', type=int, choices=(0, 1), action='append', help='Score only selected mode(s)')
    parser.add_argument('--template', help='Positive-control template ID used in the pilot')
    args = parser.parse_args()
    dataset, branch, sample, chain, _, _ = BRANCHES[args.species]
    base = ROOT / 'results' / dataset / 'simulated' / branch
    name = f'{args.species}_{sample}_{args.pairs}' + (f'_{args.template}' if args.template else '')
    pilot = ROOT / 'vdjer_pilots_20261006' / name
    represented = set()
    with (pilot / 'sub_R1.fastq').open() as handle:
        for line in handle:
            if line.startswith('@'):
                represented.add(line.split()[0][1:].removesuffix('/1').rsplit('_rc_', 1)[0])
            for _ in range(3):
                if not handle.readline():
                    raise ValueError('Incomplete pilot FASTQ')
    qc = base / '00_primary_truth' / f'{sample}_template_qc.tsv'
    junctions = set()
    eligible = set()
    with qc.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            tid = row['template_id']
            if tid in represented and row['locus'] == chain and row['verified'].lower() == 'true':
                eligible.add(tid)
                if row['junction']:
                    junctions.add(row['junction'].upper())
    trim = 20 if args.species == 'mouse' else 0
    truth_path = base / f'benchmark_unique_trim{trim}' / 'trust4' / sample / 'work/truth_vj.fasta'
    unique_vj = {seq for tid, seq in fasta(truth_path) if tid in eligible}
    report = {'species': args.species, 'sampled_pairs': sum(1 for _ in (pilot / 'sub_R1.fastq').open()) // 4,
              'template_control': args.template, 'represented_templates': len(represented),
              'represented_verified_templates': len(eligible), 'represented_verified_unique_vj': len(unique_vj),
              'represented_verified_unique_junctions': len(junctions), 'modes': {}}
    for rf in (args.rf or (1, 0)):
        path = pilot / f'rf{rf}' / 'vdj_contigs.fa'
        contigs = [seq for _, seq in fasta(path)]
        oriented, seeds = sequence_index(contigs, 25)
        vj_hits = [find_contigs(seq, oriented, seeds, 25) for seq in unique_vj]
        junction_hits = [find_contigs(seq, oriented, seeds, 25) for seq in junctions]
        report['modes'][f'rf{rf}'] = {
            'contigs': len(contigs),
            'exact_full_unique_vj': sum(bool(hit) for hit in vj_hits),
            'contigs_containing_exact_full_vj': len(set().union(*vj_hits)) if vj_hits else 0,
            'exact_junctions_in_sequence': sum(bool(hit) for hit in junction_hits),
        }
    output = pilot / 'score.json'
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
