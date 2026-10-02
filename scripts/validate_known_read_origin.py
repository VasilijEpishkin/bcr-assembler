"""Сравнивает симулированные риды напрямую с фрагментами-источниками, указанными в их именах.

Без перекартирования и без изменения исходных данных. К истинному и картированному фрагменту применяется
один критерий — редакционное расстояние по всему запросу, поэтому локальное выравнивание и обрезка
не дают картированному попаданию нечестного преимущества.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path


def fasta(path):
    with open(path) as handle:
        name, parts = None, []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().upper())
        if name is not None:
            yield name, "".join(parts)


def fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq, plus, qual = handle.readline().strip(), handle.readline(), handle.readline().strip()
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Invalid FASTQ record in {path}")
            yield header[1:].split()[0].removesuffix("/1").removesuffix("/2"), seq, qual


def fragment_id(name):
    fid, index, cpu = name.rsplit("_", 2)
    if not index.isdigit() or not cpu.isdigit() or "_rc_" not in fid:
        raise ValueError(f"Unrecognised simulated read name: {name}")
    return fid


def run(branch, sample, r1, r2, bam, output):
    import edlib
    import pysam
    branch, output = Path(branch), Path(output)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    reads = list(fastq(r1))
    mates = list(fastq(r2))
    assert len(reads) == len(mates) and len(reads) > 0
    names = {a[0] for a in reads}
    assert len(names) == len(reads), "Duplicate read-pair IDs"
    origins = {n: fragment_id(n) for n in names}
    needed = set(origins.values())
    mapped = {}
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            name = read.query_name.removesuffix("/1").removesuffix("/2")
            if name not in names or read.is_unmapped:
                continue
            key = name, 2 if read.is_read2 else 1
            assert key not in mapped, "Multiple primary alignments for one mate"
            mapped[key] = (read.reference_name, read.mapping_quality, read.get_tag("NM") if read.has_tag("NM") else None)
            needed.add(read.reference_name)
    refpath = branch / "04_read_allocation" / f"{sample}_selected_fragments.fasta"
    fragments = {n: s for n, s in fasta(refpath) if n in needed}
    assert needed <= fragments.keys(), "Missing source or alignment reference"
    print(f"Loaded {len(reads)} pairs, {len(fragments)} required fragments", flush=True)
    rc_table = str.maketrans("ACGTN", "TGCAN")
    def revcomp(seq):
        return seq.translate(rc_table)[::-1]
    def hw(query, reference):
        return edlib.align(query, reference, mode="HW", task="distance")["editDistance"]
    counters, hist, end_hist, pair_hist = Counter(), Counter(), Counter(), Counter()
    alt_comparison, pair_comparison, control_hist = Counter(), Counter(), Counter()
    mapped_mapq, locus_metrics, suspicious = Counter(), {}, []
    n_bases = n_edits = n_end_edits = 0
    qual_expected_errors = 0.0
    for idx, (a, b) in enumerate(zip(reads, mates)):
        name, s1, q1 = a
        assert b[0] == name
        s2, q2 = b[1:]
        fid = origins[name]
        ref = fragments[fid]
        locus = fid.split("_tpl_", 1)[0].rsplit("_", 1)[-1]
        lm = locus_metrics.setdefault(locus, Counter())
        own_pair, assigned_pair, mapped_refs = 0, 0, []
        for mate, seq, qual in ((1, s1, q1), (2, s2, q2)):
            counters["reads"] += 1
            oriented_ref = ref if mate == 1 else revcomp(ref)
            # SHW закрепляет ожидаемое начало и допускает небольшой сдвиг конца
            # из-за инделей до обрезки PE151->PE150. HW не ограничен
            # внутри истинного фрагмента и используется для честного сравнения попаданий.
            end_ed = edlib.align(seq, oriented_ref[:len(seq) + 12], mode="SHW", task="distance")["editDistance"]
            own_ed = hw(seq, oriented_ref)
            hist[own_ed] += 1
            end_hist[end_ed] += 1
            own_pair += own_ed
            n_bases += len(seq)
            n_edits += own_ed
            n_end_edits += end_ed
            qual_expected_errors += sum(10 ** (-(ord(q) - 33) / 10) for q in qual)
            lm["reads"] += 1
            lm["source_edits"] += own_ed
            lm["bases"] += len(seq)
            lm["reads_source_edits_le_3"] += own_ed <= 3
            hit = mapped.get((name, mate))
            if hit is None:
                counters["no_primary_mapped_hit"] += 1
            else:
                hitid, mapq, nm = hit
                mapped_mapq[mapq] += 1
                alt = fragments[hitid]
                # Проверяются обе цепи: произвольный выбор ориентации картировщиком
                # для неоднозначных попаданий не штрафуется.
                alt_ed = min(hw(seq, alt), hw(seq, revcomp(alt)))
                assigned_pair += alt_ed
                mapped_refs.append(hitid)
                counters["mapped_source_id_equal"] += hitid == fid
                counters["mapped_source_sequence_equal"] += alt == ref
                if hitid != fid:
                    relation = "source_better" if own_ed < alt_ed else "equal_full_read_edits" if own_ed == alt_ed else "assigned_better"
                    alt_comparison[relation] += 1
                    counters["different_id_but_error_free_source"] += own_ed == 0
                    counters["different_id_and_mapq_le_1"] += mapq <= 1
            if own_ed > 5 and len(suspicious) < 100:
                suspicious.append({"read": name, "mate": mate, "source_fragment": fid, "source_edits": own_ed, "end_anchored_edits": end_ed})
        pair_hist[own_pair] += 1
        if len(mapped_refs) == 2 and any(x != fid for x in mapped_refs):
            key = "source_better" if own_pair < assigned_pair else "equal_full_pair_edits" if own_pair == assigned_pair else "assigned_better"
            pair_comparison[key] += 1
            counters["mapped_mates_on_different_fragments"] += mapped_refs[0] != mapped_refs[1]
        # Детерминированный отрицательный контроль на 1% пар по другому шаблону.
        if idx % 100 == 0:
            candidate = origins[reads[(idx + 7919) % len(reads)][0]]
            if candidate.split("_rc_", 1)[0] != fid.split("_rc_", 1)[0]:
                wrong = fragments[candidate]
                wrong_ed = min(hw(s1, wrong), hw(s1, revcomp(wrong)))
                control_hist[wrong_ed] += 1
        if idx and idx % 50000 == 0:
            print(f"Checked {idx} pairs", flush=True)
    n = counters["reads"]
    result = {
        "branch": str(branch), "sample": sample, "pairs_checked": len(reads),
        "method": "Full-query Levenshtein distance to known source; end-anchored SHW and whole-fragment HW; identical HW metric for alternative mapped fragments, both strands.",
        "inputs": {"R1": str(r1), "R2": str(r2), "bam": str(bam), "source_fasta": str(refpath)},
        "counters": dict(counters), "source_edit_histogram": dict(sorted(hist.items())),
        "source_end_anchored_edit_histogram": dict(sorted(end_hist.items())),
        "source_pair_edit_histogram": dict(sorted(pair_hist.items())),
        "source_read_edit_rate": n_edits / n_bases,
        "source_end_anchored_edit_rate": n_end_edits / n_bases,
        "phred_implied_error_rate": qual_expected_errors / n_bases,
        "source_fraction_reads_le_3_edits": sum(v for k, v in hist.items() if k <= 3) / n,
        "source_fraction_reads_le_5_edits": sum(v for k, v in hist.items() if k <= 5) / n,
        "different_mapped_id_comparison": dict(alt_comparison),
        "different_mapped_pair_comparison": dict(pair_comparison),
        "mapped_mapq_histogram": dict(sorted(mapped_mapq.items())),
        "negative_control_other_template_edit_histogram": dict(sorted(control_hist.items())),
        "by_locus": {k: dict(v) for k, v in locus_metrics.items()},
        "first_high_error_examples": suspicious,
        "limitations": "Edit-distance ties establish sequence ambiguity, not unique molecular provenance. Different fragments can explain the same observed reads. The error-rate metric includes substitutions and indels and is not a calibrated Phred test."}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if "histogram" not in k and k != "first_high_error_examples"}, indent=2), flush=True)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("branch", "sample", "r1", "r2", "bam", "output"):
        p.add_argument("--" + key, required=True)
    run(**vars(p.parse_args()))
