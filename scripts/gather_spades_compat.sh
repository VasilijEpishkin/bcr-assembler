#!/usr/bin/env bash
# Обёртка rnaSPAdes для GATHeR: GATHeR читает контиги из transcripts.fasta. Если rnaSPAdes не оценил длину вставки
# и не создал этот файл, берутся финальные контиги последнего k (K*/final_contigs.fasta).
set -euo pipefail

real_spades=${GATHER_SPADES_REAL:-/data/user/epishkin/conda/envs/gather_runtime/bin/spades.py}
"$real_spades" "$@"

output_dir=
previous=
for argument in "$@"; do
    if [[ $previous == -o || $previous == --output ]]; then
        output_dir=$argument
        break
    fi
    previous=$argument
done

[[ -n $output_dir ]] || exit 0
[[ -s $output_dir/transcripts.fasta ]] && exit 0

fallback=$(find "$output_dir" -mindepth 2 -maxdepth 2 -type f \
    -path '*/K*/final_contigs.fasta' -size +0c -print | sort -V | tail -n 1)
[[ -n $fallback ]] || {
    echo "SPAdes завершился без transcripts.fasta и K*/final_contigs.fasta" >&2
    exit 1
}
cp "$fallback" "$output_dir/transcripts.fasta"
echo "GATHeR: в качестве transcripts.fasta использован $fallback" >&2
