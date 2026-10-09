#!/usr/bin/env bash
# Копирует данные основной линии с персонального тома на общий том в структуре репозитория:
#   <общий том>/bcr-assembler/results/<датасет>/simulated/<ветка>/...
# Ноутбуки, README и отчёты QC берутся из репозитория (git archive), здесь копируются только данные:
# эталон, риды, QC симуляции, таблицы бенчмарка по контигам и шаблонам, итоговые файлы сборщиков,
# post_annotation_filtered, отчёты и сводки качества аннотации, проблемные последовательности и сводки оценок. Рабочие каталоги бенчмарка, промежуточные файлы сборщиков,
# распределение ридов, аннотация и сырые FASTQ остаются на персональном томе.
#
# Все скопированные файлы и каталоги получают права на запись для всех пользователей тома (общая работа над проектом).
#
# Использование: scripts/sync_to_shared_volume.sh <каталог bcr-assembler на общем томе> [--dry-run]
set -euo pipefail
umask 000

SRC=${BCR_VOLUME:-/data/user/epishkin}/results
DST=${1:?нужен каталог bcr-assembler на общем томе}/results
shift || true
RSYNC=${RSYNC:-rsync}
S=insilicoseq_150bp_novaseq_post_annotation_filtered_random_cut_amp_
BRANCHES=(
  "ERP003950/simulated/${S}ec1x3_in5ng_rb3x_ss250"
  "ERP003950/simulated/${S}ec1x3_tails_c1400_in5ng_rb3x_ss250"
  "PRJEB30386/simulated/${S}umicons_in5ng_rb3x_ss250"
  "PRJEB30386/simulated/${S}umicons_tails_c1400_in5ng_rb3x_ss250"
)

FILTER=$(mktemp)
trap 'rm -f "$FILTER"' EXIT
cat > "$FILTER" <<'EOF'
+ /00_primary_truth/***
+ /01_pcr1/***
+ /06_fastq_pe150/***
+ /qc/***
+ /validation/***
+ /logs/***
+ /identifiability_map_*/***
+ /benchmark_unique_trim*/
+ /benchmark_unique_trim*/assemblers_*.tsv
+ /benchmark_unique_trim*/*/
+ /benchmark_unique_trim*/*/*/
+ /benchmark_unique_trim*/*/*/per_template.tsv
+ /benchmark_unique_trim*/*/*/per_contig.tsv
+ /assemblers/
+ /assemblers/*/
+ /assemblers/*/*/
+ /assemblers/*/*/contigs.fasta
+ /assemblers/*/*/*_airr.tsv
+ /assemblers/*/*/run.log
+ /assemblers/*/*/run_metadata.json
- *
EOF

for b in "${BRANCHES[@]}"; do
  mkdir -p "$DST/$b"
  echo "[$(date +%F\ %T)] $b"
  "$RSYNC" -a --chmod=a+rwX "$@" --filter="merge $FILTER" "$SRC/$b/" "$DST/$b/"
done
for d in ERP003950 PRJEB30386; do
  echo "[$(date +%F\ %T)] $d/post_annotation_filtered"
  mkdir -p "$DST/$d/post_annotation_filtered"
  "$RSYNC" -a --chmod=a+rwX "$@" "$SRC/$d/post_annotation_filtered/" "$DST/$d/post_annotation_filtered/"
done
for d in ERP003950 PRJEB30386; do
  for sub in annotation/quality annotation/report problematic_sequences; do
    [[ -d $SRC/$d/$sub ]] || continue
    echo "[$(date +%F\ %T)] $d/$sub"
    mkdir -p "$DST/$d/$sub"
    "$RSYNC" -a --chmod=a+rwX --exclude='._*' "$@" "$SRC/$d/$sub/" "$DST/$d/$sub/"
  done
done
echo "[$(date +%F\ %T)] evaluation_history"
mkdir -p "$DST/evaluation_history"
"$RSYNC" -a --chmod=a+rwX "$@" "$SRC/evaluation_history/" "$DST/evaluation_history/"
if [[ " $* " != *" --dry-run "* ]]; then
  chmod -R a+rwX "$(dirname "$DST")"
fi
echo "[$(date +%F\ %T)] готово"
