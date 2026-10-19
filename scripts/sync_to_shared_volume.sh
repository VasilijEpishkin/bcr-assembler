#!/usr/bin/env bash
# Копирует данные проекта с персонального тома на общий том в структуре репозитория:
#   <общий том>/bcr-assembler/{results,raw,seq_refs}
# Копируются все датасеты и основные ветки симуляции целиком. Ноутбуки, README, отчёты и история оценок берутся из репозитория (git archive).
#
# Все скопированные файлы и каталоги получают права на запись для всех пользователей тома (общая работа над проектом).
# Повторный запуск докопирует изменения; файлы, которых нет на источнике, на общем томе не удаляются.
#
# Использование: scripts/sync_to_shared_volume.sh <каталог bcr-assembler на общем томе> [--dry-run]
set -euo pipefail
umask 000

SRC=${BCR_VOLUME:-/data/user/epishkin}
DST=${1:?нужен каталог bcr-assembler на общем томе}
shift || true
RSYNC=${RSYNC:-rsync}
COMMON=(-a --chmod=a+rwX --exclude='._*' --exclude='.ipynb_checkpoints' --exclude='*_tails_*' --exclude='*_sat4e12' --exclude='*in100ng_rb3x*' --exclude='*_ssnone_pcr1c30_pcr2c5_seq1e-3_al1e-1')

mkdir -p "$DST/results" "$DST/raw" "$DST/seq_refs"
echo "[$(date +%F\ %T)] results"
"$RSYNC" "${COMMON[@]}" --exclude='/evaluation_history' "$@" "$SRC/results/" "$DST/results/"
# У исключённых выше веток копируются только таблицы бенчмарка: на них ссылаются отчёты.
"$RSYNC" -a --chmod=a+rwX --prune-empty-dirs --exclude='._*' --exclude='.ipynb_checkpoints' \
  --include='/*/' --include='/*/simulated/' \
  --include='/*/simulated/*in100ng_rb3x*/' --include='/*/simulated/*_ssnone_pcr1c30_pcr2c5_seq1e-3_sat4e12/' \
  --include='/*/simulated/*/benchmark_unique_*/' --include='/*/simulated/*/benchmark_unique_*/assemblers_*.tsv' \
  --exclude='*' "$@" "$SRC/results/" "$DST/results/"
echo "[$(date +%F\ %T)] raw"
"$RSYNC" "${COMMON[@]}" --exclude='/nastyas_weird_files' "$@" "$SRC/raw/" "$DST/raw/"
echo "[$(date +%F\ %T)] seq_refs"
"$RSYNC" "${COMMON[@]}" "$@" "$SRC/seq_refs/" "$DST/seq_refs/"
if [[ " $* " != *" --dry-run "* ]]; then
  chmod -R a+rwX "$DST"
fi
echo "[$(date +%F\ %T)] готово"
