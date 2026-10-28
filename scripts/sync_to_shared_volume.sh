#!/usr/bin/env bash
# Копирует данные проекта из корня проекта на персональном томе на общий том в структуре репозитория:
#   <общий том>/bcr-assembler/{results,raw,seq_refs}
# Копируются все датасеты, финальная ветка симуляции целиком и таблицы бенчмарка остальных веток. Ноутбуки, README, отчёты и история оценок берутся из репозитория (git archive).
#
# Все скопированные файлы и каталоги получают права на запись для всех пользователей тома (общая работа над проектом).
# Повторный запуск докопирует изменения; файлы, которых нет на источнике, на общем томе не удаляются.
#
# Использование: scripts/sync_to_shared_volume.sh <каталог bcr-assembler на общем томе> [--dry-run]
set -euo pipefail
umask 000

SRC=${BCR_VOLUME:-/data/user/epishkin/bcr-assembler}
DST=${1:?нужен каталог bcr-assembler на общем томе}
shift || true
RSYNC=${RSYNC:-rsync}
COMMON=(-a --chmod=a+rwX --exclude='._*' --exclude='.ipynb_checkpoints' --exclude='*_tails_*')

mkdir -p "$DST/results" "$DST/raw" "$DST/seq_refs"
echo "[$(date +%F\ %T)] results"
# Целиком копируется только финальная ветка симуляции (ss200-350h25); остальные ветки — только таблицы бенчмарка, на них ссылаются отчёты.
"$RSYNC" "${COMMON[@]}" \
  --include='/*/simulated/' --include='/*/simulated/*ss200-350h25*/' --include='/*/simulated/*ss200-350h25*/***' \
  --exclude='/*/simulated/*' "$@" "$SRC/results/" "$DST/results/"
"$RSYNC" -a --chmod=a+rwX --prune-empty-dirs --exclude='._*' --exclude='.ipynb_checkpoints' --exclude='*_tails_*' \
  --include='/*/' --include='/*/simulated/' --include='/*/simulated/*/' \
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
