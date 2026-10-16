#!/usr/bin/env bash
# setup_env.sh — единое окружение bcr_env для всех инструментов проекта.
#
# Использование:
#   source scripts/setup_env.sh            установить недостающее и добавить окружение в PATH
#   source scripts/setup_env.sh --check    только проверка: версии, совместимость, доступные обновления
#
# Состав окружения:
#   conda (bioconda/conda-forge): fastp, fastqc, bowtie2, samtools, minimap2, seqkit, blast, bcalm, salmon, STAR,
#                                 spades, trinity, gather, insilicoseq, pysam, openjdk, R и научный стек Python;
#   pip: cutadapt, multiqc, presto, changeo, abstar, airr (пакеты, которых нет в каналах или нужны новее).
#   Инструменты без conda/pip-пакета ставятся готовыми бинарниками в $TOOLS: TRUST4, MiXCR, V'DJer, IgReC, IgBLAST.
#   Отдельных окружений нет: все ноутбуки работают с одним bcr_env.

ROOT=${BCR_VOLUME:-/data/user/epishkin}
ENV_PREFIX=$ROOT/conda/envs/bcr_env
MICROMAMBA=$ROOT/bin_micromamba
TOOLS=$ROOT/tools
IGBLAST_BIN=$ROOT/igblast/bin

# Пакеты: имя=минимальная версия (проверка сравнивает установленную версию с минимальной).
CONDA_PKGS=(
  python=3.11 fastp=1.3 fastqc=0.12 bowtie2=2.5 samtools=1.2 minimap2=2.3 seqkit=2.1 blast=2.17 bcalm=2.2
  salmon=2.8 star=2.7 spades=4.3 trinity=2.15 gather=1.0 insilicoseq=2.0 pysam=0.2 openjdk=25 r-base=4.4
  rsync=3.4 biopython=1.8 numpy=2 scipy=1.1 pandas=3 polars=1.4 matplotlib=3 ipykernel=7
)
PIP_PKGS=(cutadapt=5.2 multiqc=1.35 presto=0.7 changeo=1.3 abstar=0.8 airr=2)
# Бинарники вне conda: путь относительно $TOOLS.
BINARIES=(TRUST4/run-trust4 mixcr-4.7.0/mixcr vdjer-0.12 vdjer-amplicon igrec)

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
log()  { echo -e "${GREEN}[$(date -u +%H:%M:%S)]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }

# список установленных пакетов читается один раз за проверку
installed_version() { awk -v n="$1" '$1==n{print $2; exit}' <<< "$PKG_LIST"; }

# 0, если версия $1 не ниже минимальной $2
version_ok() { [[ -n $1 && $(printf '%s\n%s\n' "$2" "$1" | sort -V | head -1) == "$2" ]] || [[ $1 == "$2"* ]]; }

check_env() {
  local bad=0 spec name min cur
  PKG_LIST=$("$MICROMAMBA" list -p "$ENV_PREFIX" 2>/dev/null | awk 'NR>3{print $1, $2}')
  log "=== Версии и совместимость ==="
  for spec in "${CONDA_PKGS[@]}" "${PIP_PKGS[@]}"; do
    name=${spec%%=*}; min=${spec#*=}; cur=$(installed_version "$name")
    if [[ -z $cur ]]; then warn "  [MISS] $name (нужна $min+)"; bad=1
    elif version_ok "$cur" "$min"; then echo "  [OK] $name $cur"
    else warn "  [OLD] $name $cur (нужна $min+)"; bad=1; fi
  done
  log "pip check (зависимости пакетов pip):"
  "$ENV_PREFIX/bin/pip" check || bad=1
  log "Бинарники в $TOOLS:"
  for spec in "${BINARIES[@]}"; do
    if [[ -e $TOOLS/$spec ]]; then echo "  [OK] $spec"; else warn "  [MISS] $spec"; bad=1; fi
  done
  if [[ -x $IGBLAST_BIN/igblastn ]]; then echo "  [OK] igblastn $(LD_LIBRARY_PATH=$ENV_PREFIX/lib "$IGBLAST_BIN/igblastn" -version | head -1)"
  else warn "  [MISS] igblastn"; bad=1; fi
  return $bad
}

outdated() {
  log "=== Доступные обновления (ничего не меняется) ==="
  log "conda:"
  "$MICROMAMBA" update -p "$ENV_PREFIX" -c bioconda -c conda-forge --all --dry-run 2>&1 | grep -E "^\s+\S+\s+[0-9]" | head -40
  log "pip:"
  "$ENV_PREFIX/bin/pip" list --outdated 2>/dev/null | grep -i -E "^($(IFS='|'; echo "${PIP_PKGS[*]%%=*}")) "
}

install_env() {
  [[ -x $MICROMAMBA ]] || { warn "micromamba не найден: $MICROMAMBA"; return 1; }
  local specs=() pspecs=() s
  for s in "${CONDA_PKGS[@]}"; do specs+=("${s%%=*}>=${s#*=}"); done
  for s in "${PIP_PKGS[@]}"; do pspecs+=("${s%%=*}>=${s#*=}"); done
  log "conda: создание или дополнение окружения"
  if [[ -x $ENV_PREFIX/bin/python ]]; then
    "$MICROMAMBA" install -y -p "$ENV_PREFIX" -c bioconda -c conda-forge "${specs[@]}" || return 1
  else
    "$MICROMAMBA" create -y -p "$ENV_PREFIX" -c bioconda -c conda-forge "${specs[@]}" || return 1
  fi
  log "pip: пакеты вне каналов conda"
  "$ENV_PREFIX/bin/pip" install --upgrade pip setuptools wheel "${pspecs[@]}" || return 1
  log "Регистрация ядра Jupyter «BCR Pipeline»"
  "$ENV_PREFIX/bin/python" -m ipykernel install --user --name bcr_env --display-name "BCR Pipeline"
  local kjson
  kjson=$("$ENV_PREFIX/bin/python" -c "import jupyter_core.paths as p,os; print(os.path.join(p.jupyter_data_dir(),'kernels','bcr_env','kernel.json'))")
  "$ENV_PREFIX/bin/python" - "$kjson" "$ENV_PREFIX" <<'PY'
import json, sys
path, env = sys.argv[1:3]
d = json.load(open(path))
d["env"] = {"PATH": f"{env}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}
json.dump(d, open(path, "w"), indent=2)
PY
}

if [[ ${1:-} == --check ]]; then
  check_env; outdated
else
  check_env >/dev/null 2>&1 || { install_env && check_env; }
fi

export PATH="$ENV_PREFIX/bin:$TOOLS/TRUST4:$TOOLS/mixcr-4.7.0:$IGBLAST_BIN:$PATH"
log "PATH обновлён: bcr_env первым"
