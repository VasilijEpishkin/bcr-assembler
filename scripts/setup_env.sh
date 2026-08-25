#!/usr/bin/env bash
# setup_env.sh — единый скрипт окружения OneQ VM (jupyter/base-notebook).
# Объединяет старые setup_vm_conda.sh (сборка bcr_env с нуля) и env_setup.sh
# (быстрая активация PATH в терминале).
#
# ЧТО ДЕЛАЕТ:
#   1. Если bcr_env отсутствует или неполон — создаёт/дополняет его:
#      fastqc, fastp, cutadapt, multiqc, presto, bowtie2, samtools,
#      insilicoseq, rsync (conda-forge/bioconda + pip).
#   2. Регистрирует ipykernel "BCR Pipeline" для Jupyter (PATH зашит в kernel.json).
#   3. Экспортирует PATH в ТЕКУЩИЙ шелл (нужно source, не bash) — быстрая
#      активация без переустановки, если всё уже стоит.
#
# ПОЧЕМУ ТАК (а не в ~/.local/bin + apt):
#   - На OneQ jupyter/base-notebook НЕТ sudo, apt недоступен.
#   - Из системных инструментов есть micromamba (/data/user/epishkin/bin_micromamba),
#     используем его — без активации conda/conda.sh, без shell-хуков.
#   - SSH-шелл на OneQ имеет более узкий cgroup (по памяти), чем Jupyter kernel —
#     тяжёлые install'ы (особенно conda-forge shard index) может убить OOM
#     через SSH. Запускай этот скрипт из Jupyter Terminal, не через SSH.
#
# ЗАПУСК:
#   В Jupyter Terminal (НЕ в обычной ячейке ноутбука):
#     source /data/user/epishkin/scripts/setup_env.sh
#   (просто "bash setup_env.sh" тоже отработает установку, но тогда PATH
#   останется только в дочернем процессе — если нужен PATH в текущем
#   терминале, используй именно "source").
#
# ПОСЛЕ ЗАПУСКА:
#   - В Jupyter: Kernel -> Change Kernel -> "BCR Pipeline"
#   - Или в ячейке ноутбука вручную:
#       import os
#       os.environ["PATH"] = "/data/user/epishkin/conda/envs/bcr_env/bin:" + os.environ.get("PATH", "")

ENV_NAME="bcr_env"
ENV_PREFIX="/data/user/epishkin/conda/envs/${ENV_NAME}"
MICROMAMBA="/data/user/epishkin/bin_micromamba"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
log()  { echo -e "${GREEN}[$(date -u +%H:%M:%S)]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }

_setup_needed=0
for tool in fastp fastqc cutadapt multiqc MaskPrimers.py bowtie2 bowtie2-build samtools iss rsync; do
    if [ ! -x "${ENV_PREFIX}/bin/${tool}" ]; then
        _setup_needed=1
        break
    fi
done

if [ "${_setup_needed}" -eq 0 ]; then
    log "bcr_env уже содержит все нужные инструменты -> ${ENV_PREFIX}, пропускаю установку"
else
    log "=== Установка/дополнение bcr_env через micromamba ==="
    if [ ! -x "${MICROMAMBA}" ]; then
        warn "micromamba не найден по ${MICROMAMBA} — поправь путь в скрипте"
        return 1 2>/dev/null || exit 1
    fi

    log "conda-forge/bioconda: fastqc fastp bowtie2 samtools rsync"
    "${MICROMAMBA}" create -y -p "${ENV_PREFIX}" -c bioconda -c conda-forge \
        python=3.11 fastqc fastp bowtie2 samtools rsync \
        || { warn "create failed, пробую install в уже существующий env"; \
             "${MICROMAMBA}" install -y -p "${ENV_PREFIX}" -c bioconda -c conda-forge \
                 fastqc fastp bowtie2 samtools rsync; }

    log "pip: cutadapt multiqc presto insilicoseq ipykernel"
    "${ENV_PREFIX}/bin/pip" install --upgrade pip setuptools wheel
    "${ENV_PREFIX}/bin/pip" install cutadapt multiqc presto insilicoseq ipykernel

    log "=== Регистрация Jupyter kernel ==="
    "${ENV_PREFIX}/bin/python3" -m ipykernel install --user --name "${ENV_NAME}" --display-name "BCR Pipeline"
    log "  kernel 'BCR Pipeline' зарегистрирован"

    KJSON=$("${ENV_PREFIX}/bin/python3" -c "import jupyter_core,os; \
print(os.path.join(jupyter_core.paths.jupyter_data_dir(),'kernels','${ENV_NAME}','kernel.json'))")
    log "  kernel.json: ${KJSON}"
    if [ -f "${KJSON}" ]; then
        "${ENV_PREFIX}/bin/python3" - "$KJSON" "${ENV_PREFIX}" <<'PY'
import json, sys
kjson, envp = sys.argv[1], sys.argv[2]
d = json.load(open(kjson))
d["env"] = {"PATH": f"{envp}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"}
json.dump(d, open(kjson, "w"), indent=2)
print("  [OK] env.PATH прописан в kernel.json")
PY
    else
        warn "  kernel.json не найден, env-cell в ноутбуке остаётся обязательным"
    fi

    log "=== Финальная проверка ==="
    for tool in fastp fastqc cutadapt multiqc MaskPrimers.py bowtie2 samtools iss rsync; do
        if [ -x "${ENV_PREFIX}/bin/${tool}" ]; then
            echo "  [OK] ${tool}"
        else
            warn "  [MISS] ${tool}"
        fi
    done
fi

# Активация PATH в текущем шелле (работает только если скрипт был source'нут).
export PATH="/data/user/epishkin/venvs/presto_env/bin:${ENV_PREFIX}/bin:${PATH}"
log "PATH обновлён (bcr_env первым). Проверь: which iss / which bowtie2"
