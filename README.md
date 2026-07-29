# RIMA — BCR Repertoire Analysis Pipeline

Проект **rima** — пайплайн для препроцессинга и сшивания paired-end
ридов B-клеточных рецепторов (BCR) из высокопроизводительного
секвенирования (HTS) на платформе Illumina MiSeq (2×300 bp).

Пайплайн выполняет полный цикл Stage A: контроль качества сырых
данных, обрезку адаптеров и технических последовательностей,
сшивание парных ридов и подготовку данных для downstream-анализа
(IgBLAST/IMGT, клонотипирование, SHM-профилирование).

## Датасеты

| Биопроект | Организм | Пары ридов | Протокол | Статья |
|-----------|----------|------------|----------|--------|
| PRJEB40348 | Человек (Homo sapiens) | 35 | Multiplex PCR, VH/VL FR1-праймеры, MiSeq 2×300 | Lomakin et al., Front. Immunol. 2022 |
| PRJNA848968 | Лошадь (Equus caballus) | 4 | EquPD v2020, phage display scFv, MiSeq 2×300 | Rosenfeld et al., Front. Immunol. 2022 |
| PRJNA900592 | Овца (Ovis aries) | 12 | 5′ RACE (SMARTer), MiSeq 2×300 | Park et al., Mol. Immunol. 2023 |

Данные скачиваются из ENA (European Nucleotide Archive) через
`wget` в `raw/<биопроект>/` на вычислительной платформе OneQ.

Сырые FASTQ-файлы **не коммитятся** в репозиторий — они хранятся
только на вычислительном volume.

## Этапы пайплайна

### 1. Контроль качества сырых ридов (QC)

**Ноутбук:** `notebooks/qc.ipynb`
**Скрипт:** `scripts/qc.py`

```python
run_qc("/data/user/epishkin", "PRJNA900592", "raw")
```

FastQC на каждый образец + MultiQC-сводка. Результат:
`results/<DS>/qc_raw/` (FastQC-отчёты + MultiQC).

### 2. Обрезка адаптеров (adapter trim)

**Ноутбуки:**
- `notebooks/adapter_trim_human.ipynb` — человек (NEBNext/TruSeq-адаптеры)
- `notebooks/adapter_trim_mouse.ipynb` — мышь (explicit adapters, PE250-like MiSeq)
- `notebooks/adapter_trim_horse.ipynb` — лошадь (эмпирические run-specific R1-адаптеры)
- `notebooks/adapter_trim_sheep.ipynb` — овца (Illumina/NEBNext read-through + quality/length)

**Инструменты:** cutadapt (явные adapter sequences → R1/R2) → fastp
(`-q 30 -l 250 --detect_adapter_for_pe`, `--cut_right` отключён).

```python
run_adapter_trim("/data/user/epishkin", "PRJNA900592")
```

Результат: `results/<DS>/trimmed/fastq/` (.trim.fastq.gz) +
`results/<DS>/trimmed/fastp_reports/` (JSON + HTML).

### 3. QC после обрезки адаптеров

```python
run_qc("/data/user/epishkin", "PRJNA900592", "trimmed")
```

Результат: `results/<DS>/trimmed/fastqc/` + `results/<DS>/trimmed/multiqc/`.

### 4. Обрезка праймеров / технических последовательностей (primer trim)

**Ноутбуки:**
- `notebooks/primer_trim_human.ipynb` — человек (V-праймеры через pRESTO MaskPrimers)
- `notebooks/primer_trim_mouse.ipynb` — мышь (Greiff 2014 multiplex primers через MaskPrimers)
- `notebooks/primer_trim.ipynb` — legacy notebook для human/horse
- `notebooks/primer_trim_sheep.ipynb` — овца (5′ RACE: SMARTer anchor + C-region reverse primers через cutadapt)

**Инструмент:** MaskPrimers.py (`align --mode cut --maxerror 0.2`) для
multiplex-PCR датасетов; anchored cutadapt для 5′ RACE.

Праймер-сетов: `seq_refs/human_primers.fasta`, `seq_refs/horse_primers.fasta`.
Для овцы праймеры встроены в ноутбук (SMARTer anchor, Sh_IGHG/IGKC/IGLC rev).

QC после primer trim — отдельным шагом:
```python
run_qc("/data/user/epishkin", "PRJNA900592", "pr_trimmed")
```

Результат: `results/<DS>/pr_trimmed/fastq/` (.pr.fastq.gz) +
`results/<DS>/pr_trimmed/fastqc/` + `results/<DS>/pr_trimmed/multiqc/`.

### 5. Сшивание парных ридов (merge)

**Ноутбуки:**
- `notebooks/presto.ipynb` — canonical merge notebook
- `notebooks/presto_mouse.ipynb` — mouse-specific merge notebook для `ERP003950`

**Инструмент:** pRESTO AssemblePairs.py align

```python
run_merge_pairs("/data/user/epishkin", "PRJNA900592", label="pr_trimmed")
```

Оптимизированные параметры:
- `--nproc 8` (8 CPU достаточно, 16 — оптимум)
- `--gzip-output` (сжатый выход)
- `--coord illumina --rc tail`
- без `--log` (verbose лог >2 GB/sample — bottleneck)
- без `--failed` (failed reads не пишем)
- GPU не нужен (pRESTO — CPU-only Python)

Результат: `results/<DS>/merged/fastq/` (_assemble-pass.fastq.gz).

### 6. QC после сшивания

```python
run_qc("/data/user/epishkin", "PRJNA900592", "merged")
```

## Структура репозитория

```
rima/
├── README.md                    этот файл
├── .gitignore                   FASTQ, OneQ-state, checkpoints, .pyc — не коммитить
├── seq_refs/                    референсные последовательности (адаптеры, праймеры)
│   ├── adapter_sequences.fasta  адаптеры (Illumina, horse run-specific)
│   ├── human_primers.fasta      36 human V-gene FR1 + JH rev (MaskPrimers)
│   ├── horse_primers.fasta      35 EquPD v2020 primers (MaskPrimers)
│   └── sheep_5prime_race.fasta  SMARTer anchor + C-region rev (cutadapt)
├── notebooks/                   Jupyter-ноутбуки — основной и canonical рабочий формат
│   ├── qc.ipynb                 QC на любой этап (raw/trimmed/pr_trimmed/merged)
│   ├── adapter_trim_human.ipynb  человек
│   ├── adapter_trim_mouse.ipynb  мышь
│   ├── adapter_trim_horse.ipynb  лошадь
│   ├── adapter_trim_sheep.ipynb  овца
│   ├── primer_trim_human.ipynb   человек
│   ├── primer_trim_mouse.ipynb   мышь
│   ├── primer_trim.ipynb         legacy human + horse
│   ├── primer_trim_sheep.ipynb   овца (5′ RACE cutadapt)
│   ├── presto.ipynb              canonical PE merge
│   └── presto_mouse.ipynb        mouse merge (`ERP003950`)
└── results/                     QC-отчёты (FastQC HTML/ZIP + MultiQC)
    ├── PRJEB40348/              человек
    ├── PRJNA848968/             лошадь
    └── PRJNA900592/             овца
```

Для каждого датасета в `results/<DS>/`:
```
results/<DS>/
├── raw/            QC сырых ридов
│   ├── fastqc/     per-sample FastQC HTML + ZIP
│   └── multiqc/    MultiQC-сводка
├── trimmed/        QC после adapter trim
│   ├── fastp_reports/  fastp JSON + HTML
│   ├── fastqc/     per-sample FastQC
│   └── multiqc/    MultiQC-сводка
├── pr_trimmed/     QC после primer trim
│   ├── fastqc/     per-sample FastQC
│   └── multiqc/    MultiQC-сводка
└── merged/         (на OneQ) результат AssemblePairs
    └── fastq/      _assemble-pass.fastq.gz
```

FASTQ-файлы (сырые, trimmed, pr_trimmed, pr_trimmed_sync, merged) хранятся **только
на OneQ volume** (`/data/user/epishkin/results/<DS>/`) и не
коммитятся в git. В репозитории — notebook'и, QC-отчёты (HTML, ZIP, JSON, MultiQC),
`seq_refs/` и стабильные metadata/summary-артефакты.

## Инструменты

| Инструмент | Версия | Назначение |
|------------|--------|------------|
| FastQC | 0.12.1 | Per-sample QC |
| MultiQC | 1.35 | Агрегация FastQC-отчётов |
| fastp | 0.23+ | Quality filtering, adapter detection (safety net) |
| cutadapt | 5.2 | Adapter trimming, 5′ RACE primer removal |
| pRESTO | 0.7.9 | AssemblePairs (merge PE), MaskPrimers (primer trim) |

Все инструменты устанавливаются в conda-окружение `bcr_env` через
`scripts/setup_vm_conda.sh`.

## Запуск на OneQ

1. **Создать task** в OneQ с образом `jupyter/base-notebook` и
   mounted volume `/data/user/epishkin/`.

2. **Установить окружение** (один раз):
   ```bash
   bash scripts/setup_vm_conda.sh
   ```

3. **Скачать данные** из ENA в `raw/<DS>/`:
   ```bash
   python3 -c "import urllib.request; \
   url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA900592&result=read_run&fields=run_accession,fastq_ftp&format=tsv'; \
   print(urllib.request.urlopen(url).read().decode())"
   # далее wget каждого fastq_ftp URL
   ```

4. **Запустить пайплайн** в Jupyter-ноутбуках (по одному шагу за раз):
   - `qc.ipynb` → `adapter_trim*.ipynb` → `qc.ipynb` →
     `primer_trim*.ipynb` → `qc.ipynb` → `presto.ipynb` → `qc.ipynb`

5. **Скачать QC-отчёты** локально и закоммитить в git.

## Особенности по датасетам

### Человек (PRJEB40348)
- Multiplex PCR с FR1-праймерами (Cheng 2011 universal set: 15 VH fwd + 4 JH rev)
- Адаптеры: Illumina TruSeq/NEBNext (`seq_refs/adapter_sequences.fasta`)
- Primer trim: MaskPrimers `align --mode cut` с `seq_refs/human_primers.fasta`
- Retention после adapter+primer trim: 88.1% (77–95% per sample)

### Лошадь (PRJNA848968)
- EquPD v2020 primer set, phage display scFv libraries
- Адаптеры: эмпирические run-specific R1 (`seq_refs/adapter_sequences.fasta`),
  НЕ копировать human TruSeq; EquPD/scFv tails НЕ резать на adapter-этапе
- Primer trim: MaskPrimers `align` с `seq_refs/horse_primers.fasta` (35 праймеров)

### Овца (PRJNA900592)
- 5′ RACE (SMARTer), универсальный anchor (не multiplex V-primer PCR)
- Адаптеры: Illumina/NEBNext read-through (`seq_refs/adapter_sequences.fasta`),
  без fastp autodetect (детектит SMARTer anchor как "адаптер")
- Primer trim: anchored cutadapt — SMARTer anchor (5′) +
  C-region reverse primers (`seq_refs/sheep_5prime_race.fasta`, 3′)
- Retention: raw → adapter trimmed 78.2% → primer trimmed ~100% (только срез тех. seqs)
- Read length: ~300 bp → ~266 bp после primer trim

## Источники

- **Лошадь:** Rosenfeld R et al. (2022) Centaur antibodies. *Front. Immunol.* 13:942317. doi:10.3389/fimmu.2022.942317. PRJNA848968.
- **Овца:** Park M et al. (2023) Exploring the sheep immunoglobulin repertoire. *Mol. Immunol.* 156:20–30. doi:10.1016/j.molimm.2023.02.008. PRJNA900592.
- **Человек:** Lomakin YA et al. (2022) Deconvolution of B cell receptor repertoire. *Front. Immunol.* 13:803229. doi:10.3389/fimmu.2022.803229. PRJEB40348 / E-MTAB-9573.
- **pRESTO:** Rosenfeld et al. horse paper: pRESTO v0.7.0, Phred Q10 filtering. *Front. Immunol.* 2022.
