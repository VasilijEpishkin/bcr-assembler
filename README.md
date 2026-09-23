# bcr-assembler

Секвенируемая V(D)J-последовательность обычно длиннее, чем покрывает пара
150bp ридов, поэтому библиотеку фрагментируют перед секвенированием, и каждый
фрагмент читается отдельно с двух концов. Задача проекта — проверить, можно ли
после такой фрагментации корректно сшить исходную последовательность обратно,
и какие инструменты (`TRUST4` и аналоги) делают это надёжно. Для честной
проверки нужен known ground truth, поэтому пайплайн сначала строит
reference-like baseline из реальных данных, а затем генерирует на его основе
синтетические фрагментированные reads с контролируемой правдой.

Актуальные simulation-ветки строятся после аннотации и post-filtering, а не
напрямую из старых merged-read simulation notebooks.

## Актуальные simulation notebooks

### Human PRJEB30386

- `results/PRJEB30386/notebooks/simulate_human_insilicoseq_150bp_novaseq.ipynb` — random-cut baseline
- `results/PRJEB30386/notebooks/simulate_human_insilicoseq_150bp_custom_umi_consensus.ipynb` — random-cut baseline
- `results/PRJEB30386/notebooks/simulate_human_insilicoseq_150bp_novaseq_ultrasonic_fragmentation.ipynb` — ultrasonic-like
- `results/PRJEB30386/notebooks/simulate_human_insilicoseq_150bp_custom_umi_consensus_ultrasonic_fragmentation.ipynb` — ultrasonic-like

### Mouse ERP003950 fastp

- `results/ERP003950_fastp_q30_u40/notebooks/simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` — random-cut baseline
- `results/ERP003950_fastp_q30_u40/notebooks/simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq_ultrasonic_fragmentation.ipynb` — ultrasonic-like

Старые `simulate_*_merged_*.ipynb` simulation notebooks удалены как legacy /
дубликаты и не являются частью production workflow.

## Пайплайн mouse fastp

| Шаг | Ноутбук | Инструмент |
|---|---|---|
| QC (на всех стадиях) | `qc.ipynb` | FastQC, MultiQC |
| Adapter trim | `adapter_trim_mouse_fastp_q30_u40.ipynb` | cutadapt, fastp |
| Primer trim | `primer_trim_mouse_fastp_q30_u40.ipynb` | cutadapt |
| Merge paired-end reads | `presto_mouse_fastp_q30_u40.ipynb` | pRESTO `AssemblePairs.py` |
| Annotation | `annotate_mouse_fastp_q30_u40.ipynb` | IgBLAST |
| Post-annotation filter | `filter_mouse_post_annotation.ipynb` | AIRR/sequence filters |
| Симуляция секвенирования | `simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` | InSilicoSeq |
| Валидация симуляции | `validate_mouse_novaseq_post_annotation_filtered.ipynb` | bowtie2, samtools |
| Сшивка фрагментов | — | `TRUST4` и аналоги |

## Симуляция: общая архитектура

Актуальная архитектура:

```text
post_annotation_filtered truth
    ↓
PCR1
    ↓
fragmentation
    ↓
PCR2
    ↓
read allocation
    ↓
InSilicoSeq
    ↓
PE150 FASTQ
```

Фрагментация в production notebooks использует random-cut baseline: одна случайная межнуклеотидная точка разрыва на выбранную library-input молекулу, оба дочерних фрагмента создаются до size-selection. Отдельные ultrasonic-like notebooks используют рекурсивные size-dependent разрывы с midpoint-centered breakpoint distribution и отдельным size-selection.

## Окружение

```bash
source scripts/setup_env.sh
```

Собирает `bcr_env` (fastqc, fastp, cutadapt, multiqc, presto, bowtie2,
samtools, insilicoseq, rsync) и регистрирует Jupyter-ядро "BCR Pipeline".

## Структура

- `notebooks/` — общие notebooks;
- `results/PRJEB30386/` — human pipeline и simulation;
- `results/ERP003950_fastp_q30_u40/` — актуальная mouse fastp production branch;
- `scripts/setup_env.sh` — окружение.
