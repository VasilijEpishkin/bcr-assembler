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

### Mouse ERP003950

- `results/ERP003950/notebooks/simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` — random-cut baseline
- `results/ERP003950/notebooks/simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq_ultrasonic_fragmentation.ipynb` — ultrasonic-like

Старые `simulate_*_merged_*.ipynb` simulation notebooks удалены как legacy /
дубликаты и не являются частью production workflow.

## Пайплайн mouse ERP003950 (fastp preprocessing)

| Шаг | Ноутбук | Инструмент |
|---|---|---|
| QC (на всех стадиях) | `qc.ipynb` | FastQC, MultiQC |
| Adapter trim | `adapter_trim_mouse.ipynb` | cutadapt, fastp |
| Primer trim | `primer_trim_mouse.ipynb` | cutadapt |
| Merge paired-end reads | `presto_mouse.ipynb` | pRESTO `AssemblePairs.py` |
| Annotation | `annotate_mouse.ipynb` | IgBLAST |
| Post-annotation filter | `filter_mouse_post_annotation.ipynb` | AIRR/sequence filters |
| Симуляция секвенирования | `simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` | InSilicoSeq |
| Валидация симуляции | `validate_mouse_novaseq_post_annotation_filtered.ipynb` | bowtie2, samtools |
| Сшивка фрагментов | — | `TRUST4` и аналоги |

## Реконструкция TRUST4

Единый runner для всех актуальных human/mouse simulation branches:

- `notebooks/run_trust4_simulated_bcr_universal.ipynb`

Он запускает TRUST4 на PE150 FASTQ, сохраняет нативные TRUST4 outputs и создаёт
нормализованные `reconstructed_contigs.fasta/tsv`, `cdr3_normalized.tsv` и
`trust4_run_summary.tsv` для последующего truth benchmark.

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

Эталон — физический ампликон PCR1 (склеенный рид с восстановленными праймерами); интервал V…J записан в `vj_start`/`vj_end` и служит эталоном для оценки сборки. У человека каждая UMI-молекула считается один раз (`umi_mode`). Порция на фрагментацию (`LIBRARY_INPUT_SCALE`), множитель глубины (`READ_BUDGET_MULTIPLIER`) и центр отбора по длине (`SIZE_SELECTION_TARGET`, по умолчанию 250 нт) входят в имя ветки.

Фрагментация в production notebooks использует random-cut baseline: одна случайная межнуклеотидная точка разрыва на выбранную library-input молекулу, оба дочерних фрагмента создаются до size-selection. QC `12a` классифицирует каждый шаблон по восстановимости V…J (перекрытие ридов / только mate-пара / физически не связан / частичное покрытие / нет ридов). Отдельные ultrasonic-like notebooks используют рекурсивные size-dependent разрывы с midpoint-centered breakpoint distribution и отдельным size-selection.

## Окружение

```bash
source scripts/setup_env.sh
```

Собирает `bcr_env` (fastqc, fastp, cutadapt, multiqc, presto, bowtie2,
samtools, insilicoseq, rsync) и регистрирует Jupyter-ядро "BCR Pipeline".

## Структура

- `notebooks/` — общие notebooks;
- `results/PRJEB30386/` — human pipeline и simulation;
- `results/ERP003950/` — mouse pipeline (fastp preprocessing) и simulation;
- `scripts/setup_env.sh` — окружение.
