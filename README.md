# bcr-assembler

Секвенируемая V(D)J-последовательность обычно длиннее, чем покрывает пара
150bp ридов, поэтому библиотеку фрагментируют перед секвенированием, и каждый
фрагмент читается отдельно с двух концов. Задача проекта — проверить, можно ли
после такой фрагментации корректно сшить исходную последовательность обратно,
и какие инструменты (`TRUST4` и аналоги) делают это надёжно. Для честной
проверки нужен known ground truth, поэтому пайплайн сначала строит
reference-like baseline из реальных данных, а затем генерирует на его основе
синтетические фрагментированные reads с контролируемой правдой.

Симуляция строится по последовательностям после аннотации и фильтрации.

## Основная линия

Симуляция: ПЦР1 (30 циклов, разброс эффективности) → случайная порция 1/10 продукта → один случайный разрез → мягкий отбор по длине
(диапазон 200–350 нт, вне него вероятность падает вдвое каждые 25 нт, фрагменты короче 152 нт не берутся) → ПЦР2 (5 циклов) →
1/1000 продукта на секвенирование → NovaSeq PE150. Ноутбуки основной линии (ветки по умолчанию —
`…_random_cut_amp_{umicons,ec1x3}_rb3x_ss200-350h25_pcr1c30_pcr2c5_seq1e-3_al1e-1`; вариант с входной массой 200 нг — `…_in200ng_rb3x_…`):

| Шаг | Человек `results/PRJEB30386/notebooks/` | Мышь `results/ERP003950/notebooks/` |
|---|---|---|
| Симуляция | `simulate_human_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` | `simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` |
| Валидация симуляции | `validate_human_novaseq_post_annotation_filtered.ipynb` | `validate_mouse_novaseq_post_annotation_filtered.ipynb` |
| Сборка | `assemble_{trust4,mixcr,gather,rnaspades,trinity}_human.ipynb` | `assemble_{trust4,mixcr,gather,rnaspades,trinity}_mouse.ipynb` |
| Оценка | `benchmark_assemblers_human.ipynb` | `benchmark_assemblers_mouse.ipynb` |

Результаты и описание метрик — `docs/evaluation_history/assembler_comparison_report.md`; параметры сценариев симуляции и их влияние на сборку — `docs/simulation_scenarios.md`.

Режимы SPAdes (`rna`, `meta`, `rnaviral`) задаются `BCR_SPADES_MODE` в `assemble_rnaspades_*.ipynb`; ещё один сборщик — `assemble_vdjer_*.ipynb` (V'DJer).

## Подготовка данных мыши ERP003950

`results/ERP003950/` — ветка с фильтрацией fastp (`-q 30 -u 40 -l 250`) под каноническим именем датасета.

| Шаг | Ноутбук | Инструмент |
|---|---|---|
| QC (на всех стадиях) | `notebooks/qc.ipynb` | FastQC, MultiQC |
| Adapter trim | `adapter_trim_mouse.ipynb` | cutadapt, fastp |
| Primer trim | `primer_trim_mouse.ipynb` | cutadapt |
| Merge paired-end reads | `presto_mouse.ipynb` | pRESTO `AssemblePairs.py` |
| Annotation | `annotate_mouse.ipynb`, `annotation_report_mouse.ipynb` | IgBLAST |
| Post-annotation filter | `filter_mouse_post_annotation.ipynb` | AIRR/sequence filters |

## Сборка и оценка сборщиков

Сборщики запускаются на PE150 FASTQ симуляций отдельными ноутбуками (ветка симуляции
задаётся `BCR_BRANCH`, число потоков — `BCR_THREADS`). Нативные результаты и нормализованный
`contigs.fasta` пишутся в `<ветка>/assemblers/<сборщик>/<образец>/`. Сравнение с эталоном —
`benchmark_assemblers_{human,mouse}.ipynb` (`<ветка>/benchmark_unique_trim{20,0,45_15}/`), отчёт по стадиям
симуляции — `notebooks/simulation_report.ipynb`, история всех оценок — `docs/evaluation_history/`.

## Симуляция: общая схема

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

Фрагментация — один случайный разрез между нуклеотидами каждой молекулы порции; оба дочерних фрагмента
проходят отбор по размеру.

## Окружение

```bash
source scripts/setup_env.sh
```

Собирает `bcr_env` (fastqc, fastp, cutadapt, multiqc, presto, bowtie2,
samtools, insilicoseq, rsync) и регистрирует Jupyter-ядро "BCR Pipeline".

## Структура

- `notebooks/` — общие ноутбуки (`qc.ipynb` — FastQC и MultiQC любой стадии, `simulation_report.ipynb` — отчёт по стадиям симуляции);
- `results/<датасет>/notebooks/` — ноутбуки стадий датасета; `results/<датасет>/<стадия>/{fastqc,multiqc}/` — отчёты QC стадии;
- `results/PRJEB30386/` — человек: подготовка, симуляция, сборка и оценка;
- `results/ERP003950/` — мышь: подготовка (fastp), симуляция, сборка и оценка;
- `results/PRJEB40348/`, `results/PRJNA848968/`, `results/PRJNA900592/`, `results/PRJNA1226555/` —
  предобработка других датасетов (человек, лошадь, овца, лёгкие цепи мыши);
- `docs/` — отчёты, история оценок, сценарии симуляции, раскладка ноутбуков;
- `scripts/` — вспомогательные скрипты (`setup_env.sh` — окружение, `collect_evaluations.py` — сводка оценок в `docs/evaluation_history/`, `igbrowser/` — R-отчёты igbrowser).
