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

Симуляция: ПЦР1 → порция **5 нг** → один случайный разрез → отбор по размеру 250±40 → ПЦР2 → NovaSeq PE150.
Ноутбуки основной линии (ветки по умолчанию — `…_random_cut_amp_{umicons,ec1x3}_in5ng_rb3x_ss250`):

| Шаг | Человек `results/PRJEB30386/notebooks/` | Мышь `results/ERP003950/notebooks/` |
|---|---|---|
| Симуляция | `simulate_human_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` | `simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb` |
| Валидация симуляции | `validate_human_novaseq_post_annotation_filtered.ipynb` | `validate_mouse_novaseq_post_annotation_filtered.ipynb` |
| Сборка | `assemble_{trust4,rnaspades,trinity}_human.ipynb` | `assemble_{trust4,rnaspades,trinity}_mouse.ipynb` |
| Оценка | `benchmark_assemblers_human.ipynb` | `benchmark_assemblers_mouse.ipynb` |

Результаты и описание метрик — `docs/evaluation_history/run2_base_report.md`.

Другие варианты симуляции: `simulate_*_pcr1_tails.ipynb` — праймеры ПЦР1 с хвостами (адаптеры, UMI) и
разбросом эффективности ПЦР1; `simulate_*_ultrasonic_fragmentation.ipynb` — ультразвуковая фрагментация
вместо одного случайного разреза. Другие сборщики: `assemble_mixcr_*.ipynb` (MiXCR), `assemble_vdjer_*.ipynb` (V'DJer).

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
`benchmark_assemblers_{human,mouse}.ipynb` (`<ветка>/benchmark_unique_trim{20,0}/`), отчёт по стадиям
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

Основная фрагментация — один случайный разрез между нуклеотидами каждой молекулы порции; оба дочерних фрагмента
проходят отбор по размеру. Вариант с ультразвуковой фрагментацией — повторные разрывы с вероятностью, зависящей
от длины, точка разрыва распределена вокруг середины фрагмента, затем отдельный отбор по размеру.

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
- `docs/` — отчёты, история оценок, раскладка ноутбуков;
- `scripts/` — вспомогательные скрипты (`setup_env.sh` — окружение, `igbrowser/` — R-отчёты igbrowser).
