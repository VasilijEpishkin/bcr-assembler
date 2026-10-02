# ERP003950 — мышь, IGH

Ноутбуки и отчёты QC датасета ERP003950 (тяжёлая цепь IgG, шесть образцов ERR346596–ERR346601).

## Порядок обработки

1. удаление адаптеров cutadapt и фильтрация пар fastp (`-q 30 -u 40 -l 250`), без обрезки концов по качеству — `adapter_trim_mouse.ipynb`;
2. удаление праймеров — `primer_trim_mouse.ipynb`;
3. сборка пар ридов pRESTO — `presto_mouse.ipynb`;
4. аннотация IgBLAST и отчёт — `annotate_mouse.ipynb`, `annotation_report_mouse.ipynb`;
5. фильтрация после аннотации — `filter_mouse_post_annotation.ipynb`;
6. симуляция NovaSeq PE150 (основная линия: порция 5 нг продукта ПЦР1, случайный разрез) —
   `simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb`, валидация —
   `validate_mouse_novaseq_post_annotation_filtered.ipynb`;
7. сборка (TRUST4, rnaSPAdes, Trinity) и сравнение с эталоном симуляции — `assemble_*_mouse.ipynb`,
   `benchmark_assemblers_mouse.ipynb`.

Варианты симуляции, ещё не проверенные как основная линия: `*_pcr1_tails.ipynb`, `*_ultrasonic_fragmentation.ipynb`;
сборщики вне основной линии: `assemble_mixcr_mouse.ipynb`, `assemble_vdjer_mouse.ipynb`.
Диагностика: `read_geometry_mouse.ipynb` (длины ридов по стадиям).

## Содержимое

- `notebooks/` — исполняемые ноутбуки;
- `<стадия>/qc/fastqc/`, `<стадия>/qc/multiqc/` — отчёты FastQC и MultiQC (`trimmed`, `pr_trimmed`, `merged`);
- `merged/qc/assembly_qc.tsv` — сводка сборки пар ридов.

FASTQ, полные AIRR-таблицы, индексы SQLite и результаты симуляции лежат на томе данных и в git не хранятся.
