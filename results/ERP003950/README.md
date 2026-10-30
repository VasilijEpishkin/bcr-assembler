# ERP003950 — мышь, IGH

Ноутбуки и отчёты QC датасета ERP003950 (тяжёлая цепь IgG, шесть образцов ERR346596–ERR346601).

## Порядок обработки

1. удаление адаптеров cutadapt и фильтрация пар fastp (`-q 30 -u 40 -l 250`) — `adapter_trim_mouse.ipynb`;
2. удаление праймеров — `primer_trim_mouse.ipynb`;
3. сборка пар ридов pRESTO — `presto_mouse.ipynb`;
4. аннотация IgBLAST и отчёт — `annotate_mouse.ipynb`, `annotation_report_mouse.ipynb`;
5. фильтрация после аннотации — `filter_mouse_post_annotation.ipynb`;
6. симуляция NovaSeq PE150 (основная линия: ПЦР1 30 циклов, порция 1/10 продукта, случайный разрез, отбор 200–350 нт) —
   `simulate_mouse_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb`, валидация —
   `validate_mouse_novaseq_post_annotation_filtered.ipynb`;
7. сборка (TRUST4, rnaSPAdes, Trinity) и сравнение с эталоном симуляции — `assemble_*_mouse.ipynb`,
   `benchmark_assemblers_mouse.ipynb`.

Другие сборщики:
`assemble_mixcr_mouse.ipynb` (MiXCR), `assemble_gather_mouse.ipynb` (GATHeR). Распределения длин ридов по
стадиям — `read_geometry_mouse.ipynb`.

## Содержимое

- `notebooks/` — исполняемые ноутбуки;
- `<стадия>/fastqc/`, `<стадия>/multiqc/` — отчёты FastQC и MultiQC (`trimmed`, `pr_trimmed`, `merged`);
- `merged/assembly_qc.tsv` — сводка сборки пар ридов.
