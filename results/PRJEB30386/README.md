# PRJEB30386 — человек, IGH/IGK/IGL

Ноутбуки и отчёты QC датасета PRJEB30386 (библиотеки 5′MTPX с UMI: IgM, IgG, IgK, IgL).

## Порядок обработки

1. удаление адаптеров cutadapt и фильтрация пар fastp (`-q 30 -u 40 -l 250`), без обрезки концов по качеству — `adapter_trim_human.ipynb`;
2. удаление праймеров — `primer_trim_human.ipynb`;
3. сборка пар ридов pRESTO — `presto_human.ipynb`;
4. аннотация IgBLAST и отчёт — `annotate_human.ipynb`, `annotation_report_human.ipynb`;
5. фильтрация после аннотации — `filter_human_post_annotation.ipynb`;
6. симуляция NovaSeq PE150 (основная линия: порция 5 нг продукта ПЦР1, случайный разрез) —
   `simulate_human_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb`, валидация —
   `validate_human_novaseq_post_annotation_filtered.ipynb`;
7. сборка (TRUST4, rnaSPAdes, Trinity) и сравнение с эталоном симуляции — `assemble_*_human.ipynb`,
   `benchmark_assemblers_human.ipynb`.

Варианты симуляции, ещё не проверенные как основная линия: `*_pcr1_tails.ipynb`, `*_ultrasonic_fragmentation.ipynb`;
сборщики вне основной линии: `assemble_mixcr_human.ipynb`, `assemble_vdjer_human.ipynb`.
Диагностика: `read_geometry_human.ipynb` (длины ридов по стадиям), `align_simulated_reads_bowtie2_human.ipynb`,
`audit_human_v_primer_trimming.ipynb`.

## Содержимое

- `notebooks/` — исполняемые ноутбуки;
- `<стадия>/fastqc/`, `<стадия>/multiqc/` — отчёты FastQC и MultiQC (`raw`, `trimmed`, `pr_trimmed`, `merged`);
- `trimmed/fastp_reports/`, `*/filter_summary.json`, `merged/qc/assembly_qc.tsv` — сводки стадий;
- `post_annotation_filtered/qc/` — QC последовательностей после фильтрации;
- `qc/read_geometry/` — распределения длин ридов.

FASTQ, полные AIRR-таблицы и результаты симуляции лежат на томе данных и в git не хранятся.
