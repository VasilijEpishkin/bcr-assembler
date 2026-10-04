# PRJEB30386 — человек, IGH/IGK/IGL

Ноутбуки и отчёты QC датасета PRJEB30386 (библиотеки 5′MTPX с UMI: IgM, IgG, IgK, IgL).

## Порядок обработки

1. удаление адаптеров cutadapt и фильтрация пар fastp (`-q 30 -u 40 -l 250`) — `adapter_trim_human.ipynb`;
2. удаление праймеров — `primer_trim_human.ipynb`;
3. сборка пар ридов pRESTO — `presto_human.ipynb`;
4. аннотация IgBLAST и отчёт — `annotate_human.ipynb`, `annotation_report_human.ipynb`;
5. фильтрация после аннотации — `filter_human_post_annotation.ipynb`;
6. симуляция NovaSeq PE150 (основная линия: порция 5 нг продукта ПЦР1, случайный разрез) —
   `simulate_human_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb`, валидация —
   `validate_human_novaseq_post_annotation_filtered.ipynb`;
7. сборка (TRUST4, rnaSPAdes, Trinity) и сравнение с эталоном симуляции — `assemble_*_human.ipynb`,
   `benchmark_assemblers_human.ipynb`.

Другие варианты симуляции: `*_pcr1_tails.ipynb` — праймеры ПЦР1 с хвостами (адаптеры, UMI) и разбросом
эффективности ПЦР1; `*_ultrasonic_fragmentation.ipynb` — ультразвуковая фрагментация. Другие сборщики:
`assemble_mixcr_human.ipynb` (MiXCR), `assemble_vdjer_human.ipynb` (V'DJer). Распределения длин ридов по стадиям —
`read_geometry_human.ipynb`; выравнивание симулированных ридов на шаблоны bowtie2 —
`align_simulated_reads_bowtie2_human.ipynb`; проверка обрезки V-праймеров — `audit_human_v_primer_trimming.ipynb`.

## Содержимое

- `notebooks/` — исполняемые ноутбуки;
- `<стадия>/fastqc/`, `<стадия>/multiqc/` — отчёты FastQC и MultiQC (`raw`, `trimmed`, `pr_trimmed`, `merged`);
- `trimmed/fastp_reports/`, `*/filter_summary.json`, `merged/assembly_qc.tsv` — сводки стадий;
- `post_annotation_filtered/fastqc/`, `post_annotation_filtered/multiqc/` — QC последовательностей после фильтрации;
- `read_geometry/` — распределения длин ридов.
