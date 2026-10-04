# PRJNA900592 — овца

Предобработка датасета PRJNA900592: удаление адаптеров cutadapt и фильтрация пар fastp (`-q 30 -u 40 -l 250`) —
`adapter_trim_sheep.ipynb`, удаление праймеров — `primer_trim_sheep.ipynb`, сборка пар ридов pRESTO —
`presto_sheep.ipynb`, сравнение аннотаторов IgBLAST и AbStar — `annotator_compare_sheep.ipynb`.

- `notebooks/` — ноутбуки;
- `<стадия>/fastqc/`, `<стадия>/multiqc/` — отчёты FastQC и MultiQC (`raw`, `trimmed`, `pr_trimmed`, `merged`); `merged/assembly_qc.tsv`, `merged/sheep_merge_dashboard.html` — сводки сборки пар ридов;
- `annotator_compare/` — результаты и отчёты сравнения аннотаторов.
