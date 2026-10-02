# PRJEB40348 — человек

Предобработка датасета PRJEB40348: удаление адаптеров cutadapt и фильтрация пар fastp (`-q 30 -u 40 -l 250`) —
`adapter_trim_human.ipynb`, удаление праймеров — `primer_trim_human.ipynb`. Вне основной линии сборки.

- `notebooks/` — ноутбуки;
- `<стадия>/fastqc/`, `<стадия>/multiqc/` — отчёты FastQC и MultiQC (`raw`, `trimmed`, `pr_trimmed`);
- `pr_trimmed/fastp_reports/` — отчёты fastp и маскирования праймеров.

FASTQ лежат на томе данных и в git не хранятся.
