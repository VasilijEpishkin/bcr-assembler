# PRJNA848968 — лошадь

Предобработка датасета PRJNA848968: удаление адаптеров cutadapt и фильтрация пар fastp (`-q 30 -u 40 -l 250`) —
`adapter_trim_horse.ipynb`, удаление праймеров — `primer_trim_horse.ipynb`. Вне основной линии сборки.

- `notebooks/` — ноутбуки;
- `<стадия>/fastqc/`, `<стадия>/multiqc/` — отчёты FastQC и MultiQC (`raw`, `trimmed`);
- `trimmed/fastp_reports/` — отчёты fastp.

FASTQ лежат на томе данных и в git не хранятся.
