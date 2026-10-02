# ERP003950 mouse workflow

Canonical notebook-first workflow for the ERP003950 mouse IGH dataset.

## Processing sequence

1. explicit adapter removal with cutadapt;
2. whole-read filtering with fastp (`-q 30 -u 40 -l 250`), without quality-end trimming;
3. constant-primer removal;
4. paired-read assembly;
5. mouse IgBLAST AIRR annotation;
6. post-annotation filtering;
7. NovaSeq PE150 simulation (main line: 5 ng PCR1 aliquot, random cut), validation and QC;
8. assembly (TRUST4, rnaSPAdes, Trinity) and benchmark against the simulation truth
   (`assemble_*_mouse.ipynb`, `benchmark_assemblers_mouse.ipynb`).

Simulation variants not yet validated as the main line: `*_pcr1_tails.ipynb`, `*_ultrasonic_fragmentation.ipynb`;
assembler notebooks outside the main line: `assemble_mixcr_mouse.ipynb`, `assemble_vdjer_mouse.ipynb`.

## Repository contents

- `notebooks/` — executable workflow notebooks;
- `<stage>/qc/fastqc/` — FastQC reports;
- `<stage>/qc/multiqc/` — MultiQC report and supporting data;
- lightweight manifests and summary tables.

FASTQ, full AIRR TSV, SQLite indexes and other large runtime artifacts remain on the data volume and are not stored in Git.
