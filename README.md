# bcr-assembler

Секвенируемая V(D)J-последовательность обычно длиннее, чем покрывает пара
150bp ридов, поэтому библиотеку фрагментируют перед секвенированием, и каждый
фрагмент читается отдельно с двух концов. Задача проекта — проверить, можно ли
после такой фрагментации корректно сшить исходную последовательность обратно,
и какие инструменты (`TRUST4` и аналоги) делают это надёжно. Для честной
проверки нужен known ground truth, поэтому пайплайн сначала строит
reference-like baseline из реальных данных, а затем генерирует на его основе
синтетические фрагментированные reads с контролируемой правдой.

Рабочий кейс — mouse `ERP003950` (Greiff 2014, IgG heavy-chain репертуар,
6 сэмплов, MiSeq 2×250). Остальные виды в репозитории — эксплораторный
материал, не описаны здесь.

## Пайплайн

| Шаг | Ноутбук | Инструмент |
|---|---|---|
| QC (на всех стадиях) | `qc.ipynb` | FastQC, MultiQC |
| Adapter trim | `adapter_trim_mouse.ipynb` | cutadapt |
| Primer trim (только constant-region; V-region primer сохраняется) | `primer_trim_mouse.ipynb` | cutadapt |
| Merge paired-end reads | `presto_mouse.ipynb` | pRESTO `AssemblePairs.py` |
| Annotation | `annotate_mouse.ipynb` | IgBLAST |
| Quality summary | `mouse_full_quality_summary_6samples.ipynb` | — |
| Сравнение аннотаторов | `annotator_compare_mouse.ipynb`, `igblast_abstar_problematic_2seq.ipynb` | IgBLAST vs abstar |
| Симуляция секвенирования | `simulate_mouse_merged_insilicoseq_150bp.ipynb` | InSilicoSeq, bowtie2 |
| Валидация симуляции | `align_simulated_reads_bowtie2.ipynb` | bowtie2, samtools |
| Сшивка фрагментов в исходную последовательность | — | `TRUST4` и аналоги |

## Симуляция: как устроена

`simulate_mouse_merged_insilicoseq_150bp.ipynb` — dedup merged reads в
уникальные templates → error model (KDE, обучена на реальных ридах, обрезанных
до 150bp) → **PCR#1** (branching, на целых templates, до фрагментации) →
**фрагментация** (templates явно нарезаются на отслеживаемые кандидат-фрагменты,
не эфемерно внутри InSilicoSeq) → **PCR#2** (независимый branching на каждом
фрагменте, library-prep после фрагментации) → мультиномиальная аллокация read
budget по фрагментам → `iss generate --sequence_type amplicon` (сиквенирует
фрагмент как есть, без повторной нарезки).

Известное упрощение: PCR#2 стартует с полного `pcr_copies` родителя на каждый
кандидат-фрагмент, а не тайлит одну амплифицированную молекулу физически.

## Окружение

```bash
source scripts/setup_env.sh
```

Собирает `bcr_env` (fastqc, fastp, cutadapt, multiqc, presto, bowtie2,
samtools, insilicoseq, rsync), регистрирует Jupyter-ядро "BCR Pipeline". Только
в OneQ Jupyter Terminal — через SSH более тесный memory cgroup убивает тяжёлые
conda-установки.

## Структура

- `notebooks/` — пайплайн (mouse core + эксплораторные ноутбуки по другим видам)
- `results/` — QC/аннотации/summary; большие бинарники (FASTQ, IgBLAST TSV,
  симуляция) в git не хранятся, живут только на OneQ
- `scripts/setup_env.sh` — окружение
