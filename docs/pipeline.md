# Конвейер: от сырых ридов до оценки сборщиков

Для каждого ноутбука: что он берёт на вход, что делает с данными и что кладёт на выход. Пути даны от корня проекта
(`bcr-assembler/`); `<DS>` — датасет (`ERP003950` мышь, `PRJEB30386` человек), `<S>` — образец (run),
`<ветка>` — ветка симуляции (`results/<DS>/simulated/<ветка>/`). Ноутбуки датасета лежат в `results/<DS>/notebooks/`,
общие — в `notebooks/`. Подробности о параметрах — в `docs/report/project_report.md` (разделы 2–13).

## Схема целиком

```text
raw/<DS>/<S>_{1,2}.fastq.gz                                 сырые пары MiSeq (мышь 2×250, человек 2×300)
  │ adapter_trim_*        cutadapt (адаптеры) → fastp (фильтр пар)
  ▼
results/<DS>/trimmed/fastq/<S>_{1,2}.trim.fastq.gz
  │ primer_trim_*         мышь: cutadapt (C-праймер); человек: MaskPrimers (V- и C-праймеры, UMI → заголовок) + PairSeq
  ▼
results/<DS>/pr_trimmed/fastq/<S>_{1,2}.pr.fastq.gz
  │ presto_*              AssemblePairs: R1 + R2 → один склеенный рид
  ▼
results/<DS>/merged/fastq/<S>_assemble-pass.fastq.gz
  │ annotate_*            IgBLAST → таблица AIRR на каждый склеенный рид
  ▼
results/<DS>/annotation/igblast/<S>.airr.tsv
  │ filter_*_post_annotation   отбор полных продуктивных V…J нужного локуса
  ▼
results/<DS>/post_annotation_filtered/{fastq,airr_pass}/
  │ simulate_*            эталон → ПЦР1 → порция → фрагментация → ПЦР2 → секвенирование → InSilicoSeq
  ▼
results/<DS>/simulated/<ветка>/00_primary_truth/ (эталон)  +  06_fastq_pe150/ (симулированные PE150)
  │ validate_*            проверка симуляции          → <ветка>/validation/
  │ assemble_<сборщик>_*  сборка симулированных ридов → <ветка>/assemblers/<сборщик>/<образец>/contigs.fasta
  ▼
  │ benchmark_assemblers_*   контиги против эталона
  ▼
results/<DS>/simulated/<ветка>/benchmark_unique_trim*/assemblers_summary.tsv
  │ scripts/collect_evaluations.py, scripts/make_report_figures.py
  ▼
docs/evaluation_history/*.tsv, docs/report/figures/*.png
```

После каждой стадии с FASTQ (raw, trimmed, pr_trimmed, merged, post_annotation_filtered) запускается
`notebooks/qc.ipynb`: FastQC и MultiQC кладутся внутрь стадии, в `<стадия>/fastqc/` и `<стадия>/multiqc/`
(для сырых ридов — `results/<DS>/raw/`).

## 1. Тримминг адаптеров и фильтр по качеству — `adapter_trim_{mouse,human}.ipynb`

| | |
|---|---|
| вход | `raw/<DS>/<S>_1.fastq.gz`, `<S>_2.fastq.gz` |
| что происходит | cutadapt вырезает адаптеры Illumina TruSeq с 3′-конца R1 и R2 (там, где рид дочитал до адаптера); fastp удаляет пару целиком, если в одном из ридов больше 40% оснований ниже Q30 или он короче 250 нт (`-q 30 -u 40 -l 250`); сами риды fastp не обрезает |
| выход | `results/<DS>/trimmed/fastq/<S>_{1,2}.trim.fastq.gz` — те же пары без адаптеров, плохие пары отброшены |
| отчёты | `trimmed/logs/` (cutadapt), `trimmed/fastp_reports/<S>.fastp.{json,html}`, `trimmed/filter_summary.json` |

```text
R1: [праймер V]──────── V(D)J ────────[адаптер]      →   [праймер V]──────── V(D)J ────────
пара, где >40% оснований < Q30                        →   удалена
```

## 2. Тримминг праймеров — `primer_trim_{mouse,human}.ipynb`

| | |
|---|---|
| вход | `results/<DS>/trimmed/fastq/<S>_{1,2}.trim.fastq.gz` |
| что происходит, мышь | cutadapt вырезает C-праймер IgG (вырожденный `CARKGGATRRRCHGATGGGG`) в том риде, где он стоит; V-праймер остаётся (он лежит на FR1, внутри V…J); риды не отбрасываются |
| что происходит, человек | pRESTO MaskPrimers ищет и вырезает праймер leader-области V в R1 и праймер константной области в R2, имя праймера пишет в заголовок (`PRIMER=`); из R2 вырезается 20-нт UMI и пишется в заголовок (`BARCODE=`); риды без праймера отбрасываются; PairSeq восстанавливает пары и копирует `BARCODE` на R1 |
| выход | `results/<DS>/pr_trimmed/fastq/<S>_{1,2}.pr.fastq.gz` |
| отчёты | `pr_trimmed/logs/` (человек) |

```text
человек  R2: [UMI 20][C-праймер]── C ── J ── CDR3 …   →   ── C ── J ── CDR3 …   (заголовок: BARCODE=UMI, PRIMER=…)
мышь     R2: [C-праймер]── J ── CDR3 …                →   ── J ── CDR3 …
```

## 3. Склейка пар — `presto_{mouse,human}.ipynb`

| | |
|---|---|
| вход | `results/<DS>/pr_trimmed/fastq/<S>_{1,2}.pr.fastq.gz` |
| что происходит | pRESTO `AssemblePairs.py align`: R2 переводится в обратный комплемент, перекрытие конца R1 с началом R2 ищется выравниванием (8–1000 нт, ≤30% несовпадений, p ≤ 1e−5); пара превращается в один рид — весь ампликон; несклеенные пары отбрасываются; у человека `BARCODE` переносится в склеенный рид |
| выход | `results/<DS>/merged/fastq/<S>_assemble-pass.fastq.gz` (+ `*_assemble-fail.fastq.gz` у человека) |
| отчёты | `merged/assembly_qc.tsv`, `merged/logs/` |

```text
R1  ─────────────────────►
R2               ◄─────────────────────   (revcomp)
    ───────────────────────────────────   склеенный рид = ампликон без вырезанных праймеров
```

## 4. Аннотация — `annotate_{mouse,human}.ipynb`, отчёт — `annotation_report_{mouse,human}.ipynb`

| | |
|---|---|
| вход | `results/<DS>/merged/fastq/<S>_assemble-pass.fastq.gz` |
| что происходит | FASTQ → FASTA → IgBLAST (`-outfmt 19`, AIRR): для каждого рида V, D, J гены, границы FR/CDR, идентичность, продуктивность. Мышь: одинаковые риды сначала схлопываются (CollapseSeq), аннотируются уникальные, аннотация переносится на все риды; базы `mouse_gl_V/D/J`. Человек: аннотируется каждый рид, базы AIRR-C/OGRDB |
| выход | `results/<DS>/annotation/igblast/<S>.airr.tsv` (одна строка на склеенный рид) |
| промежуточное (мышь) | `annotation/unique/<S>.collapse-unique.fasta`, `<S>.unique.airr.tsv`, `annotation/index/<S>.unique_airr.sqlite` |
| отчёт | `annotation/report/` (сводки, распределения, HTML с igbrowser) |

## 5. Фильтр после аннотации — `filter_{mouse,human}_post_annotation.ipynb`

| | |
|---|---|
| вход | `results/<DS>/merged/fastq/<S>_assemble-pass.fastq.gz` + `results/<DS>/annotation/igblast/<S>.airr.tsv` (построчно, по `sequence_id`) |
| что происходит | рид остаётся, если: V с 1-й позиции germline и полный V…J до 3′-конца J; V…J ≥ 300 нт; локус ожидаемый (мышь IGH; человек IGH/IGK/IGL по библиотеке); продуктивный, без стоп-кодона; идентичность V ≥ 80%, J ≥ 85%; e-value V и J ≤ 1e−5 |
| выход | `post_annotation_filtered/fastq/<S>_filtered.fastq.gz` (прошедшие риды), `airr_pass/<S>.airr.tsv[.gz]` (их аннотация), отклонённые — `airr_reject/` (человек) или `reject_index/` (мышь) |
| сводка | `post_annotation_filtered/filter_summary.json` (правила и счётчики) |

## 6. Симуляция — `simulate_{mouse,human}_post_annotation_filtered_insilicoseq_150bp_novaseq.ipynb`

Один ноутбук, внутри — несколько стадий; каждая пишет свою папку в `results/<DS>/simulated/<ветка>/`.
Сценарий «фиксированные кратности ПЦР» — тот же ноутбук с суффиксом `_fixed_pcr_fold`. Модель и параметры —
`docs/simulation_library.md`, сценарии — `docs/simulation_scenarios.md`.

| стадия | вход | что происходит | выход |
|---|---|---|---|
| 00 эталон | `post_annotation_filtered/{fastq,airr_pass}`, `filter_summary.json`; человек — ещё `pr_trimmed/fastq` (имена праймеров) | риды одной молекулы сводятся в один шаблон: человек — консенсус по UMI, мышь — схлопывание ошибок (1 замена, родитель ≥3× обильнее). Шаблон = ампликон ПЦР1: V…J + фланги + праймеры, возвращённые на концы. Число исходных копий: человек — число UMI, мышь — риды после схлопывания. Проверенный шаблон (колонка `verified` в `template_qc.tsv`): человек — консенсус UMI по ≥2 ридам, мышь — ≥2 рида после схлопывания | `00_primary_truth/<образец>_templates.fasta`, `<образец>_template_qc.tsv`, `source_provenance.json` |
| 01 ПЦР1 | шаблоны и их копии | 30 циклов, у каждого шаблона своя эффективность | `01_pcr1/<образец>_pcr1_pool.tsv` (шаблон → копии) |
| 02 порция и фрагментация | пул ПЦР1 | случайная 1/10 пула; каждая молекула рвётся в одной случайной точке; куски проходят отбор по длине (200–350 нт) | `02_fragmentation/<образец>_library_summary.tsv`, `_fragment_length_model.tsv` |
| 03 ПЦР2 | куски | 5 циклов, у каждого вида куска своя эффективность | временный пул, после записи 04 удаляется (`03_pcr2/` пустая) |
| 04 секвенирование | пул ПЦР2, сырые счётчики пар (`raw/<DS>`) | на секвенатор 1/1000 продукта; бюджет пар = 3 × сырые пары (`rb3x`, у человека по локусам); пары раздаются по видам кусков пропорционально копиям | `04_read_allocation/<образец>_selected_fragments.fasta`, `_allocation.tsv`, `_read_counts.tsv` |
| 05–06 риды | фрагменты и число пар на каждый | InSilicoSeq (модель NovaSeq) читает с каждого фрагмента PE150; `05_fastq_native/` — рабочие файлы InSilicoSeq | `06_fastq_pe150/<образец>_R{1,2}.fastq.gz` |
| сводки | всё выше | контроль чисел по стадиям | `qc/final_qc.tsv`, `library_summary.tsv`, `template_summary*.tsv`, `*_template_reconstructability.tsv`; отчёт по стадиям — `notebooks/simulation_report.ipynb` |

Образец: у мыши — каждый из шести run отдельно; у человека — один объединённый образец `PRJEB30386_all_chains` (IGH+IGK+IGL).

```text
шаблон ── ПЦР1 ──► копии ── 1/10 ──► молекулы ── разрез ──► куски ── отбор 200–350 ──► ПЦР2 ── 1/1000 ──► пары PE150
[Vпр]──V…J──C──[Cпр]                       [Vпр]──V…    …J──C──[Cпр]                      R1 ──►      ◄── R2
```

## 7. Проверка симуляции — `validate_{mouse,human}_novaseq_post_annotation_filtered.ipynb`

| | |
|---|---|
| вход | `<ветка>/06_fastq_pe150/`, `04_read_allocation/`, `00_primary_truth/`, `qc/final_qc.tsv`; реальные пары — `pr_trimmed/fastq` (только прошедшие фильтр аннотации) |
| что происходит | равномерная подвыборка пар (человек 250 000, мышь по 50 000 на образец); симулированные риды выравниваются bowtie2 на свои фрагменты (должны лечь на 100%); реальные и симулированные — на одни шаблоны (ошибки, длина вставки); сравниваются качество, GC, дубликаты, доли локусов |
| выход | `<ветка>/validation/`: `validation_summary.json`, `mixture_qc.tsv`, `internal_invariants.tsv`, `subset_fastq/`, `fragment_alignment/`, `template_alignment/` |
| FastQC для рис. 11.1 | `validation/fastqc/*_sim_R{1,2}` (симуляция) и `*_raw150_R{1,2}` (сырые реальные пары тех же молекул, первые 150 нт) |

## 8. Сборка — `assemble_{trust4,mixcr,rnaspades,trinity,gather}_{mouse,human}.ipynb`

| | |
|---|---|
| вход | `<ветка>/06_fastq_pe150/<образец>_R{1,2}.fastq.gz` |
| что происходит | сборщик собирает контиги из PE150; результат приводится к одному формату — `contigs.fasta`. rnaSPAdes-ноутбук запускает три режима (`BCR_SPADES_MODE`: `rna`, `meta`, `rnaviral`) |
| выход | `<ветка>/assemblers/<сборщик>/<образец>/contigs.fasta`, рядом — родные файлы сборщика, `run.log`, `run_metadata.json` |
| папки сборщиков | `trust4`, `mixcr`, `rnaspades`, `spades_meta`, `spades_rnaviral`, `trinity`, `gather` |

## 9. Оценка сборки — `benchmark_assemblers_{mouse,human}.ipynb`

| | |
|---|---|
| вход | `<ветка>/00_primary_truth/` (эталон и признак «проверенный»), `<ветка>/assemblers/*/<образец>/contigs.fasta`, `06_fastq_pe150/` (поддержка ридами) |
| что происходит | эталон = V…J каждого шаблона с обрезанными краями (`trim20` мышь, `trim0` человек, `trim45_15` оба); одинаковые после обрезки V…J объединяются. Для каждой эталонной V…J ищется контиг, содержащий её целиком и точно (`exact_full`), с одной заменой, почти целиком, частично или никак; для каждого контига — что он такое (точный, мозаика, химера, ошибка) |
| выход | `<ветка>/benchmark_unique_trim{0,20,45_15}/assemblers_summary.tsv` (доли классов по сборщикам), `assemblers_precision.tsv`, `assemblers_structure_by_support.tsv`, `assemblers_support_thresholds.tsv`, `assemblers_summary_all_templates.tsv`; по сборщику и образцу — `<сборщик>/<образец>/per_template.tsv`, `per_contig.tsv` |

## 10. Сводка и отчёт — скрипты

| скрипт | вход | выход |
|---|---|---|
| `scripts/collect_evaluations.py` | `results/<DS>/simulated/*/benchmark_*` | `docs/evaluation_history/*.tsv`, отчёты в `docs/evaluation_history/` |
| `scripts/make_report_figures.py` | `results/<DS>/<стадия>/multiqc`, `annotation/report`, `simulated/<ветка>/{qc,validation/fastqc,benchmark_*}` | `docs/report/figures/*.png` |

## Вспомогательные ноутбуки (вне основной цепочки)

| ноутбук | что делает | выход |
|---|---|---|
| `read_geometry_{mouse,human}.ipynb` | распределения длин R1/R2, вставки и перекрытия по стадиям | `results/<DS>/read_geometry/` |
| `align_simulated_reads_bowtie2_human.ipynb` | выравнивание всех симулированных ридов на их фрагменты | `<ветка>/alignment/` |
| `audit_human_v_primer_trimming.ipynb` | что именно вырезал MaskPrimers из R1 и как это соотносится с границей V | `results/PRJEB30386/annotation/qc/` |
| `annotator_compare_mouse.ipynb`, `igblast_abstar_problematic_2seq.ipynb` | сравнение IgBLAST и abstar | `results/ERP003950/annotator_compare/`, `problematic_sequences/` |
| `assemble_vdjer_{mouse,human}.ipynb` | V'DJer; исключён из сравнения | `<ветка>/assemblers/vdjer/` |
