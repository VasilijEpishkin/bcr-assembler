# bcr-assembler

`bcr-assembler` — исследовательский репозиторий для разработки и проверки пайплайна
воссоздания и последующей сборки BCR-последовательностей из paired-end reads.

Сейчас репозиторий описывает и отрабатывает только тот сценарий, который уже удалось
провести последовательно от начала до конца на одном датасете:

- датасет: `ERP003950`
- организм: мышь
- тип данных: BCR heavy-chain amplicon reads
- текущая подтверждённая траектория: raw FASTQ → QC → adapter trim → QC → primer trim
  → QC → merge paired-end reads в reference-like dataset → annotation/compare/quality
  analysis → simulation

Финальный шаг проекта — запуск repertoire assemblers уровня `TRUST4` на
симулированных/фрагментированных reads и оценка качества сборки относительно
референсного датасета — ещё в работе.

Важно различать две разные задачи:

1. `merge paired-end reads` — техническое восстановление одной последовательности из
   пары R1/R2 для формирования эталонного набора reads/templates.
2. `BCR assembly` — более поздняя задача реконструкции V(D)J/C-последовательностей из
   reads с помощью инструментов уровня `TRUST4`.

В этом репозитории первый шаг уже реализован для мыши, второй — исследуется.

## Основная идея проекта

Проект строится вокруг такой цепочки:

1. взять сырые paired-end BCR reads;
2. определить тип эксперимента и характеристики секвенирования;
3. убрать технический шум, не потеряв полезный иммунорепертуарный сигнал;
4. восстановить merged dataset, который можно использовать как reference-like baseline;
5. аннотировать этот baseline несколькими аннотаторами;
6. на его основе сгенерировать синтетические датасеты с контролируемым ground truth;
7. запускать assemblers и сравнивать их результат с известной reference truth.

Идея не в том, чтобы привязаться к одному инструменту, а в том, чтобы получить
проверяемый workflow для сравнения разных решений на BCR reads.

## Текущий подтверждённый объект: mouse dataset `ERP003950`

Именно мышиный датасет `ERP003950` сейчас является главным рабочим кейсом репозитория,
потому что по нему уже собрана непрерывная экспериментальная цепочка:

- сырые reads;
- QC на нескольких стадиях;
- adapter trimming;
- primer trimming;
- merge paired-end reads;
- annotation через IgBLAST;
- сравнение с abstar;
- quality summaries по annotation output;
- подготовка reference-like merged dataset;
- simulation через InSilicoSeq;
- подготовка входа для будущего тестирования `TRUST4` и подобных assembler'ов.

Другие датасеты и species-специфичные ноутбуки в репозитории остаются как exploratory
материал, но этот README описывает только мышиный пайплайн, потому что именно он сейчас
лучше всего подтверждён на практике.

## Пайплайн для мыши: от начала до конца

### 1. Сырые данные и первичное понимание эксперимента

Исходная точка — paired-end FASTQ мышиного датасета `ERP003950`.

Для этого датасета в репозитории уже зафиксированы рабочие предпосылки:

- mouse IgG heavy-chain dataset;
- paired-end geometry соответствует MiSeq/PE250-like данным;
- это не финальная задача de novo repertoire assembly, а сначала задача корректного
  восстановления reference-like merged dataset из исходных reads.

Смысл этого этапа — до любых преобразований понять:

- тип протокола;
- где likely adapters;
- где biological signal, который нельзя случайно срезать;
- какие технические праймеры нужно удалять отдельно.

### 2. QC сырых reads

Инструменты:

- `FastQC`
- `MultiQC`
- helper script: `scripts/qc.py`
- notebook: `notebooks/qc.ipynb`

Что делается:

- per-sample FastQC;
- сводка через MultiQC;
- проверка качества по циклам, длинам, adapter contamination и общему профилю run.

Цель:

- зафиксировать baseline качества до любых обрезок;
- понять, какие именно технические последовательности реально мешают дальше.

### 3. Определение стратегии adapter trimming

Для `ERP003950` стратегия не сводится к blind autodetect.

Подтверждённый текущий подход:

- notebook: `notebooks/adapter_trim_mouse.ipynb`
- делается только `adapter/quality trimming`;
- FR1/CH1 technical primers на этом этапе не режутся;
- после smoke-test и ручной проверки для датасета захардкожены явные
  `Illumina/TruSeq-style` адаптеры;
- фактическая обрезка выполняется `cutadapt`.

Ключевая идея этого этапа: сначала убрать sequencing/library adapters и length/quality
artefacts, не смешивая это с primer-trimming логикой.

### 4. Adapter trimming и техническая фильтрация

Основной рабочий mouse notebook:

- `notebooks/adapter_trim_mouse.ipynb`

Инструменты:

- `cutadapt`

Что делается:

- обрезаются R1/R2 Illumina-style adapters;
- применяется quality trimming;
- применяется length filtering;
- выход сохраняется как `*.trim.fastq.gz`.

Для общего/раннего workflow в репозитории есть и helper script
`scripts/adapter_trim.py`, где показан общий шаблон `cutadapt + fastp`, но для текущего
mouse-пайплайна подтверждённая логика зафиксирована именно в
`notebooks/adapter_trim_mouse.ipynb` и построена вокруг `cutadapt`.

### 5. QC после adapter trimming

После adapter/quality trimming снова запускается QC:

- notebook: `notebooks/qc.ipynb`
- helper script: `scripts/qc.py`

Зачем:

- проверить, что adapter contamination действительно ушла;
- оценить распределение длин и качество после фильтрации;
- убедиться, что данные готовы к следующему этапу.

### 6. Primer trimming

Для мыши текущая подтверждённая рабочая политика такая:

- notebook: `notebooks/primer_trim_mouse.ipynb`
- рабочий helper: `scripts/primer_trim_mouse_constant_only.py`
- обрезается только праймер constant-региона IgG;
- V-region forward primers сохраняются и не включаются в trim reference для этого шага.

Это важный результат проекта: primer trimming для mouse оказался не просто
"запустить MaskPrimers на всё подряд", а отдельной исследовательской задачей с выбором,
что именно нужно резать без потери полезной информации.

Текущая реализация:

- `cutadapt` используется для constant-region primer trimming;
- для разных samples заранее зафиксировано, в каком mate находится соответствующий
  constant-region primer;
- второй mate при необходимости копируется без изменений;
- выход сохраняется как `*.pr.fastq.gz`.

В репозитории остаётся и более общий `scripts/primer_trim.py` на `MaskPrimers.py`, но
mouse-пайплайн сейчас опирается на `primer_trim_mouse_constant_only.py`.

### 7. QC после primer trimming

После primer trim снова выполняется QC:

- `FastQC`
- `MultiQC`
- `notebooks/qc.ipynb`

Это нужно, чтобы убедиться, что технические последовательности убраны аккуратно и без
нежелательной деградации данных.

### 8. Merge paired-end reads для формирования reference-like dataset

Это центральный подтверждённый шаг текущего проекта.

Цель:

- не запустить ещё repertoire assembler,
- а сначала восстановить merged sequences из R1/R2, чтобы получить максимально
  надёжный reference-like dataset для дальнейшей annotation и simulation.

Основной notebook:

- `notebooks/presto_mouse.ipynb`

Основной инструмент:

- `pRESTO AssemblePairs.py`

Что в нём зафиксировано:

- используется оптимизированный mouse-specific merge workflow;
- `--nproc` задаётся явно, а не берётся дефолт pRESTO;
- verbose per-read `--log` по умолчанию выключен, чтобы не раздувать I/O;
- запись failed reads по умолчанию тоже выключена;
- long-running merge запускается с heartbeat и отдельными stdout/stderr log files.

Для экспериментов в репозитории также лежит альтернативный merge comparator:

- `scripts/merge_mouse_pairs_custom.py`

Этот скрипт реализует собственную overlap-based сборку и полезен как дополнительный
контрольный/исследовательский инструмент, но основной подтверждённый merge pipeline для
мыши сейчас — через `pRESTO`.

### 9. QC merged reads

После merge выполняется ещё один QC-этап:

- `FastQC`
- `MultiQC`
- `notebooks/qc.ipynb`

Цель:

- проверить профиль уже reconstructed reads;
- убедиться, что merged dataset пригоден как baseline/reference для следующего анализа.

### 10. Аннотация reference-like merged dataset

После получения merged reads запускается annotation.

Основной notebook:

- `notebooks/annotate_mouse.ipynb`

Основной инструмент:

- `IgBLAST`

Что делает notebook:

1. конвертирует merged FASTQ в FASTA;
2. запускает `igblastn` по каждому sample;
3. делает safe resume;
4. пропускает только действительно complete outputs;
5. архивирует stale/incomplete TSV и пересчитывает их;
6. использует `nproc=4` по умолчанию, потому что `nproc=8` уже приводил к OOM на
   `ERR346598`.

Это уже не просто preprocessing: здесь формируется аннотированный baseline, который
можно использовать как reference truth для последующих сравнений.

### 11. Quality analysis по аннотации

После IgBLAST в репозитории есть отдельный notebook:

- `notebooks/mouse_full_quality_summary_6samples.ipynb`

Он считает по всем 6 samples mouse dataset `ERP003950`:

- stop codon rate;
- долю reads с плохой V/J amino-acid identity;
- долю reads с плохими `v_support/j_support` значениями;
- combined quality metrics;
- дополнительные annotation-level показатели.

Это важный слой проекта: качество оценивается не только по FASTQ/QC-графикам, но и по
содержимому аннотации.

### 12. Сравнение аннотаторов

Чтобы не полагаться на один annotation stack, в репозитории есть сравнение:

- `notebooks/annotator_compare_mouse.ipynb`
- `notebooks/igblast_abstar_problematic_2seq.ipynb`

Инструменты:

- `IgBLAST`
- `abstar`

Что сравнивается:

- productivity;
- V/D/J calls;
- CDR3 AA;
- per-read agreement/disagreement;
- problematic sequences отдельным notebook'ом.

То есть reference-like dataset не только строится, но и интерпретируется через несколько
аннотаторов.

### 13. Подготовка reference dataset для simulation

После merge + annotation полученный mouse dataset используется как основа для
симуляции.

Смысл:

- сформировать controlled ground truth;
- получить synthetic reads, на которых потом можно честно сравнивать assemblers.

### 14. Simulation sequencing reads

В репозитории уже есть несколько mouse simulation notebooks:

- `notebooks/simulate_mouse_merged_insilicoseq.ipynb`
- `notebooks/simulate_mouse_merged_insilicoseq_150bp.ipynb`
- `notebooks/simulate_mouse_merged_insilicoseq_150bp_pcr_fragmented.ipynb`

Основные инструменты:

- `InSilicoSeq`
- `bowtie2` для подготовки/обучения error model и привязки реальных reads к merged truth

Что уже делается:

- строится общий reference FASTA из merged mouse reads;
- используются реальные mouse merged sequences как templates;
- моделируется sequencing error profile;
- генерируются synthetic paired-end datasets;
- отдельно исследуется 150bp сценарий;
- отдельно исследуется PCR-like amplification + fragmentation сценарий.

Это уже прямая подготовка к следующему шагу — тестированию repertoire assemblers.

### 15. Следующий этап: assemblers уровня TRUST4

Это то, к чему проект идёт сейчас.

Цель следующего шага:

- взять simulated fragmented reads;
- прогнать `TRUST4` и сопоставимые инструменты;
- сравнить восстановленные последовательности с известной reference truth,
  полученной на предыдущих этапах.

Именно этот шаг в проекте ещё не завершён, поэтому текущий README честно фиксирует:

- reference-like merged mouse dataset уже построен;
- annotation и simulation уже идут;
- полноценная оценка assemblers относительно reference dataset — следующий основной
  milestone.

## Что уже сделано и подтверждено

Для mouse dataset `ERP003950` в репозитории уже есть подтверждённые рабочие этапы:

1. raw QC;
2. adapter trimming;
3. QC после adapter trimming;
4. primer trimming с текущей mouse-specific политикой;
5. QC после primer trimming;
6. merge paired-end reads через `pRESTO`;
7. QC merged reads;
8. full-dataset annotation через `IgBLAST`;
9. quality summaries по annotation output;
10. сравнение `IgBLAST` vs `abstar`;
11. problematic-sequence analysis;
12. simulation на основе reference-like merged dataset.

## Что ещё исследуется

1. насколько текущая mouse-specific primer trimming политика оптимальна;
2. какие merge/annotation/quality thresholds лучше держать как канонические;
3. как лучше моделировать realistic fragmented BCR reads;
4. как сравнивать `TRUST4` и другие assemblers на controlled ground truth;
5. какие метрики лучше всего отражают качество конечной BCR assembly.

## Основные файлы и ноутбуки

### Core notebooks для мыши

- `notebooks/qc.ipynb` — QC на разных стадиях
- `notebooks/adapter_trim_mouse.ipynb` — adapter/quality trimming для `ERP003950`
- `notebooks/primer_trim_mouse.ipynb` — текущая mouse-specific primer trim логика
- `notebooks/presto_mouse.ipynb` — основной merge paired-end workflow
- `notebooks/annotate_mouse.ipynb` — full IgBLAST annotation
- `notebooks/mouse_full_quality_summary_6samples.ipynb` — quality summary по annotation
- `notebooks/annotator_compare_mouse.ipynb` — IgBLAST vs abstar
- `notebooks/igblast_abstar_problematic_2seq.ipynb` — разбор проблемных sequences
- `notebooks/simulate_mouse_merged_insilicoseq.ipynb` — базовая simulation
- `notebooks/simulate_mouse_merged_insilicoseq_150bp.ipynb` — 150bp fragmented simulation
- `notebooks/simulate_mouse_merged_insilicoseq_150bp_pcr_fragmented.ipynb` — PCR-like
  fragmented simulation

### Helper scripts

- `scripts/qc.py` — FastQC + MultiQC runner
- `scripts/adapter_trim.py` — общий adapter trim helper
- `scripts/primer_trim.py` — общий MaskPrimers-based helper
- `scripts/primer_trim_mouse_constant_only.py` — текущий mouse-specific primer trim helper
- `scripts/merge_mouse_pairs_custom.py` — альтернативный custom merge comparator

## Используемые инструменты

Подтверждённые по текущему mouse workflow инструменты:

- `FastQC`
- `MultiQC`
- `cutadapt`
- `pRESTO`
- `IgBLAST`
- `abstar`
- `InSilicoSeq`
- `bowtie2`

Целевые/следующие инструменты для основной assembly-задачи:

- `TRUST4`
- другие repertoire assemblers для сравнения

## Что хранится в репозитории

В git хранятся:

- notebooks;
- helper scripts;
- reference sequences;
- QC/summary/report artefacts;
- стабильные вспомогательные таблицы и отчёты.

Смысл репозитория — сохранить воспроизводимую логику исследования и проверяемые
артефакты по mouse pipeline.

## Кратко в одной фразе

`bcr-assembler` — это репозиторий, в котором на mouse dataset `ERP003950`
последовательно строится и проверяется пайплайн от сырых BCR paired-end reads до
reference-like merged/annotated baseline, на основе которого затем можно честно
тестировать `TRUST4` и другие BCR assemblers.