# Результаты симуляции PRJEB30386: NovaSeq PE150

Ветка: `insilicoseq_150bp_novaseq_post_annotation_filtered`

## Итоговый набор данных

| Метрика | Значение |
|---|---:|
| Статус итоговой проверки | `valid = true` |
| Сгенерировано read pairs | 4 552 714 |
| R1 reads | 4 552 714 |
| R2 reads | 4 552 714 |
| Длина R1/R2 | 150 / 150 нт |
| Уникальные исходные templates | 549 820 |
| Молекулы после PCR1 | 10 942 016 500 292 |
| Fragment species | 2 749 100 |
| Fragment species, получившие reads | 1 110 858 |
| Fragment sampling dropout | 59,59% |
| Сохранение числа молекул при фрагментации | `true` |

## Состав по локусам

| Локус | Уникальные templates | Ожидаемые read pairs | Полученные read pairs | Доля | Отклонение |
|---|---:|---:|---:|---:|---:|
| IGH | 314 476 | 2 509 309 | 2 509 309 | 55,12% | 0 |
| IGK | 134 406 | 958 261 | 958 261 | 21,05% | 0 |
| IGL | 100 938 | 1 085 144 | 1 085 144 | 23,84% | 0 |
| **Всего** | **549 820** | **4 552 714** | **4 552 714** | **100%** | **0** |

## Проверка выравниванием

Validation выполнена на детерминированной подвыборке из 250 000 read pairs.

| Референс для проверки | Mapped | Properly paired | Error rate |
|---|---:|---:|---:|
| Точные simulated fragments | 100,00% | 100,00% | 0,001618 |
| Исходные V–J templates | 100,00% | 100,00% | 0,001521 |

Дополнительные метрики идентификации происхождения read:

| Метрика | Значение |
|---|---:|
| Проверено primary reads | 500 000 |
| Успешный разбор origin из read name | 100,00% |
| Exact fragment-record origin | 39,22% |
| Sequence-equivalent fragment origin | 39,36% |

Последние две метрики отражают разрешение конкретной fragment-записи среди повторяющихся/эквивалентных фрагментов и не заменяют показатели mapped/properly paired.

## FastQC / MultiQC

| Read | Reads | GC | Средняя длина | Медианная длина | Duplicate reads | Доля failed FastQC modules |
|---|---:|---:|---:|---:|---:|---:|
| R1 | 4,553 млн | 53% | 150 | 150 | 89,02% | 10% |
| R2 | 4,553 млн | 53% | 150 | 150 | 85,21% | 10% |

Высокая duplication ожидаема для AIRR repertoire с PCR-like abundance model и не означает автоматически низкое качество симуляции.

## Аудит входных библиотек

| Run | Library | Locus | Прошло post-annotation filter | Допущено к симуляции | Исключено из-за отсутствующего barcode | Другие нарушения |
|---|---|---|---:|---:|---:|---:|
| ERR3004229 | IgM | IGH | 1 090 122 | 1 079 373 | 10 748 | 1 sequence с `N` |
| ERR3004230 | IgG | IGH | 782 331 | 774 016 | 8 315 | 0 |
| ERR3004231 | IgK | IGK | 651 862 | 646 730 | 5 132 | 0 |
| ERR3004232 | IgL | IGL | 712 731 | 704 601 | 8 130 | 0 |

## Где смотреть результаты

- Итоговые FASTQ: `06_fastq_pe150/PRJEB30386_all_chains_R{1,2}.fastq.gz`
- Итоговая проверка: `validation/validation_summary.json`
- Проверка смеси локусов: `validation/mixture_qc.tsv`
- Метрики стадий: `qc/pcr1_fragmentation_pcr2_summary.tsv`
- MultiQC: `qc/multiqc/PRJEB30386_simulated_insilicoseq_150bp_novaseq_post_annotation_filtered_multiqc.html`
