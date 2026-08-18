# ERP003950 mouse IgBLAST quality statistics

Source TSV files:

- `mouse_quality_summary_by_sample.tsv`
- `mouse_quality_summary_overall.tsv`
- `mouse_vj_quality_summary.tsv`
- `mouse_stop_codon_summary.tsv`

Thresholds:

- AA identity bad: `< 85%`
- IgBLAST support/e-value bad: `> 1`
- e-value fields: `IgBLAST_AIRR_v_support_j_support` (`v_support` / `j_support` in IgBLAST AIRR/Change-O output)

## TL;DR overall

| sample | annotated reads | stop codon | bad V/J AA identity <85% | bad V/J e-value >1 | bad V/J combined | productive=False |
| --- | --- | --- | --- | --- | --- | --- |
| ALL 6 samples | 4 612 410 | 221 509 (4.802%) | 185 740 (4.027%) | 35 082 (0.761%) | 209 783 (4.548%) | 256 674 (5.565%) |

Interpretation:

- **stop codon**: reads with `stop_codon=True`.
- **bad V/J AA identity**: V or J amino-acid alignment identity `<85%`.
- **bad V/J e-value**: `v_support` or `j_support >1`.
- **combined**: bad by AA identity **OR** by e-value/support.

## Main metrics by sample

| sample | annotated reads | stop codon | bad V/J AA identity <85% | bad V/J e-value >1 | bad V/J combined | productive=False |
| --- | --- | --- | --- | --- | --- | --- |
| ERR346596 | 658 771 | 30 671 (4.656%) | 24 061 (3.652%) | 5 035 (0.764%) | 27 473 (4.170%) | 35 762 (5.429%) |
| ERR346597 | 627 705 | 29 393 (4.683%) | 23 137 (3.686%) | 4 484 (0.714%) | 26 095 (4.157%) | 33 898 (5.400%) |
| ERR346598 | 1 325 428 | 66 338 (5.005%) | 57 725 (4.355%) | 12 287 (0.927%) | 66 632 (5.027%) | 77 543 (5.850%) |
| ERR346599 | 722 057 | 34 820 (4.822%) | 30 579 (4.235%) | 5 961 (0.826%) | 34 721 (4.809%) | 40 576 (5.620%) |
| ERR346600 | 580 241 | 27 397 (4.722%) | 21 006 (3.620%) | 2 717 (0.468%) | 22 800 (3.929%) | 30 964 (5.336%) |
| ERR346601 | 698 208 | 32 890 (4.711%) | 29 232 (4.187%) | 4 598 (0.659%) | 32 062 (4.592%) | 37 931 (5.433%) |

## V/J breakdown by sample

| sample | bad V AA | bad J AA | bad V or J AA | bad V e-value | bad J e-value | bad V or J e-value | combined V | combined J | combined V or J |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ERR346596 | 20 661 (3.136%) | 3 596 (0.546%) | 24 061 (3.652%) | 888 (0.135%) | 4 411 (0.670%) | 5 035 (0.764%) | 20 824 (3.161%) | 7 961 (1.208%) | 27 473 (4.170%) |
| ERR346597 | 19 846 (3.162%) | 3 481 (0.555%) | 23 137 (3.686%) | 809 (0.129%) | 3 909 (0.623%) | 4 484 (0.714%) | 19 963 (3.180%) | 7 332 (1.168%) | 26 095 (4.157%) |
| ERR346598 | 47 184 (3.560%) | 11 512 (0.869%) | 57 725 (4.355%) | 1 268 (0.096%) | 11 363 (0.857%) | 12 287 (0.927%) | 47 391 (3.576%) | 22 731 (1.715%) | 66 632 (5.027%) |
| ERR346599 | 25 086 (3.474%) | 5 945 (0.823%) | 30 579 (4.235%) | 881 (0.122%) | 5 316 (0.736%) | 5 961 (0.826%) | 25 244 (3.496%) | 11 184 (1.549%) | 34 721 (4.809%) |
| ERR346600 | 17 783 (3.065%) | 3 373 (0.581%) | 21 006 (3.620%) | 782 (0.135%) | 2 128 (0.367%) | 2 717 (0.468%) | 17 909 (3.086%) | 5 489 (0.946%) | 22 800 (3.929%) |
| ERR346601 | 23 876 (3.420%) | 5 816 (0.833%) | 29 232 (4.187%) | 1 200 (0.172%) | 3 699 (0.530%) | 4 598 (0.659%) | 24 092 (3.451%) | 9 480 (1.358%) | 32 062 (4.592%) |
| ALL 6 samples | 154 436 (3.348%) | 33 723 (0.731%) | 185 740 (4.027%) | 5 828 (0.126%) | 30 826 (0.668%) | 35 082 (0.761%) | 155 423 (3.370%) | 64 177 (1.391%) | 209 783 (4.548%) |

## Additional diagnostics

| sample | complete_vdj=False | missing V call | missing J call | V AA identity available | J AA identity available | V support available | J support available |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ERR346596 | 308 890 (46.889%) | 227 (0.034%) | 2 603 (0.395%) | 658 474 | 656 168 | 658 544 | 656 168 |
| ERR346597 | 293 813 (46.807%) | 220 (0.035%) | 2 459 (0.392%) | 627 429 | 625 246 | 627 485 | 625 246 |
| ERR346598 | 619 395 (46.732%) | 374 (0.028%) | 7 376 (0.556%) | 1 324 960 | 1 318 052 | 1 325 054 | 1 318 052 |
| ERR346599 | 334 891 (46.380%) | 257 (0.036%) | 3 610 (0.500%) | 721 718 | 718 447 | 721 800 | 718 447 |
| ERR346600 | 273 357 (47.111%) | 214 (0.037%) | 2 039 (0.351%) | 579 959 | 578 202 | 580 027 | 578 202 |
| ERR346601 | 319 932 (45.822%) | 350 (0.050%) | 3 264 (0.467%) | 697 757 | 694 944 | 697 858 | 694 944 |

## Files produced by notebook

- `mouse_quality_summary_by_sample.tsv` — rich per-sample table.
- `mouse_quality_summary_overall.tsv` — one-row overall table across all processed samples.
- `mouse_quality_summary.json` — full structured report.
- `mouse_bad_vj_reads.tsv` / `mouse_vj_quality_bad_reads.tsv` — capped bad-read detail table.
- `mouse_input_status.tsv` — input existence/size/processed status.

All 6 expected samples were processed: `ERR346596`, `ERR346597`, `ERR346598`, `ERR346599`, `ERR346600`, `ERR346601`.
