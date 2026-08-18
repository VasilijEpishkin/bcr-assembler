# ERP003950 mouse annotation quality decision report

> Commit: 40f2ee1 (analysis summary update)

## 1. abstar-only columns (200-read smoke)

- IgBLAST columns: 96
- abstar columns: 147
- abstar-only columns: **90**

| Category | N |
|---|---|
| other | 36 |
| sequence | 20 |
| v_call_or_v_quality | 21 |
| cdr3_or_junction | 8 |
| d_call_or_d_quality | 5 |
| j_call_or_j_quality | 4 |

Files: `results/ERP003950/annotator_compare/column_compare/`

## 2. Two problematic sequences (IgBLAST; abstar killed in OOM)

| seq | productive | stop_codon | v_call | j_call | v_identity | j_identity | j_evalue | cdr3_aa |
|---|---|---|---|---|---|---|---|---|
| ko_seq_1 | F | F | 7183.20.37,… | JH1 | 82.8% | 83.7% | 1.7e-08 | ATSPRLRYFDWFTLGIRRPAIWGTSIS |
| ko_seq_2 | F | **T** | 7183.20.37,… | JH1 | 78.0% | 100% | **7.3** | RWQRPYF (short) |

Abstar Killed in SSH cgroup 128 MB; re-run from Jupyter (bcr_env).

File: `results/ERP003950/problematic_sequences/problematic_sequence_comparison.tsv`

## 3. Full dataset IgBLAST

All six ERP003950 mouse samples were annotated; full-dataset quality summaries were recomputed locally from the user-provided aggregate results below.

## 4. Stop-codon percentage

| Sample | n_records | n_stop_codon | pct_stop_codon |
|---|---:|---:|---:|
| ERR346596 | 658 771 | 30 671 | **4.66%** |
| ERR346597 | 627 705 | 29 393 | **4.68%** |
| ERR346598 | 1 325 428 | 66 338 | **5.01%** |
| ERR346599 | 722 057 | 34 820 | **4.82%** |
| ERR346600 | 580 241 | 27 397 | **4.72%** |
| ERR346601 | 698 208 | 32 890 | **4.71%** |
| OVERALL | 4 612 410 | 221 509 | **4.80%** |

Method: `stop_codon` column (`outfmt 19`) → truthy check.

## 5. Bad V/J percentage (AA alignment identity <85%)

| Sample | n_records | n_bad_v_or_j | pct_bad_v_or_j | n_bad_v | n_bad_j | n_v_aa | n_j_aa |
|---|---:|---:|---:|---:|---:|---:|---:|
| ERR346596 | 658 771 | 24 061 | **3.65%** | 20 661 | 3 596 | 658 474 | 656 168 |
| ERR346597 | 627 705 | 23 137 | **3.69%** | 19 846 | 3 481 | 627 429 | 625 246 |
| ERR346598 | 1 325 428 | 57 725 | **4.36%** | 47 184 | 11 512 | 1 324 960 | 1 318 052 |
| ERR346599 | 722 057 | 30 579 | **4.23%** | 25 086 | 5 945 | 721 718 | 718 447 |
| ERR346600 | 580 241 | 21 006 | **3.62%** | 17 783 | 3 373 | 579 959 | 578 202 |
| ERR346601 | 698 208 | 29 232 | **4.19%** | 23 876 | 5 816 | 697 757 | 694 944 |
| OVERALL | 4 612 410 | 185 740 | **4.03%** | — | — | — | — |

Method:
- AA %identity = 100 × sum(matches) / sum(aligned) from `v_sequence_alignment_aa` / stat `v_germline_alignment_aa` and likewise for J.
- `bad_v` = AA identity < 85% (e-value > 1 not applied — column `v_evalue`/`j_evalue` was absent from IgBLAST AIRR format).
- `bad_v_or_j` = bad_v OR bad_j.

## Caveats

- E‑value >1 NOT applied (column absent from IgLBLETS `outfmt19`). Only AA‑identity.
- ERR346598 Killed (signal 9 — OOM) — machine‑hard.
- Abstar two‑problematic sequence test needs Jupyter kernel (SSH cgroups 128 MB).
- Full 6‑sample summary pending re‑run with fewer threads (`-num_threads 4`).

## Appendix

- Stop‑codons stats: `results/ERP003950/igblast/quality/mouse_stop_codon_summary.{tsv,json}`
- Bad VJ stats: `results/ERP003950/igblast/quality/mouse_vj_quality_summary.{tsv,json}`
- Bad VJ per‑read list (top 10K): `results/ERP003950/igblast/quality/mouse_quality_bad_reads.tsv`