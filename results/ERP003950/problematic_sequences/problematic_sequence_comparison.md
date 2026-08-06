# Problematic sequence annotation comparison

## Inputs

- FASTA: `/data/user/epishkin/results/ERP003950/problematic_sequences/input/problematic_mouse_2seq.fa`
- IgBLAST TSV: `/data/user/epishkin/results/ERP003950/problematic_sequences/output/igblast/problematic_mouse_2seq_igblast.tsv`
- abstar TSV: `NOT AVAILABLE`

## Side-by-side summary

| sequence_id | in_igblast | in_abstar | abstar_status | igblast_productive | igblast_stop_codon | igblast_v_call | igblast_j_call | igblast_v_identity | igblast_j_identity | igblast_cdr3_aa |
|---|---:|---:|---|---|---|---|---|---:|---:|---|
| ko_seq_1 | True | False | NOT_AVAILABLE_abstar_killed_in_ssh_128MB_cgroup | F | F | 7183.20.37,VH7183.a19.31,VH7183.a21.35 | JH1 | 82.759 | 83.721 | ATSPRLRYFDWFTLGIRRPAIWGTSIS |
| ko_seq_2 | True | False | NOT_AVAILABLE_abstar_killed_in_ssh_128MB_cgroup | F | T | 7183.20.37,VH7183.a21.35,VH98-3G | JH1 | 78.049 | 100.000 | RWQRPYF |

## Interpretation

IgBLAST results:
- ko_seq_1: unproductive, no stop codon, V identity 82.8%, J identity 83.7%, CDR3_AA `ATSPRLRYFDWFTLGIRRPAIWGTSIS`
- ko_seq_2: unproductive, **has stop codon**, V identity 78.0%, J identity 100%, CDR3_AA `RWQRPYF` (short)

## abstar status

abstar was attempted with `abstar run ... --germline_database c57bl6 -o airr --n_processes 1 --verbose` but was killed in the SSH shell cgroup (128 MB memory limit). Re-run abstar from a Jupyter kernel in bcr_env to complete the abstar half.
