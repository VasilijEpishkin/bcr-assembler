# ERP003950 problematic sequences: updated IgBLAST vs abstar annotation report

Input FASTA is preserved with original mixed-case sequence letters. IgBLAST/abstar outputs still report uppercase sequences internally, so the annotation calls are unchanged after preserving case.

## TL;DR

- **IgBLAST annotated both sequences**.
- **abstar annotated only `ko_seq_1`**; `ko_seq_2` has no abstar row.
- `ko_seq_1`: non-productive, no stop codon; abstar reason: `junction does not end with conserved W`.
- `ko_seq_2`: non-productive, stop codon, V frameshift, low V identity, weak J e-value/support.

## Side-by-side calls

| sequence_id | IgBLAST found | IgBLAST productive | IgBLAST stop | IgBLAST V frameshift | IgBLAST V call | IgBLAST J call | IgBLAST CDR3_AA | abstar found | abstar productive | abstar stop | abstar V call | abstar J call | abstar CDR3_AA |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | yes | F | F | F | 7183.20.37,VH7183.a19.31,VH7183.a21.35 | JH1 | ATSPRLRYFDWFTLGIRRPAIWGTSIS | yes | false | false | IGHV0-D4D7*00 | IGHJ0-32C2*00 | ATSPRLRYFDWFTLGIRRPAIWGTSI |
| ko_seq_2 | yes | F | T | T | 7183.20.37,VH7183.a21.35,VH98-3G | JH1 | RWQRPYF | no | — | — | — | — | — |

## Identity and support/e-value

| sequence_id | IgBLAST V identity | IgBLAST J identity | IgBLAST V support/e-value | IgBLAST J support/e-value | abstar V identity NT | abstar V identity AA | abstar J identity NT | abstar J identity AA | abstar V support | abstar J support |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | 82.759% | 83.721% | 1.34e-82 | 1.72e-08 | 82.534% | 77.320% | — | — | 3.22e-83 | 2.61e-09 |
| ko_seq_2 | 78.049% | 100.000% | 2.43e-34 | 7.35 | — | — | — | — | — | — |

## Flags

| sequence_id | flags |
| --- | --- |
| ko_seq_1 | IgBLAST: non-productive<br>IgBLAST: V identity <85%<br>IgBLAST: J identity <85%<br>abstar: non-productive<br>abstar: junction does not end with conserved W |
| ko_seq_2 | IgBLAST: non-productive<br>IgBLAST: stop codon<br>IgBLAST: V frameshift<br>IgBLAST: V identity <85%<br>IgBLAST: J support/e-value >1<br>abstar: no output row |

## Conclusions

### ko_seq_1

Both tools mark it as non-productive and without stop codon. Support/e-values are strong in both tools. abstar gives the most interpretable cause: `junction does not end with conserved W`. CDR3 calls differ by one amino acid at the end: IgBLAST `ATSPRLRYFDWFTLGIRRPAIWGTSIS`, abstar `ATSPRLRYFDWFTLGIRRPAIWGTSI`.

### ko_seq_2

IgBLAST marks it as non-productive with stop codon and V frameshift. CDR3 is very short (`RWQRPYF`). V identity is low (`78.049%`). J identity is nominally `100%`, but J support/e-value is poor (`7.349 > 1`). abstar emits no row for this sequence.

Overall: IgBLAST remains the complete annotator for these problematic records; abstar is useful for extra diagnostics when it emits a row, but does not replace IgBLAST.
