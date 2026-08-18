# ERP003950 problematic sequences: IgBLAST vs abstar with FWR/CDR regions

Input FASTA preserved as mixed-case sequence; TSV outputs contain the annotator-normalized uppercase sequence.

## Main calls

| sequence_id | IgBLAST found | IgBLAST productive | IgBLAST stop | IgBLAST V frameshift | IgBLAST V call | IgBLAST J call | abstar found | abstar productive | abstar stop | abstar V call | abstar J call |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | yes | F | F | F | 7183.20.37,VH7183.a19.31,VH7183.a21.35 | JH1 | yes | false | false | IGHV0-D4D7*00 | IGHJ0-32C2*00 |
| ko_seq_2 | yes | F | T | T | 7183.20.37,VH7183.a21.35,VH98-3G | JH1 | no | — | — | — | — |

## NT regions

| sequence_id | tool | fwr1 | cdr1 | fwr2 | cdr2 | fwr3 | cdr3 | fwr4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | IgBLAST | GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCT | GGATTCACCTTTAGCAACTATGCC | ATGAGCTGGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTCTCAGCT | ATTACCGGTGGGGGTAGAAGGACA | TACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTACTTCTGT | GCGACCTCCCCCCGATTACGATATTTTGACTGGTTCACTCTTGGAATTCGCCGCCCCGCCATATGGGGTACTTCGATCTC | TGGGGCCGTGGCACCCTGGTCACCGTCTC |
| ko_seq_1 | abstar | GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCT | GGATTCACCTTTAGCAACTATGCC | ATGAGCTGGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTCTCAGCT | ATTACCGGTGGGGGTAGAAGGACA | TACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTACTTCTGT | GCGACCTCCCCCCGATTACGATATTTTGACTGGTTCACTCTTGGAATTCGCCGCCCCGCCATATGGGGTACTTCGATCTC | TGGGGCCGTGGCACCCTGGTCACCGTCTCGAGT |
| ko_seq_2 | IgBLAST | GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCTGTGCAGCCTCT | GGATTCTCTTTACGGACATGCC | ATGAGCTGGGTCCGCAGGCTCCCGGGAGGGGCTGGAGTGGGTCTCAGCT | ATTACGGGGGGGTAAGAGGACA | TACTACGCAGACTCGTGAGGGGCGGTTACATCTCAGAGACATCAGAAACCCTGTATCTGCAAGACGCTGAAACAAGAACGCTGTTACTCTGTGCGAAGAGGACTCGGATGACT | ATCGCTGGCAGAGACCGTACTTC | — |
| ko_seq_2 | abstar | — | — | — | — | — | — | — |

## AA regions

| sequence_id | tool | fwr1_aa | cdr1_aa | fwr2_aa | cdr2_aa | fwr3_aa | cdr3_aa | fwr4_aa |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | IgBLAST | EVQLLESGGGLVQPGGSLRLSCAAS | GFTFSNYA | MSWVRQAPGKGLEWVSA | ITGGGRRT | YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYFC | ATSPRLRYFDWFTLGIRRPAIWGTSIS | GAVAPWSPS |
| ko_seq_1 | abstar | EVQLLESGGGLVQPGGSLRLSCAAS | GFTFSNYA | MSWVRQAPGKGLEWVSA | ITGGGRRT | YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYFC | ATSPRLRYFDWFTLGIRRPAIWGTSI | WGRGTLVTVSS |
| ko_seq_2 | IgBLAST | EVQLLESGGGLVQPGGSLRLSVQPL | DSLYGHA | MSWVRRLPGGAGVGLS | YYGGVRGH | TTQTREGRLHLRDIRNPVSARR*NKNAVTLCEEDSDD | RWQRPYF | — |
| ko_seq_2 | abstar | — | — | — | — | — | — | — |

## Identity / support

| sequence_id | tool | V identity | J identity | V support/e-value | J support/e-value |
| --- | --- | --- | --- | --- | --- |
| ko_seq_1 | IgBLAST | 82.759% | 83.721% | 1.34e-82 | 1.72e-08 |
| ko_seq_1 | abstar | 82.534% | — | 3.22e-83 | 2.61e-09 |
| ko_seq_2 | IgBLAST | 78.049% | 100.000% | 2.43e-34 | 7.35 |
| ko_seq_2 | abstar | — | — | — | — |

## Interpretation

| sequence_id | notes |
| --- | --- |
| ko_seq_1 | IgBLAST: non-productive<br>abstar: junction does not end with conserved W<br>abstar: non-productive |
| ko_seq_2 | IgBLAST: non-productive<br>IgBLAST: stop codon<br>IgBLAST: V frameshift<br>abstar: no output row |

Notes:
- Empty `—` means the field is absent/empty in that annotator output.
- `ko_seq_2` is absent from abstar output, so all abstar region fields are `—`.
- `ko_seq_2` in IgBLAST has partial framework extraction: FWR1/2/3 are present, FWR4 is empty.
