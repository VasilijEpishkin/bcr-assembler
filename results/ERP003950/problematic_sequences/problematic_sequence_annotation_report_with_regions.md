# ERP003950 problematic sequences: IgBLAST vs abstar with FR/CDR regions

Input FASTA preserved the original mixed-case sequence letters. IgBLAST/abstar outputs internally display uppercase sequences, but the original FASTA file was refreshed before annotation.

Source artifacts:

- `input/problematic_mouse_2seq.fa`
- `output/igblast/problematic_mouse_2seq_igblast.tsv`
- `output/abstar/problematic_mouse_2seq_abstar.tsv`
- `problematic_sequence_comparison_full.tsv`

## Input sequences

- `ko_seq_1`: 401 nt
- `ko_seq_2`: 357 nt

## Annotation calls

| sequence_id | IgBLAST found | IgBLAST productive | IgBLAST stop | IgBLAST V frameshift | IgBLAST V call | IgBLAST J call | IgBLAST junction_aa | abstar found | abstar productive | abstar stop | abstar V call | abstar J call | abstar junction_aa |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | yes | F | F | F | 7183.20.37,VH7183.a19.31,VH7183.a21.35 | JH1 | CATSPRLRYFDWFTLGIRRPAIWGTSISG | yes | false | false | IGHV0-D4D7*00 | IGHJ0-32C2*00 | CATSPRLRYFDWFTLGIRRPAIWGTSIS |
| ko_seq_2 | yes | F | T | T | 7183.20.37,VH7183.a21.35,VH98-3G | JH1 | YRWQRPYFR | no | — | — | — | — | — |

## FR/CDR amino-acid regions

| sequence_id | IgBLAST fwr1_aa | IgBLAST cdr1_aa | IgBLAST fwr2_aa | IgBLAST cdr2_aa | IgBLAST fwr3_aa | IgBLAST cdr3_aa | IgBLAST fwr4_aa | abstar found | abstar fwr1_aa | abstar cdr1_aa | abstar fwr2_aa | abstar cdr2_aa | abstar fwr3_aa | abstar cdr3_aa | abstar fwr4_aa |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | EVQLLESGGGLVQPGGSLRLSCAAS | GFTFSNYA | MSWVRQAPGKGLEWVSA | ITGGGRRT | YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYFC | ATSPRLRYFDWFTLGIRRPAIWGTSIS | GAVAPWSPS | yes | EVQLLESGGGLVQPGGSLRLSCAAS | GFTFSNYA | MSWVRQAPGKGLEWVSA | ITGGGRRT | YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYFC | ATSPRLRYFDWFTLGIRRPAIWGTSI | WGRGTLVTVSS |
| ko_seq_2 | EVQLLESGGGLVQPGGSLRLSVQPL | DSLYGHA | MSWVRRLPGGAGVGLS | YYGGVRGH | TTQTREGRLHLRDIRNPVSARR*NKNAVTLCEEDSDD | RWQRPYF | — | no | — | — | — | — | — | — | — |

## FR/CDR nucleotide regions

| sequence_id | IgBLAST fwr1 | IgBLAST cdr1 | IgBLAST fwr2 | IgBLAST cdr2 | IgBLAST fwr3 | IgBLAST cdr3 | IgBLAST fwr4 | abstar found | abstar fwr1 | abstar cdr1 | abstar fwr2 | abstar cdr2 | abstar fwr3 | abstar cdr3 | abstar fwr4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCT | GGATTCACCTTTAGCAACTATGCC | ATGAGCTGGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTCTCAGCT | ATTACCGGTGGGGGTAGAAGGACA | TACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTACTTCTGT | GCGACCTCCCCCCGATTACGATATTTTGACTGGTTCACTCTTGGAATTCGCCGCCCCGCCATATGGGGTACTTCGATCTC | TGGGGCCGTGGCACCCTGGTCACCGTCTC | yes | GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCT | GGATTCACCTTTAGCAACTATGCC | ATGAGCTGGGTCCGCCAGGCTCCCGGGAAGGGGCTGGAGTGGGTCTCAGCT | ATTACCGGTGGGGGTAGAAGGACA | TACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTACTTCTGT | GCGACCTCCCCCCGATTACGATATTTTGACTGGTTCACTCTTGGAATTCGCCGCCCCGCCATATGGGGTACTTCGATCTC | TGGGGCCGTGGCACCCTGGTCACCGTCTCGAGT |
| ko_seq_2 | GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCTGTGCAGCCTCT | GGATTCTCTTTACGGACATGCC | ATGAGCTGGGTCCGCAGGCTCCCGGGAGGGGCTGGAGTGGGTCTCAGCT | ATTACGGGGGGGTAAGAGGACA | TACTACGCAGACTCGTGAGGGGCGGTTACATCTCAGAGACATCAGAAACCCTGTATCTGCAAGACGCTGAAACAAGAACGCTGTTACTCTGTGCGAAGAGGACTCGGATGACT | ATCGCTGGCAGAGACCGTACTTC | — | no | — | — | — | — | — | — | — |

## Identity / support metrics

| sequence_id | IgBLAST V identity | IgBLAST J identity | IgBLAST V support/evalue | IgBLAST J support/evalue | abstar V identity NT | abstar V identity AA | abstar J identity NT | abstar J identity AA | abstar V support | abstar J support | abstar productivity_issues |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ko_seq_1 | 82.759% | 83.721% | 1.34e-82 | 1.72e-08 | 82.534% | 77.320% | — | — | 3.22e-83 | 2.61e-09 | junction does not end with conserved W |
| ko_seq_2 | 78.049% | 100.000% | 2.43e-34 | 7.35 | — | — | — | — | — | — | — |

## Interpretation

### ko_seq_1

- IgBLAST and abstar both found it.
- Both call it non-productive, without stop codon and without V frameshift.
- FR1/CDR1/FR2/CDR2/FR3 regions are present in both tools.
- FWR4 differs: IgBLAST reports `GAVAPWSPS`, abstar reports `WGRGTLVTVSS`.
- abstar gives an explicit reason: `junction does not end with conserved W`.

### ko_seq_2

- IgBLAST found it; abstar emitted no output row.
- IgBLAST marks it as non-productive with stop codon and V frameshift.
- FR3 contains a stop codon symbol (`*`) and truncated sequence; FWR4 is absent.
- J identity is nominal `100%`, but J support/e-value is poor (`7.349 > 1`).
