
### Full mouse IgBLAST quality

Full IgBLAST annotation outputs are being computed for all six mouse samples. As of now:

- ERR346596: complete (2.0 GB, 658,771 reads)
- ERR346597: complete (1.9 GB, 627,705 reads)
- ERR346598: in progress (30 MB so far)
- ERR346599, ERR346600, ERR346601: not yet started

A background summary computation is running on task for the two complete samples. The full-dataset stop-codon and V/J quality percentages will be appended to this report once all six samples finish.

## Decision

Do **not** drop abstar yet based only on the current evidence.

Reason: abstar has many abstar-only diagnostic/alignment fields (90 abstar-only columns), but the decisive two-problematic-sequence abstar test could not complete in the SSH memory cgroup. The current status is: IgBLAST full-dataset quality is being quantified; abstar usefulness remains pending a Jupyter/high-memory rerun on `problematic_mouse_2seq.fa`.

## Caveats

- V/J exact string agreement can be misleading due to database nomenclature differences.
- The current bad V/J statistic uses AA alignment identity with a provisional 85% threshold. If you want a stricter or looser definition of "маленькое identity", rerun the summary with the chosen threshold.
- E-value criterion `>1` could not be applied because `v_evalue`/`j_evalue` columns were absent from IgBLAST AIRR `outfmt 19` outputs for the problematic sequences (the IgBLAST outfmt 19 output includes `v_evalue`/`j_evalue` for full dataset but the two problematic sequences show them; the full-dataset summary will confirm whether they are present).
- abstar problematic-sequence run should be repeated in Jupyter/bcr_env rather than SSH.
- Full-dataset quality summaries are still in progress for ERR346598-ERR346601.
