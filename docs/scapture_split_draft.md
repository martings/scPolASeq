# SCAPTURE discovery split: status and real-data gate

This branch introduces `SCAPTURE_DISCOVERY` as a separate DSL2 subworkflow and
connects it through the opt-in `--scapture_execution_mode staged` path in
`APA_CORE`. `SCAPTURE_FILTER` remains the default production path while the
complete staged PBMC1k run and its `-resume` equivalence gate are pending. A
stub run cannot establish this equivalence.

## Stage contracts

| Stage | Inputs | Output boundary | Cache benefit |
| --- | --- | --- | --- |
| `SCAPTURE_PREPARE` | library metadata, GTF | collision-free GTF and gene map | GTF rewrite survives downstream errors |
| `SCAPTURE_ANNOTATION` | prepared GTF, FASTA/FAI, chromosome sizes | `${library_id}_annot.*`, log | annotation survives DeepPASS and merge failures |
| `SCAPTURE_PASCALL` | indexed BAM, FASTA/FAI, annotation files, optional poly(A) DB | three evaluated BEDs, counts, log | DeepPASS survives merge and quantification failures |
| `SCAPTURE_PASMERGE` | one evaluated BED and category | Integrated.bed plus explicit status | Each category resumes independently |
| `SCAPTURE_PASREFERENCE` | three category Integrated.bed files | sorted BED12 reference and row counts | PASquant can rerun without PASmerge |

The process modules have `meta.yml` contracts and `versions.yml` output. In the
container, SCAPTURE's Git commit is recorded when available, alongside a script
SHA256 fingerprint. The current container recipe clones a moving upstream HEAD;
pin that commit before proposing these as portable nf-core modules. Their
`meta.library_id` key is also a local pipeline contract; portable nf-core
modules should use `meta.id` or `task.ext.prefix`.
The GTF helper and patched DeepPASS predictor are explicit staged process
inputs, so changes to either invalidate the relevant Nextflow task hash.

`SCAPTURE_PASCALL` fails if no category has a valid evaluated BED row or if no
row has a non-`NA` DeepPASS label. The label is not a numeric probability. An
absent 3primeExtended category is represented by an empty BED and noted in the
log; the other categories still have to pass validation.

## Frozen PBMC1k reference for the real gate

- Historical session: `c3d2f281-2d73-4d34-99c6-cc876610548a`.
- Historical monolithic task: `85/5407a938ac61cb168479ca6cc8a528` under
  `/scratch/work/scPolASeq_runs/pbmc1k_20261024_replicate`.
- Successful monolith yielded 28,464 exonic, 70,568 intronic, and 2,366
  3primeExtended scored rows. The source task and published output are read-only
  evidence, not a replacement for a fresh split run.

Before promoting the opt-in staged path to the default, run the split on the same frozen
inputs and compare evaluated BED row counts and content hashes, the normalized
catalog, PASquant matrix, and manifest. Deliberately fail a downstream stage,
then rerun with the same `-work-dir` and `-resume`; verify the PAScall task hash
is cached and no DeepPASS prediction executes again. Use a new output directory.
Track PASmerge, PASquant, and orchestration in issues #27, #28, and #29.

PASquant in the current monolith uses `|| true`; the split must define whether
quantification is mandatory, then fail on missing or corrupt output rather than
silently reporting success. The integration must preserve any intentional
optional mode explicitly in the manifest.

## Verification completed for this increment

`tests/scapture_discovery_contract.nf` compiles the full DSL2 graph with an empty
channel. This only verifies wiring and parsing; it creates no biological files.
`tests/scapture_prepare_real.nf` executes the real GTF normalizer on two genes
sharing a name. It produced an explicit two-gene collision report and distinct
gene-id tokens. A second `nextflow run ... -resume` reported the same
`4d/652cc8...` task as `Cached process`, demonstrating hash reuse for the
preparation stage. That first test used host Python. A later bounded run used
the validated `/scratch/cache/apptainer/scpolaseq-scapture.sif` on a 1,200 bp
synthetic chromosome. `tests/scapture_annotation_real.nf` completed PREPARE at
`f1/d2f528...` and ANNOTATION at `50/006eff...`; both were cached on the next
`-resume`. The SIF used Python 3.8.19 and its SCAPTURE checkout was
`8845e84aa1b1f7c66759074437fe4df774f89318`. Annotation emitted six
nonempty files (`tiny_annot.*`). This is an actual SCAPTURE annotation run,
not a stub and not PBMC1k equivalence.

`tests/scapture_discovery_fail_closed.nf` supplies the same tiny reference with
an intentionally invalid BAM. PREPARE and ANNOTATION completed, while PAScall
exited 1 at `31/4a6b8b...` because no valid scored peaks existed. This checks
the subworkflow joins, predictor staging and Apptainer bind, and failure
propagation. The deliberately invalid BAM cannot test DeepPASS prediction or
scientific concordance. A full real-input PAScall run remains issue #29.

## Bounded real PASmerge gate (#27)

`tests/scapture_merge_historical.nf` takes `--scapture_evaluated_dir` and uses
the three already-evaluated BEDs from the read-only historical task; it does not
rerun PAScall or DeepPASS. With the validated SCAPTURE SIF, each category ran
separately and the combined BED12 reference was written. The category row
counts were 28,452 exonic, 70,568 intronic, and 2,366 3primeExtended; total
101,386. SHA256 matched the historical monolith output **exactly for each of
the four BED files**. The combined BED12 SHA256 was
`80afa97bdfe64f0b01f1200905ccbe97036f91087f86bace9a0e37614724e45d`.
The first attempt failed closed because GNU Parallel needed a task-local
`TMPDIR`; the module now creates one, matching the monolith's environment.

An unchanged `-resume` reused all four task hashes: exonic `bd/96b2ed`,
intronic `3f/81ec61`, 3primeExtended `49/abe153`, combined `a3/3c4796`.
With the final test harness revision, switching only 3primeExtended to an empty
BED produced `empty_input` status and a valid 99,020-row reference; switching
back cached exonic `07/32e96f` and intronic `24/8c5ef4`, reran only the changed
3primeExtended category (`c5/44642a`) and combined reference (`6f/97d998`).
The final combined SHA256 again matched the historical monolith. When all three
categories were empty, the category tasks recorded `empty_input` and the
reference task failed with `all SCAPTURE PASmerge categories were empty`.

Independent review found two portability/cache refinements. The merge
subworkflow now orders categories explicitly before passing BEDs to the
reference task, so asynchronous completion cannot perturb its task hash. The
reference process now selects the same SCAPTURE Conda/SIF environment as the
merge processes instead of depending on host GNU `sort`. A subsequent real
PBMC1k rerun produced combined reference hash `7b/f51bf2`, with SHA256 still
`80afa97bdfe64f0b01f1200905ccbe97036f91087f86bace9a0e37614724e45d`;
the next explicit `-resume ba59fd49-47ff-492c-949c-71c2b7d43e69` cached
all four PASmerge/reference tasks, including `7b/f51bf2`.

## Bounded real PASquant gate (#28)

`tests/scapture_quant_historical.nf` takes a completed SCAPTURE task directory
and runs the catalog normalizer and PASquant as separate processes on frozen
PBMC1k inputs. Its final gate can intentionally fail after PASquant. It does
not run PAScall/DeepPASS or change the historical output directory. The
historical source is session `c3d2f281-2d73-4d34-99c6-cc876610548a`, task
`85/5407a938ac61cb168479ca6cc8a528`.

The catalog was computed in Apptainer with 21,906 rows and SHA256
`1621066b0cd864a589b45f7d882ff20cbbae594e5a537bb5adbf00a22d507dc4`,
identical to the historical published catalog. The first PASquant attempt
completed the actual SCAPTURE command and wrote a 101,386-site by 1,139-cell
matrix. Its compressed SHA256 was
`b16c11fe3e71f8a85a475a4c5535f9035930106b0f51cdce5aa54fdba79cafa9`;
the uncompressed content SHA256 was
`882ff8ef36b6abf7333a00fdd34379b2e642587aae996e2af7ae20f8d7ba1e95`.
Both matched the historical monolith byte for byte. This is a real scientific
equivalence result, but the Nextflow task failed **after** SCAPTURE completed:
an AWK metric-writing statement expanded a newline into its program string.
The module was corrected to emit metrics with `print` and a final production
run is required to establish a successful Nextflow task and cache reuse.
A stub cannot satisfy this gate.

### PASquant execution and resume evidence (2026-10-27 UTC)

Nextflow 25.10.4 used the validated
`/scratch/cache/apptainer/scpolaseq-scapture.sif` (SCAPTURE checkout
`8845e84aa1b1f7c66759074437fe4df774f89318`). The test work directory was
`/tmp/scpolaseq_pasquant_20261027_work`, independent of the historical work
and publication directories. For this test, `-resume` takes the **session UUID**
`e8c8637d-60ad-49cc-b6f0-e80385801dc8`; task hashes are used only to
inspect which processes were reused.

```bash
nextflow run tests/scapture_quant_historical.nf \
  -profile deepthought,apptainer \
  --scapture_task_dir /scratch/work/scPolASeq_runs/pbmc1k_20261024_replicate/85/5407a938ac61cb168479ca6cc8a528 \
  --scapture_test_fail_after_quant true \
  --outdir /tmp/scpolaseq_pasquant_20261027_out \
  --apptainer_cache_dir /scratch/cache/apptainer \
  -work-dir /tmp/scpolaseq_pasquant_20261027_work \
  -resume e8c8637d-60ad-49cc-b6f0-e80385801dc8
```

After PASquant completed at task `60/5c5bd1` with exit 0, the deliberately
failing `SCAPTURE_TEST_GATE` exited 42 at `06/ee25b8`. Changing only
`--scapture_test_fail_after_quant false` and rerunning the command with the
same `-resume` UUID completed in 5.8 seconds: `SCAPTURE_CATALOG` `9d/4fee47`
was cached, `SCAPTURE_PASQUANT` `60/5c5bd1` was cached, and the gate completed
at `4b/ab93a6`. The final matrix gzip SHA256 was
`b16c11fe3e71f8a85a475a4c5535f9035930106b0f51cdce5aa54fdba79cafa9`,
identical to the historical published matrix. The PASquant manifest reports
`status=complete`, 1,139 barcodes, 101,386 site rows, and 1,139 cell columns.

A later fix preserved an explicitly supplied DeepPASS threshold of zero. A
catalog-only real-input run with `--scapture_deeppass_threshold 0` completed
at `ab/a71f72` and recorded `deeppass_threshold=0.0` with 99,024 catalog
rows. Under the default 0.5 threshold, the final catalog task `36/f0a4d5`
produced 21,906 rows and the historical SHA256 above; PASquant still resumed
from `60/5c5bd1`. An invalid `--scapture_test_sample_id bad/name` was
rejected before PASquant launched, exercising the filename/shell guard.

This proves the split catalog and PASquant stages on frozen PBMC1k inputs,
including scientific equality and reuse after a later failure. It does not
prove the staged top-level `APA_CORE` orchestration, multi-sample grouping, or
portable nf-core module compliance; those remain separate integration gates.
