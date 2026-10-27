# SCAPTURE discovery split: status and real-data gate

This branch introduces `SCAPTURE_DISCOVERY` as a separate DSL2 subworkflow. It
is not wired into `APA_CORE` yet. `SCAPTURE_FILTER` remains the production path
until PASmerge, PASquant, catalog conversion, and a real `-resume` equivalence
test are complete. A stub run cannot establish this equivalence.

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

Before enabling the subworkflow in `APA_CORE`, run the split on the same frozen
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
