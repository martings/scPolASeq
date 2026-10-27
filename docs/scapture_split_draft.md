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

The process modules have `meta.yml` contracts and `versions.yml` output. In the
container, SCAPTURE's Git commit is recorded when available, alongside a script
SHA256 fingerprint. The current container recipe clones a moving upstream HEAD;
pin that commit before proposing these as portable nf-core modules. Their
`meta.library_id` key is also a local pipeline contract; portable nf-core
modules should use `meta.id` or `task.ext.prefix`.

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
preparation stage. Annotation and PAScall still require a real, bounded run.
