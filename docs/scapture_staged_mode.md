# Opt-in staged SCAPTURE path

`--enable_scapture true --scapture_execution_mode staged` selects the new
cache-separated SCAPTURE path. The default remains `monolith`, so existing
runs keep their current task graph. Both paths emit the same `APA_CORE`
channels: `scapture_catalog` (library metadata and normalized TSV),
`scapture_quant` (library metadata and gzipped UMI matrix), and
`scapture_manifest` (path to a field/value TSV). The staged manifest includes
the monolith's discovery fields plus PASquant completion and matrix dimensions.
The historical monolith left `site_catalog_rows` blank due to an AWK bug; the
staged manifest records the measured row count.

The staged path runs GTF preparation, annotation, PAScall with DeepPASS,
per-category PASmerge, combined BED12 reference, catalog conversion, and
PASquant as separate tasks. PAScall and PASquant have their own 14 CPU, 64 GB,
48 hour limits. A later stage failure does not force DeepPASS to run again.
Quantification is mandatory in this mode: a missing or invalid matrix fails
the pipeline. An explicit catalog-only mode remains future work in issue #28.
Staged mode currently requires each `scapture_group_id` to contain exactly one
`sample_id`; it fails before SCAPTURE when a group spans samples. The monolith
path retains its existing behavior for such groups until multi-sample barcode
selection is validated in the staged PASquant module.
The `scapture_group_id` used for output filenames must start with a letter or
digit and contain only letters, digits, `.`, `_`, or `-`; unsafe IDs fail before
any SCAPTURE task starts.
The catalog, evaluated and merged BEDs, combined BED12 reference, PASquant
matrix and manifest, metrics, and stage logs are published under
`${params.outdir}/scapture`; `versions.yml` files remain in task workdirs.

On Deepthought, PAScall carries the `gpu` label and requests one accelerator.
The host's broad GPU label config also sets `containerOptions=--nv`, which
would hide the task's DeepPASS predictor bind. The name-specific PAScall
selector keeps that bind, while `apptainer.runOptions` supplies `--nv`.
A bounded real invalid-BAM probe on 27 October 2026 generated
`/scratch/tmp/scapture_orch_gpu_probe_override_work/db/780867b9ffc3dec6a31e401b89e237/.command.run`
with both `--nv` and `--bind DeepPASS_Predict.py:/opt/scapture/DeepPASS/Predict.py`.
That probe verifies the launcher and expected failure propagation; it does not
measure DeepPASS GPU computation.

Use a persistent work directory and explicit Nextflow session UUID when
resuming. For example, after the first run records its UUID:

```bash
nextflow run . -profile deepthought,apptainer,test_pbmc1k_deepthought \
  --enable_scapture true --scapture_execution_mode staged \
  --outdir /scratch/results/scPolASeq/pbmc1k_staged \
  -work-dir /scratch/work/scPolASeq_runs/pbmc1k_staged \
  -resume SESSION_UUID
```

Keep the same inputs, parameters, containers, and work directory for a cache
check. Changed helper scripts are staged process inputs and invalidate the
relevant task hash. `--scapture_resume_pascall_dir` is a monolith-only recovery
option; the staged path uses Nextflow `-resume` instead.

`tests/scapture_staged_contract.nf` checks DSL2 wiring with an empty input
channel. This is a parse test only. Real-data release evidence should compare
evaluated BEDs, the combined reference, normalized catalog, PASquant matrix,
and manifest against the frozen PBMC1k monolith outputs. It should then
deliberately fail a downstream task and demonstrate that PAScall is cached on
an explicit `-resume` run. The historical run and hashes are recorded in
`docs/scapture_split_draft.md`.

## Bounded integration checks (27 October 2026)

The integrated `SCAPTURE_STAGED` graph parsed with an empty input channel, and
the repository's full local stub suite passed (23 pipeline processes, three
module contract runs, and the manifest helper test). On frozen PBMC1k
evaluated BEDs, a real Apptainer `SCAPTURE_CATALOG` task published its TSV
and metrics to a new outdir. The published TSV SHA256
`1621066b0cd864a589b45f7d882ff20cbbae594e5a537bb5adbf00a22d507dc4`
matched the historical catalog exactly. A separate 1,200 bp real annotation
fixture published only `tiny.annotation.log` to its fresh outdir, confirming
the selective log rule without copying annotation intermediates. These checks
cover publication and conversion in the integrated branch; they do not exercise a complete staged
PBMC1k PAScall/DeepPASS run or publication from every later process.

The separate real PASquant gate, including exact matrix parity and cached
`-resume`, is documented in `docs/scapture_split_draft.md`. The remaining
release gate is a complete staged PBMC1k execution through `APA_CORE`, with
output and manifest comparison and a downstream failure/resume check after
DeepPASS. Until that gate, this mode remains opt-in and this PR remains draft.

## Full frozen-input staged gate (28 October 2026)

`tests/scapture_staged_historical.nf` invokes the complete staged graph on the
read-only historical SCAPTURE task's BAM, references, and cell annotations. It
starts with a deliberate failure gate after PASquant so a later explicit resume
can test cache reuse without repeating DeepPASS. The run uses Nextflow 25.10.4,
`deepthought,apptainer`, session UUID
`0ab7d746-6dfe-48d4-afb7-8c680c9ce639`, work directory
`/tmp/scapture_staged_full_20261027_work`, and a separate output directory
`/tmp/scapture_staged_full_20261027_out`.

```bash
nextflow run tests/scapture_staged_historical.nf -profile deepthought,apptainer \
  --scapture_task_dir /scratch/work/scPolASeq_runs/pbmc1k_20261024_replicate/85/5407a938ac61cb168479ca6cc8a528 \
  --scapture_library_id pbmc_1k_v3 --scapture_test_fail_after_quant true \
  --outdir /tmp/scapture_staged_full_20261027_out \
  --apptainer_cache_dir /scratch/cache/apptainer \
  -work-dir /tmp/scapture_staged_full_20261027_work \
  -resume 0ab7d746-6dfe-48d4-afb7-8c680c9ce639
```

The first attempt found that SCAPTURE annotation inherited `TMPDIR=/scratch/tmp`,
which did not exist inside the Apptainer task. Commit `e234013` sets a task-local
TMPDIR. The same explicit `-resume` then cached PREPARE (`61/ad622e`), completed
ANNOTATION (`83/19c80e`), and submitted PAScall (`9a/f6ae8f`). That task was
interrupted externally without a Nextflow error or `.exitcode` while calling
peaks. A persistent session resumed the same UUID on 27 October and completed
PAScall (`68/72cd8d`) in 4h 42s, all three PASmerge tasks, PASreference,
catalog, PASquant, and manifest. The deliberate post-quantification gate exited
42 as designed. Repeating the command with only
`--scapture_test_fail_after_quant false` cached all ten scientific tasks and
finished with exit 0; the final gate ran in a new task (`28/a6923b`).

The published outputs are in `/tmp/scapture_staged_full_20261027_out/scapture`.
The manifest reports `status=complete`, 101,340 valid and scored peaks, 21,904
catalog sites, and a 101,330-site by 1,139-cell PASquant matrix. `gzip -t`
passed; the decompressed matrix has one header and 101,330 data lines with
1,139 barcode columns. The historical monolith used the same 1,139 barcodes,
but its catalog had 21,906 sites and its matrix had 101,386 rows. Catalogs
share 21,774 exact site IDs with identical shared records (130 staged-only,
132 historical-only). Matrices share 101,006 row IDs (324 staged-only, 380
historical-only); 1,564 shared rows have different values. This establishes
end-to-end operation, not byte-for-byte scientific parity.

The PAScall launcher contains `apptainer exec --nv`. Each DeepPASS log reports
TensorFlow creating `/device:GPU:0` on the NVIDIA GeForce RTX 3080, loading
CUDA/cuDNN libraries and completing prediction. The logs also contain a
nonfatal `cublasGetSmCountTarget` symbol warning; all three predictions and
PAScall returned exit 0. This gate validates the staged SCAPTURE subworkflow
on frozen PBMC1k inputs. An integrated `APA_CORE` run remains a separate
release gate.

### Investigation of the row differences

Both PAScall tasks use the same resolved BAM and reference, the six generated
SCAPTURE annotation files are byte-identical, and PAScall uses the same read
length (90), peak width (400), species (human), and thread count (14). The
first divergence is inside SCAPTURE's own peak caller, before DeepPASS:
all 25,931 shared per-gene `exonic.peak.isoform.bed` files are byte-identical,
but 2,502 `exonic.peak.isoform.normality.bed` files differ. The raw global
exonic and intronic peak BEDs then contain 36,639 versus 36,636 and 75,101
versus 75,037 rows, respectively (historical versus staged).

In the SIF, `/opt/scapture/scripts/Normality.R` samples read coverage with
R's `sample()` for the Anderson-Darling, symmetry, and dip tests (lines
39–57), without `set.seed()`. `/opt/scapture/scapture_callpeak` uses these
p-values to keep or reject exonic peaks at `0.0001` (line 150). This is a
directly observed source of run-to-run variation, even with identical inputs
and parameters. All 704 matrix row IDs present in only one run belong to
genes whose normality-filtered BED changed. Of the 262 catalog site IDs
present in only one run, 261 belong to such genes; the last is the same
genomic interval assigned to a different overlapping gene by the catalog.

The two matrices have exactly the same 1,139 barcodes. Their total UMI counts
are 9,079,774 historical and 9,089,259 staged (0.10% difference). Among
shared row IDs, 1,564 rows differ in at least one cell: 1,450 belong to genes
with a changed normality-filtered BED, while 114 belong to 102 genes whose own
combined reference BED is unchanged. The latter may reflect assignments
influenced by changed neighboring sites, but that mechanism has not been
isolated. A separate bounded PASquant test on identical frozen historical
inputs produced byte-identical matrix output, so this full-run comparison
does not implicate the split PASquant module itself. Exact reproducibility
would require a deliberate, documented random-seed policy in SCAPTURE and a
fresh comparison; changing that policy could change peak calls and must be
treated as a scientific-method change.
