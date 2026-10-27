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
matched the historical catalog exactly. This checks catalog publication and
conversion in the integrated branch; it does not exercise a complete staged
PBMC1k PAScall/DeepPASS run or publication from every later process.

The separate real PASquant gate, including exact matrix parity and cached
`-resume`, is documented in `docs/scapture_split_draft.md`. The remaining
release gate is a complete staged PBMC1k execution through `APA_CORE`, with
output and manifest comparison and a downstream failure/resume check after
DeepPASS. Until that gate, this mode remains opt-in and this PR remains draft.
