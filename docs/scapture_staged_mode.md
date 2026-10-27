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
