# Workflow smoke tests

Expected commands once `nextflow` is available:

```bash
nextflow run . -profile test -params-file tests/config/params_full.json -stub-run
nextflow run . -profile test -params-file tests/config/params_from_bam.json -stub-run
nextflow run . -profile test -params-file tests/config/params_feedback.json -stub-run
```

Recommended non-stub integration test:

```bash
nextflow run . -profile test -params-file tests/config/params_full.json
```
