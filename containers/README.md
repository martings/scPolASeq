# Reproducible execution environments

This directory is the source of the scPolASeq runtime images. The layout follows
the per-tool Docker/Apptainer recipes and validation guide used in
`martings/nnclinssoap`, adapted to the five environments declared in `envs/`.

## Evidence from Deepthought

The documented May 16–17 Nextflow session used
`-profile deepthought,apptainer,dengue_day1_wk2_deepthought`, the top-level
shared SIF cache selected by the Deepthought profile, and the same SCAPTURE
stage that this pipeline uses.
The completed workflow reported 99 succeeded, 9 cached, and 0 failed tasks. The
trace included completed Sierra and SCAPTURE tasks. Direct import/version probes
also passed for the top-level SCAPTURE SIF (Python 3.8.19, NumPy 1.19.5,
pandas 1.1.5, TensorFlow 2.4.1). This is evidence that Apptainer can run the
SCAPTURE workload on Deepthought, including its TensorFlow/DeepPASS stack.

The older nested per-project cache is not equivalent: its SCAPTURE image has a
NumPy/pandas binary ABI failure. Keep using the profile-selected cache until
replacement images pass validation.

## Build and run

The tracked Conda YAML files are the environment specifications referenced by
Nextflow's `conda` profile; that profile has not passed end-to-end validation.
In particular, the historical SCAPTURE lockfile has a known NumPy/pandas ABI
incompatibility. The `.def` files are the Apptainer build recipes.
Build the Apptainer images in a writable local/cache directory, then validate
them before using them in a run:

```bash
bash containers/build.sh --cache-dir /path/to/sif-cache
```

The builder should retain existing SIFs unless `--force` is specified. SCAPTURE
is an opt-in Nextflow module (`--enable_scapture true`) and its image must be
present when that module is enabled. On Deepthought, use
`-profile deepthought,apptainer`; the profile selects the cache used by the
documented completed run. On other hosts, set `--apptainer_cache_dir` to the
directory containing the SIF images.

For Conda execution, use `-profile conda` (and the relevant machine profile).
The Conda environment files are module-specific so Nextflow creates only the
environments required by scheduled processes. SCAPTURE's historical lockfile
is preserved as provenance. The checked-in `scpolaseq-scapture.def` is not yet
proven to reproduce the working SIF: it clones SCAPTURE without a commit pin
and applies a text-anchored patch. The old Python 3.7 lock also has an
incompatible NumPy/pandas ABI. The smoke checks below validate the existing
working SIFs; they do not validate a fresh rebuild from the current recipes.

## Docker and Nextflow status

`conf/docker.config` currently only turns on Docker. Module `container`
directives currently resolve to local SIF paths when
`--apptainer_cache_dir` is set, and to no image otherwise. Therefore the Docker
profile does **not yet select runnable OCI image names**. A Docker recipe alone
would not make `-profile docker` operational; the module image mapping and
container build/push names must be added together, then tested with a small
Nextflow run. Do not treat the current Docker profile as validated.

## uv assessment

`uv` is a good fit for locking and installing the ordinary Python-only
environment (for example the `python` module), and can produce a reproducible
`uv.lock` for a Docker/Apptainer build. It is not a replacement for Conda across
this workflow: the pipeline also needs Bioconductor/R, STAR/samtools, and
SCAPTURE's legacy Python 3.8 + TensorFlow 2.4/CUDA stack. Keep those in Conda or
the proven SCAPTURE image. A uv pilot should target a Python-only module, compare
its imports and command-line behavior against the existing Conda/SIF image, and
only then become a build input. No uv lock is checked in yet because the
repository's current Python modules share dependencies supplied by Conda and a
partial pip lock would give a misleading reproducibility claim.

## Release gate

Before publishing new Docker, Apptainer, or Conda recipes:

1. Pin every Git-sourced tool to a commit and record base-image digests.
2. Build into a clean cache, preserving build logs and image digests.
3. Run each image's import/version smoke checks, including Sierra and
   SCAPTURE/DeepPASS.
4. Run the synthetic Nextflow profile through Apptainer and Docker/Conda as
   applicable.
5. Run the PBMC1k end-to-end profile and verify expected task completion,
   output manifests, and APA result tables before publication.

The current recipes have not yet passed that complete release gate.
