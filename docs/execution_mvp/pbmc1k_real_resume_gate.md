# PBMC1k real execution and resume gate

This runbook records a **real** PBMC1k Nextflow session and the checks required
before accepting changes to coverage, Sierra, or SCAPTURE/DeepPASS. A stub can
check channel wiring, but it cannot verify BigWig conversion, DeepPASS inference,
PASquant, container dependencies, or meaningful output. See the [execution
dossier](pbmc1k_execution_2026-10-24_v1.md) for the original measurements.

## Frozen baseline (observed, not a passing release gate)

| Item | Value |
|---|---|
| Git source commit | `7eb27942c446ea3d962f03f7961bc7c6e0e1ca99` |
| Nextflow revision reported by log | `6b4f463f817f00f9e528c8e437e180c5` (do not treat as a Git commit) |
| Session ID for `-resume` | `c3d2f281-2d73-4d34-99c6-cc876610548a` |
| Initial / resumed run | `evil_lovelace` / `romantic_einstein` |
| Work directory | `/scratch/work/scPolASeq_runs/pbmc1k_20261024_replicate` |
| Historical output | `/scratch/results/scPolASeq/pbmc1k_20261024_replicate` |
| Nextflow history directory | `/home/megas/agents/.worktrees/scpolaseq-pbmc-fixes/.nextflow` |
| Resumed terminal status | 2026-10-25 11:09:59 UTC: 24 completed, 14 cached, 0 failed |

The historical run completed both APA routes. Its nine Sierra TSVs were nonempty;
SCAPTURE published 21,906 site-catalog data rows and 101,398 scored peaks, and
PASquant published a nonempty `.KeepCell.UMIs.tsv.gz` matrix. These are comparison
figures, **not biological acceptance thresholds**. All 18 published BigWigs were
zero bytes while paired bedGraphs had data. Thus the historical result fails the
coverage gate despite Nextflow reporting `OK`.

The task hash identifies a cache/work entry; it is **not** an argument to
`-resume`. Nextflow resumes a session ID, then decides per task whether its hash
and outputs still match. Baseline entries worth checking are:

| Process | Historical task hash | Reason to inspect |
|---|---|---|
| `SCAPTURE_FILTER` | `85/5407a9` (full work suffix `5407a938ac61cb168479ca6cc8a528`) | Original single-task boundary; a split should expose new phase hashes. |
| `PAS_REFERENCE_BUILD` | `81/713111` | Stable upstream comparison. |
| `SIERRA_QUANT` | `7e/ace907` (one of nine) | Check real Sierra outputs, not exit status alone. |
| `COVERAGE_TRACKS` | `7b/cd0c46` (one of nine) | Historical cached BigWig is invalid; corrected conversion must execute. |

From the **historical launch directory**, retrieve all task entries with:

```bash
cd /home/megas/agents/.worktrees/scpolaseq-pbmc-fixes
nextflow log c3d2f281-2d73-4d34-99c6-cc876610548a \
  -f 'hash,process,status,exit,workdir'
```

Keep this worktree, its `.nextflow` history, original work directory, and
historical outputs intact. It currently has unrelated uncommitted container
files; do not reset or check it out to another commit.

## Preconditions for a new real run

1. Use an isolated candidate source checkout and record `git rev-parse HEAD`,
   `git status --short`, `nextflow -version`, `apptainer --version`, and the
   exact command. Recheck that the input sample sheet and reference paths from
   `test_pbmc1k_deepthought` still exist; retain their hashes privately.
2. Build `scpolaseq-python.sif` from the candidate recipe in a **new** staging
   directory. Run `bash containers/validate.sh STAGING_DIR --env python` there.
   That check executes `bedGraphToBigWig` and verifies a real BigWig header.
   Record `sha256sum` for all five SIFs. Populate the staging directory with
   the other four validated images; do not overwrite the shared cache.
3. Check scratch capacity before starting. On 2026-10-27, `/scratch` had only
   88 GB available. A code or image change may invalidate enough task hashes
   to require substantially more space than a small resume. Estimate from
   historical work sizes and stop if there is no safe margin.
4. Choose a new output directory and preserve the original outputs. Use the
   same historical work directory **only** for the controlled resume; expect
   selective recomputation when source, parameters, inputs, or containers
   change. A clean run in a new work directory is a separate later gate.

The historical Python SIF in `/scratch/cache/apptainer` has SHA-256
`cc0010d9f80e7b8ee60569a7bf690482c38d5f1334f8cf5050020933efa841f7`
and lacks `bedGraphToBigWig`. It cannot pass this gate. A corrected image was
built and checked in isolated staging on 2026-10-27, as recorded below. The
other four images have not been staged with it. **No candidate resume has yet
run.**

## Focused real-data converter probe — 2026-10-27

`containers/build.sh --env python` built a new image at
`/scratch/tmp/scpolaseq-pr19-python-20261027/scpolaseq-python.sif` with private
Apptainer build cache and temporary directory. Its SHA-256 is
`37a1627ccad4b71f06cf28f443ad15457635476417b5ce1bf2378cdfc53f5ee1`.
`containers/validate.sh ... --env python` passed its synthetic conversion and
BigWig header check. The built image reports NumPy 2.4.6, pandas 2.2.3,
and pysam 0.22.1; keep this exact image digest with the run record.

Next, the image converted both **actual** PBMC1k cluster 5 bedGraphs from
historical coverage task `7b/cd0c46` using that task's chromosome sizes.
The isolated outputs were 1,157,896 bytes (forward) and 1,044,833 bytes
(reverse), both nonempty; both files start with `26 fc 8f 88`, the BigWig
magic. Their SHA-256 values
were `a34eb2d0a00a38eb629895e63ae0a3744ad14d67dc4a59e750026fc48722f3b8`
and `b2553229492dd02c63e73bc1843d2ec96ef1bf8387cd09573d12e8536f0b5ebc`.
The converted files remain under the isolated staging directory, separate
from the historical outputs. This probe validates the real converter on one
group; it does not test Nextflow's new coverage process or all 18 tracks.

## Candidate resume command (execute only after preconditions)

Run from the historical launch directory so Nextflow sees its original
`.nextflow` history. Replace `CANDIDATE` and `STAGING_DIR` with recorded absolute
paths from the isolated checkout and image build. Confirm the profile and paths
with `nextflow config` before launching; the command below has not been tried.

```bash
cd /home/megas/agents/.worktrees/scpolaseq-pbmc-fixes
nextflow run CANDIDATE/main.nf \
  -profile deepthought,apptainer,test_pbmc1k_deepthought \
  --enable_scapture true \
  --outdir /scratch/results/scPolASeq/pbmc1k_candidate_validation \
  --apptainer_cache_dir STAGING_DIR \
  -work-dir /scratch/work/scPolASeq_runs/pbmc1k_20261024_replicate \
  -resume c3d2f281-2d73-4d34-99c6-cc876610548a \
  -ansi-log false
```

Changing the code, image, output path, or process boundary can invalidate
cache entries. Do not claim that a task was rerun or reused solely from its
historical hash; capture the new trace/log. The candidate source and launch
directory are deliberately distinct, so confirm `projectDir`-based inputs and
Apptainer binds before committing compute time.

## Acceptance evidence to record after execution

Record the new session ID, run name, source SHA, five image SHA-256 values,
command, trace, report, process hashes, exit codes, and actual output files in
a dated companion report. Verify each point below using the **candidate output
directory**, never the historical one:

- Nextflow status `OK`, zero failed/ignored tasks. Trace must show which tasks
  were `CACHED` and which executed; changes to coverage and SCAPTURE phases must
  have the expected new hashes. Check the corresponding `.command.sh`, `.exitcode`,
  and `.command.log` if any phase is surprising.
- Each emitted BigWig has nonzero size and a valid BigWig header, and its
  paired bedGraph has data. The corrected `COVERAGE_TRACKS` tasks must not reuse
  the 2026-10-25 zero-byte BigWigs. Run a real reader such as `bigWigInfo` or
  `pyBigWig` on representative forward and reverse tracks.
- All nine Sierra groups have parseable, nonempty result TSVs. Compare group IDs
  to the grouping manifest and report row counts; a changed count needs review.
- SCAPTURE/DeepPASS has parseable evaluated peak tables, positive and negative
  predictions where expected, a valid 21,906-row baseline catalog comparison,
  and a parseable PASquant `.KeepCell.UMIs.tsv.gz` matrix with cells and nonzero
  signal. The 101,398 scored-peak baseline is a diagnostic comparison. Confirm
  the GPU selection in DeepPASS logs without inferring utilization from device
  visibility alone.
- For a split SCAPTURE subworkflow, inspect a task hash/work directory for each
  actual phase: preparation, annotation, PAScall/DeepPASS, merge, and PASquant.
  Deliberately failing a late phase in a **separate small real fixture run**
  should allow earlier successful phases to be cached on resume while the late
  phase is retried. The production run must fail if PASquant fails; no `|| true`,
  fabricated empty outputs, or unreported optional matrix.
- Store representative hashes and row counts for downstream comparison. Treat
  a changed biological result as a review item, not an automatic failure based
  on an unapproved numerical threshold.

## Issue order and release decision

1. **P0: #6 / PR #19, coverage and image gate.** Stage and validate the fixed
   Python SIF, then pass real BigWig conversion and candidate resume. Keep PR
   draft until this evidence exists.
2. **P0: #16 and #25–29, SCAPTURE phases.** Make phase boundaries resumable,
   remove silent PASquant success, add real output contracts and provenance,
   then test with a small real fixture and PBMC1k resume.
3. **P1: Sierra contract and runtime parity.** Pass explicit group metadata,
   pin Sierra/SCAPTURE source revisions, align Conda and SIF behavior, and
   document module inputs/outputs and versions. Keep portability claims
   limited to tested profiles.
4. **P1: downstream and larger datasets.** Connect only outputs that pass the
   above gates; defer more expensive runs until image, space, and phase checks
   are recorded.

This is an integration gate for this repository. The [nf-core module
guidelines](https://nf-co.re/docs/guidelines/components/modules) and
[subworkflow testing specification](https://nf-co.re/docs/specifications/components/subworkflows/testing)
also call for stable input/output contracts, tool versions, metadata, and
nf-test coverage. Their minimum stub test does not replace this PBMC1k real
execution gate. Nextflow's [execution management
guide](https://training.nextflow.io/latest/nextflow_run/03_manage_executions/)
explains the session-based cache used by `-resume`.
