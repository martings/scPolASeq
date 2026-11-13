# Sierra quantification contract (issue #32)

`SIERRA_QUANT` (`subworkflows/local/sierra_quant/main.nf`) composes an
`nf-core`-compliant atomic module `SIERRA_COUNTPEAKS`
(`modules/nf-core/sierra/countpeaks/main.nf`) between two thin `scPolASeq`
local adapters:

1. `SIERRA_QUANT_PREPARE` (`modules/local/apa/sierra_quant_prepare/main.nf`) — filters
   `cell_annotations` for the scoped `(library_id, group_level, group_id)` barcode
   whitelist and converts `pas_reference.tsv` into Sierra `peak_sites.tsv` +
   canonical `site_id_map.tsv`.
2. `SIERRA_COUNTPEAKS` (`modules/nf-core/sierra/countpeaks/main.nf`) — reusable atomic
   module wrapping `Sierra::CountPeaks(..., strandMode = ...)` over standard
   `meta` + BAM/BAI/whitelist + peak sites BED/TSV + GTF inputs, emitting standard
   10x MEX (`matrix.mtx.gz`, `barcodes.tsv.gz`, `sitenames.tsv.gz`), log, and
   `versions.yml`. Strandedness (`forward` -> `strandMode=1`, `reverse` ->
   `strandMode=2`, `unstranded` -> `strandMode=0`) is configured via `ext.args`
   (`conf/modules.config`).
3. `SIERRA_QUANT_TABLE` (`modules/local/apa/sierra_quant_table/main.nf`) — converts the
   MEX triplet + `site_id_map.tsv` back into the group-scoped
   `*.sierra_quant.tsv`, `*.sierra_quant.manifest.tsv`, `*.sierra_quant.log`,
   and `versions.yml` contract consumed by `APA_CORE`.

Because `SIERRA_QUANT` fans out each `(group_level, group_id)` group into three
process stages (`PREPARE` -> `COUNTPEAKS` -> `TABLE`), full pipeline runs with 2
default Sierra groups execute 4 additional tasks (+4 processes) compared to the
legacy monolithic `SIERRA_QUANT` process.

The `SIERRA_QUANT` subworkflow still accepts one explicit `(meta, group_level, group_id, grouped_bam, grouped_bai, pas_reference, cell_annotations, gtf)` tuple. Inputs are
Nextflow `path` values; the module does not construct dataset paths. The
adapter in `APA_CORE` obtains `group_id` and BAM filename from the grouping
manifest and checks that every manifest row has a BAM and index. It no longer
infers group identity from filename prefixes. The manifest produced for the
October PBMC1k cluster run contains eight clusters and `unassigned`, so this
adapter preserves nine Sierra tasks for that profile. An empty cell-type group
still yields no Sierra task.

The original `group_id` remains the value used to filter cells and is recorded
in the manifest. File names use a sanitized group token matching the grouped
BAM writer's allowed characters. Shell arguments quote the original value;
labels containing spaces, `/`, `+`, or apostrophes therefore do not become
paths or shell syntax. The BAM writer fails before opening files when two
labels sanitize to the same name. Input validation rejects distinct library
IDs that sanitize to the same filename token. The adapter rejects duplicate
group IDs/BAM filenames and requires the manifest to cover the exact
BAM/index set.

The R wrapper records the Sierra and R versions before quantification. The
module checks nonempty staged inputs and requires both its TSV and version
file. The group-scoped version filename prevents concurrent tasks from
publishing to one `versions.yml` path. A TSV with only a header remains a
possible scientific outcome for a group without matching cells/reads; this
contract does not treat that as a quantification success claim. A corrupt BAM,
`scanBam` or `CountPeaks` error, or a missing matrix after matching barcodes
fails the task. An empty scoped whitelist or zero barcode overlap produces an
explicitly logged header-only TSV.

## Phase B2 validation (nf-core module + local adapter split)

- Real Nextflow `tests/modules/apa/test_sierra_quant.nf` run (`-profile deepthought,apptainer`, `--strandedness reverse`) on `Adult_Basal.cluster_1.grouped.bam` executed `SIERRA_QUANT_PREPARE` -> `SIERRA_COUNTPEAKS` (`ghcr.io/martings/sierra:84b71c4` / `scpolaseq-r.sif`) -> `SIERRA_QUANT_TABLE` and reproduced exact biological parity: `rows=28122`, `total_umis=48850`, `Krt14=51`, `Acta2=114`.
- `modules/nf-core/sierra/countpeaks/tests/main.nf.test` provides `nf-test` coverage for `SIERRA_COUNTPEAKS`.

## Validation (historical monolith baseline)

On 2026-10-27, the `full` Nextflow stub scenario passed with 23 processes.
A separate real `SIERRA_QUANT` task ran the existing PBMC1k `cluster_5` BAM,
index, PAS reference, annotations and GTF through the existing R SIF. It
completed with 9,764 TSV lines; SHA-256 of its TSV and the October output was
identical: `329fc77cf5e737fe0ea9ed1a9ee9bc8d0bb76b898db6601b7ac28dbbd029f6c4`.
The version output reports Sierra and R runtime versions (now pinned to `martings/Sierra` with native `strandMode` support via issue #31; `bin/sierra_quant.R` delegates directly to `Sierra::CountPeaks(..., strandMode = strand_mode)`).

After review fixes, the real `cluster_5` task ran again at hash `59/a19585`
in session `1b21f9d5-e12e-4615-93c3-b05dc718d91c`, retained the exact TSV
SHA256 above, and was cached on the next explicit `-resume`. Two additional
real module invocations with unmatched labels `T cell/A+B` and `T'cell`
completed exit 0, wrote safe output names `T_cell_A_B` and `T_cell`, and kept
the original labels in their manifests. They correctly produced header-only
TSVs because no PBMC1k cells had those labels; they test quoting and naming,
not quantification biology.

The portable module test is `tests/modules/apa/test_sierra_quant.nf`. It
requires explicit file parameters for a real run; no Deepthought paths are
committed to the test or module.

After the fail-closed fixes, the real cluster_5 rerun completed at hash
`01/663c82` in the same session, retained the historical SHA-256, and a
second explicit `-resume` cached `01/663c82`. A nonempty corrupt BAM caused
the real Apptainer module test to exit 1 at hash `ad/30a0db`; its log shows
the `scanBam` failure and no successful quantification was published. The
full local stub suite (23 processes and two module tests) and 12 barcode
identity/regression checks passed. These tests do not replace a fresh
nine-group PBMC1k run.
