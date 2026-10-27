# Sierra quantification contract (draft, issue #32)

`SIERRA_QUANT` accepts one explicit `(meta, group_level, group_id, grouped_bam,
grouped_bai, pas_reference, cell_annotations, gtf)` tuple plus the staged
`bin/sierra_quant.R` script. Its output keeps the same group-scoped TSV,
manifest and log names and adds a group-scoped `versions.yml`. Inputs are
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
paths or shell syntax. The adapter rejects duplicate group IDs/BAM filenames
and requires the manifest to cover the exact BAM/index set.

The R wrapper records the Sierra and R versions before quantification. The
module checks nonempty staged inputs and requires both its TSV and version
file. The group-scoped version filename prevents concurrent tasks from
publishing to one `versions.yml` path. A TSV with only a header remains a
possible scientific outcome for a group without matching cells/reads; this
contract does not treat that as a quantification success claim.

## Validation

On 2026-10-27, the `full` Nextflow stub scenario passed with 23 processes.
A separate real `SIERRA_QUANT` task ran the existing PBMC1k `cluster_5` BAM,
index, PAS reference, annotations and GTF through the existing R SIF. It
completed with 9,764 TSV lines; SHA-256 of its TSV and the October output was
identical: `329fc77cf5e737fe0ea9ed1a9ee9bc8d0bb76b898db6601b7ac28dbbd029f6c4`.
The version output reports Sierra 0.99.27 and R 4.3.3. This is one group,
not a new nine-group run or scientific validation. The SIF is the historical
unversioned local image; pinning/build reproducibility is tracked in #31.

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
