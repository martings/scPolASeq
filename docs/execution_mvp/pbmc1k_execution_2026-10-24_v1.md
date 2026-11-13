# PBMC1k execution dossier — 2026-10-24/25 (v1)

## Assessment

The resumed PBMC1k Nextflow invocation reached `Status: OK` on 2026-10-25 at 11:09:59 UTC. It reports 24 completed and 14 cached tasks, with no failed, retried, or aborted tasks. This is evidence of technical completion for the recorded workflow tasks. It is **not a scientifically accepted reference run**: all 18 published BigWig files were empty, although their matching bedGraph files contained data. The workflow version used returned success after silently substituting empty files when BigWig conversion was unavailable or failed. PR #19 adds a fail-closed fix and regression coverage; this run predates that fix and must not be represented as validating it.

Sierra, SCAPTURE/DeepPASS, APA feature calling, and the random-forest scoring path have run and produced nonempty outputs. These results are execution evidence only. There is no independent ground truth review in the evidence collected for this dossier. Human scientific review is pending.

## Identity and provenance

| Field | Observed value | Confidence / qualification |
|---|---|---|
| Dataset | PBMC1k test dataset | High; run parameters and output records identify this dataset. |
| Workflow repository | `martings/scPolASeq` | High. |
| Source commit | `7eb27942c446ea3d962f03f7961bc7c6e0e1ca99` | High; checkout identity in run evidence. |
| Nextflow script revision | `6b4f463f81` | Observed, but not independently mapped to a full Git SHA. |
| Nextflow | 25.10.4 | High. |
| Apptainer | 1.4.5 | High. |
| Profiles | `deepthought,apptainer,test_pbmc1k_deepthought` | High. The named Deepthought profile is host-specific. |
| Session | `c3d2f281-2d73-4d34-99c6-cc876610548a` | High; both invocations share this session. |
| First invocation | `evil_lovelace`, 2026-10-24 20:51 UTC, no `-resume` | Launch evidence present; its terminal cause/status was not established. |
| Resumed invocation | `romantic_einstein`, 2026-10-25 06:14 UTC, `-resume` | High; terminal summary below. |
| Terminal status | `OK`, 2026-10-25 11:09:59 UTC | High for the resumed invocation. |
| Elapsed time | 4h 55m 28s for the resumed invocation | High; does not include the first invocation or cached wall time from earlier work. |
| Work / output locations | Deliberately omitted from this public-safe dossier | The run used private Deepthought and scratch locations. |
| Run provenance type | Human-operated | No controller-authenticated observer evidence is asserted. |

The two invocations are reported separately because the first invocation's terminal state is unknown. The completed resumed invocation reports 24 `COMPLETED`, 14 `CACHED`, 0 `FAILED`, 0 `RETRIED`, and 0 `ABORTED`; peak concurrency was 4 tasks, peak requested CPU 16, and peak memory 80 GB. Task-duration totals include cached work and should not be interpreted as wall-clock elapsed time. The minimum observed free scratch space was approximately 87 GB; after completion it was approximately 94 GB free (about 92%).

## Execution and output review

| Output family | Execution evidence | Output evidence and interpretation |
|---|---|---|
| STARsolo alignment | Completed. STARsolo reported 66,601,887 input reads; 87.74% uniquely mapped and 95.44% unique plus multiple mapping; mismatch rate 0.55%. | Nonempty BAM/BAI and filtered matrix were present. Filtered MEX: 62,700 features × 1,139 cells, 2,633,085 nonzero entries. Median 7,109 UMIs/cell and 2,098 genes/cell. These are observed QC values, not independently benchmarked acceptance thresholds. |
| Barcode identity | Registry, filtered matrix, and annotation IDs were checked. | 1,139 cell barcodes agree across the inspected artifacts. |
| Barcode filter | Completed. | Its `total_reads` is 73,849,140 BAM alignment records, with 67,889,904 valid and filter rate 0.9193. This count is not comparable to STAR's input-read count. |
| Clustering / labels | K-means completed. | Eight clusters and two unassigned groups; cell types are `unlabeled`. Condition and batch labels were absent. |
| Sierra | Nine tasks completed (clusters 1–8 and unassigned). | Nonempty TSV outputs were observed for the nine task groups. This confirms execution and files, not biological correctness. |
| APA feature/calling | `APA_CALLING` completed with exit 0. | `apa_events.tsv`: 147,568 data rows; `scored_apa_events.tsv`: 7,053,948 rows; 30,129 predictions at threshold 0.5. APA feature table has 7,053,948 rows. |
| Random-forest model | Training completed on 7,053,948 samples and 8 features; 6,651 positive labels. | Cross-validation AUC 0.7613 ± 0.0055; training accuracy 0.9966. Accuracy is misleading under this class imbalance; no independent validation set or scientific acceptance was recorded. |
| PAS reference | Reference build completed. | 783,772 catalog rows before merge; 741,392 merged PAS references. Manifest records 214,767 terminal exon rows and 569,005 known polyA rows. |
| SCAPTURE / DeepPASS | `SCAPTURE_FILTER` completed exit 0, including PAScall peak calling, annotation, DeepPASS evaluation, PASmerge, and PASquant. | 101,398 scored peaks: 28,464 exonic, 70,568 intronic, 2,366 extended. Site catalog: 21,906 data rows. `callpas_jobs` is `NA` and `site_catalog_rows` is blank in the manifest: provenance fields need correction. |
| Coverage tracks | Coverage tasks reported exit 0. | All 18 published `.bw` outputs were 0 bytes, while paired `.bedGraph` files had data. These BigWigs are invalid and must not be consumed. The old production module used `touch` if `bedGraphToBigWig` was absent and `|| touch` if conversion failed. |
| GPU / DeepPASS | Three DeepPASS logs report TensorFlow GPU:0 on an RTX 3080 with 8,935 MB visible and “Finished Prediction”. | GPU was exposed and detected; no `nvidia-smi` sample covers the short inference window, so utilization cannot be established. Logs also report `ptxas` exit code 256 and `libcublas.so.11` undefined symbol `cublasGetSmCountTarget`; prediction finished despite these warnings. Later PASquant and random-forest stages used CPU; this is expected for the recorded stages. |
| Reports | Nextflow trace, report, timeline, and DAG were published. | Fontconfig reported no writable cache directories during SCAPTURE/PASquant; those tasks still exited 0. |

## Runtime artifacts and observed hashes

The SIF hashes below were observed from the run environment. They identify the bytes present there; they are not independently expected/release-approved hashes.

| Image role | SHA-256 |
|---|---|
| Alignment | `b7ecc90e5e1fd17cf6436f2f4952fcbcb7f6e7cdcc7eb3fc38519277b3ab911e` |
| Clustering | `d16d4fa79f29820f11ccd6003255d1f51f847a9652a63ea811b075e8386d45e1` |
| Python | `cc0010d9f80e7b8ee60569a7bf690482c38d5f1334f8cf5050020933efa841f7` |
| R | `78efbf3b49ba4490733ff7b9417a99336854e5f9f3e0272d370b3b2a2a41bdc2` |
| SCAPTURE | `fb087f77f3c19083008aa83086669f417426a43b2c52d4901d0732d89aa13a63` |

Selected reference hashes observed in the run environment:

| Reference artifact | SHA-256 |
|---|---|
| `reference.genome.fa` | `e49b92b3e4f321bf254c042f25b726d9931c4d74c7523e8b6bb530e63b0cfd4b` |
| `reference.annotation.gtf` | `46fc6e733fc73b236ffa68b770c5ee0cedd998797107cb4691bd8ca5c104df16` |
| `known_polya.reference.tsv` | `3466af7496d343aed4dcec704790c91f7bbb045e18da0005b1c520a2a7d1a6bc` |
| STAR index `Genome` | `be7ab7de7a9552474a7ae8128326b302a66762ab1295e1ac1970e4704874b0fe` |
| STAR index `SA` | `d49e7d547a7c13bc7987d31f45e9912e4cad4bd216379124bb5a12ecb24046f3` |
| STAR index `SAindex` | `dda81433f64d620b2d043f9f3e2b4d50ac59b4297f21ac347fdad612398da02f` |
| STAR index `Genome` parameters | `096e8df336f42657d5e629b71ee0aff3ca6f827bc68bf589dcf7e07971ceb772` |

Reference index metadata reports STAR 2.7.4a. The PR #19 image recipe uses STAR 2.7.11b, so reuse of this index across those versions remains unvalidated. No FASTQ content was read or hashed. The input manifest records two validated rows, but no independent input content hashes were available. The reference files have observed hashes but no independently expected hashes.

## Historical comparison and confidence

| Evidence period | What is established | What is not established | Confidence |
|---|---|---|---|
| PBMC1k, May 7, 2026 | A historical Deepthought Apptainer run reported 36 successful tasks and nonempty STARsolo alignment outputs. | Independent scientific QC, complete Sierra + SCAPTURE provenance, and equivalence to this October run. | High for technical completion as recorded; low for scientific acceptance. |
| PBMC1k, May 9, 2026 | Historical Sierra and SCAPTURE result files exist. | A single complete execution log binding both result families to one run. File existence/timestamps do not prove common provenance. | High for historical files; low for run-level linkage. |
| Dengue, May 16–17, 2026 | Separate run reported 99 successful and 9 cached tasks, including 35 Sierra and 4 SCAPTURE tasks. | Any PBMC1k validation. | High as separate workflow capability evidence; not PBMC1k evidence. |
| PBMC1k, Oct 24–25, 2026 | One Nextflow session with an initial invocation and a resumed invocation; resumed invocation completed and produced both Sierra and SCAPTURE/DeepPASS outputs. | First invocation terminal cause; hashes of FASTQ contents; scientific acceptance; BigWig validity; GPU utilization during inference; controller-authenticated provenance. | High for resumed technical completion and observed outputs; partial for full-run provenance; no scientific sign-off. |

Similar file or row counts across historical and current runs would not demonstrate scientific equivalence. Input identity, exact parameters, reference identity, software versions, and independent QC acceptance would be needed.

## Repository review and disposition of PR #19

The generic `test_pbmc1k` configuration does not hard-code Deepthought host paths. Remaining `/scratch` references belong to explicitly named Deepthought profiles and host runbooks. Legacy sample sheets in `tests/data` are referenced by host-specific Deepthought profiles; configurable helper/notebook defaults also retain `/cluster` paths and should remain configurable. This audit found no reason to remove the known-good Apptainer profile.

PR #19 preserves the existing Deepthought Apptainer execution route. The run used `7eb27942c446ea3d962f03f7961bc7c6e0e1ca99`, not the PR branch, and therefore does not test the PR changes. Existing PR checks had passed on the earlier branch head. This audit adds a fail-closed coverage converter dependency and negative/positive output tests. The Python Apptainer image recipe now pins `ucsc-bedgraphtobigwig=482`, and image validation checks the executable. The current run's Python SIF does not contain a verified working converter; a newly built PR image has not yet been installed into the shared Deepthought cache.

The old STAR index version mismatch remains an explicit compatibility question. `--nv` was present in historical launch configuration but does not establish that workflow resources or every task require a GPU. No workflow pin is promoted by this dossier; scientific review is still required before changing the Megas Infra pin.

## Missing evidence and next reproducibility steps

1. Preserve the exact samplesheet and parameter file under controlled storage and publish a sanitized copy plus SHA-256; do not publish private paths or FASTQ hashes.
2. Map the observed Nextflow revision token to a Git object and capture a source tree fingerprint.
3. Record independent expected image/reference hashes and retain image build metadata.
4. Fix SCAPTURE manifest fields for `callpas_jobs` and `site_catalog_rows`.
5. Build the updated Python SIF in a temporary location, verify `bedGraphToBigWig`, then run the focused real coverage conversion integration test before populating a shared cache.
6. Re-run PBMC1k only under separate authorization after image and scientific-review gates; do not count the October outputs as BigWig-valid.
7. Have a domain reviewer assess alignment/QC, APA model validity, Sierra, and SCAPTURE/DeepPASS outputs before calling this a scientific reference case.

The machine-readable companion is `pbmc1k_execution_2026-10-24_v1.json`. Its field names are a proposal for handoff, not an approved Controller/API schema.
