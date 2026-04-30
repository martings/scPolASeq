/*
 * tests/modules/reference/test_prepare_polya_references.nf
 *
 * Standalone module-level test for PREPARE_POLYA_REFERENCES.
 *
 * Purpose:
 *   Validates the output contract of prepare_polya_references/main.nf in stub
 *   mode — verifying that both output files are produced with the expected
 *   columns and manifest fields.
 *
 * Usage (stub — fast, no real tools needed):
 *   nextflow run tests/modules/reference/test_prepare_polya_references.nf -stub-run
 *
 * Exit 0 = all assertions passed.
 */

nextflow.enable.dsl = 2

include { PREPARE_POLYA_REFERENCES } from '../../../modules/local/reference/prepare_polya_references/main'

// ── Minimal params required by publishDir / module script ────────────────────
params.outdir               = 'tests/results/module_prepare_polya_references'
params.publish_dir_mode     = 'copy'
params.apptainer_cache_dir  = null

// ── Test fixtures ─────────────────────────────────────────────────────────────
// Use the existing NO_FILE sentinel for the legacy known_polya and polya_db
// (no legacy file, no PolyA_DB source).  The polyasite local-file path is
// reused from the standard test data so the process gets a real staged file.
def KNOWN_POLYA    = file("${launchDir}/assets/NO_FILE",          checkIfExists: true)
def POLYA_DB_FILE  = file("${launchDir}/assets/NO_FILE_POLYA_DB", checkIfExists: true)
def POLYASITE_FILE = file("${launchDir}/tests/data/known_polya.tsv", checkIfExists: true)

// ── Expected output schema ────────────────────────────────────────────────────
def REQUIRED_COLUMNS = [
    'gene_id', 'gene_name', 'chrom', 'start', 'end', 'strand', 'score', 'source'
]
def REQUIRED_MANIFEST_FIELDS = ['mode']

workflow {
    PREPARE_POLYA_REFERENCES(
        Channel.value(KNOWN_POLYA),
        Channel.value('skip'),                                        // polya_db val
        Channel.value(POLYA_DB_FILE),
        Channel.value("${launchDir}/tests/data/known_polya.tsv"),    // polyasite val
        Channel.value(POLYASITE_FILE)
    )

    // ── Assert prepared_known_polya.tsv header ────────────────────────────────
    PREPARE_POLYA_REFERENCES.out.known_polya
        .map { f ->
            assert f.exists()   : "prepared_known_polya.tsv does not exist"
            assert f.size() > 0 : "prepared_known_polya.tsv is empty"

            def lines  = f.readLines()
            def header = lines.first().split('\t').collect { it.trim() }
            REQUIRED_COLUMNS.each { col ->
                assert header.contains(col) :
                    "prepared_known_polya.tsv missing column '${col}'. Found: ${header}"
            }
            assert lines.size() >= 2 : "prepared_known_polya.tsv must have at least one data row"
            log.info "[PASS] prepared_known_polya.tsv — header OK: ${header}"
            f
        }
        .view { "prepared_known_polya.tsv: ${it.name} (${it.size()} bytes)" }

    // ── Assert manifest is non-empty and contains expected fields ─────────────
    PREPARE_POLYA_REFERENCES.out.manifest
        .map { f ->
            assert f.exists() : "manifest TSV does not exist"
            assert f.size() > 0 : "manifest TSV is empty"

            def fields = f.readLines()
                .findAll { it.trim() && !it.startsWith('field') }
                .collect { it.split('\t')[0].trim() }

            REQUIRED_MANIFEST_FIELDS.each { field ->
                assert fields.contains(field) :
                    "manifest missing field '${field}'. Found: ${fields}"
            }
            log.info "[PASS] manifest — all required fields present: ${fields}"
            f
        }
        .view { "manifest: ${it.name} (${it.size()} bytes)" }
}
