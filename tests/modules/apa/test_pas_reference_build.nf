/*
 * tests/modules/apa/test_pas_reference_build.nf
 *
 * Standalone module-level test for PAS_REFERENCE_BUILD.
 *
 * Purpose:
 *   Validates the output contract of pas_reference_build/main.nf both while
 *   the module runs the scaffold shell commands AND after bin/build_pas_reference.py
 *   is wired in as the real implementation.
 *
 * Usage (stub — fast, no real tools needed):
 *   nextflow run tests/modules/apa/test_pas_reference_build.nf -stub-run
 *
 * Usage (real — exercises the actual script):
 *   nextflow run tests/modules/apa/test_pas_reference_build.nf
 *
 * Exit 0 = all assertions passed.
 */

nextflow.enable.dsl = 2

include { PAS_REFERENCE_BUILD } from '../../../modules/local/apa/pas_reference_build/main'

// ── Minimal params required by publishDir etc ─────────────────────────────────
params.outdir           = 'tests/results/module_pas_reference_build'
params.publish_dir_mode = 'copy'

// ── Test fixtures (reuse workflow test data) ─────────────────────────────────
// launchDir = directory from which `nextflow run` was called (repo root).
// projectDir = directory containing this script (tests/modules/apa/).
def SITE_CATALOG   = file("${launchDir}/tests/data/site_catalog.tsv",     checkIfExists: true)
def TERMINAL_EXONS = file("${launchDir}/tests/data/terminal_exons.tsv",   checkIfExists: true)
def KNOWN_POLYA    = file("${launchDir}/tests/data/known_polya.tsv",      checkIfExists: true)

// ── Required header columns in pas_reference.tsv ────────────────────────────
def REQUIRED_COLUMNS = [
    'pas_reference_id', 'site_id', 'gene_id',
    'chrom', 'start', 'end', 'strand', 'reference_source'
]

// ── Required fields in pas_reference_build.manifest.tsv ─────────────────────
def REQUIRED_MANIFEST_FIELDS = ['site_catalog', 'terminal_exons', 'known_polya']

workflow {
    PAS_REFERENCE_BUILD(
        Channel.value(SITE_CATALOG),
        Channel.value(TERMINAL_EXONS),
        Channel.value(KNOWN_POLYA)
    )

    // ── Assert pas_reference.tsv header ─────────────────────────────────────
    PAS_REFERENCE_BUILD.out.pas_reference
        .map { f ->
            assert f.exists()         : "pas_reference.tsv does not exist"
            assert f.size() > 0       : "pas_reference.tsv is empty"

            def lines = f.readLines()
            def header = lines.first().split('\t').collect { it.trim() }
            REQUIRED_COLUMNS.each { col ->
                assert header.contains(col) :
                    "pas_reference.tsv missing column '${col}'. Found: ${header}"
            }
            assert lines.size() >= 2 : "pas_reference.tsv must include at least one data row"
            log.info "[PASS] pas_reference.tsv — header OK: ${header}"
            f
        }
        .view { "pas_reference.tsv: ${it.name} (${it.size()} bytes)" }

    // ── Assert manifest lists expected input fields ──────────────────────────
    PAS_REFERENCE_BUILD.out.manifest
        .map { f ->
            assert f.exists() : "manifest TSV does not exist"

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

    // ── Assert log file exists and is non-empty ──────────────────────────────
    PAS_REFERENCE_BUILD.out.log
        .map { f ->
            assert f.exists() : "log file does not exist"
            assert f.size() > 0 : "log file is empty"
            log.info "[PASS] log file present (${f.size()} bytes)"
            f
        }
        .view { "log: ${it.name} (${it.size()} bytes)" }
}
