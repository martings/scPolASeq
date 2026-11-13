nextflow.enable.dsl = 2

include { PREPARE_REFERENCE_BUNDLE } from '../../../modules/local/reference/prepare_reference_bundle/main'

params.fixture = null

workflow {
    if (!params.fixture) {
        error "Missing required --fixture directory"
    }
    PREPARE_REFERENCE_BUNDLE(
        Channel.value(file("${params.fixture}/genomes/GRCh38/genome.fa", checkIfExists: true)),
        Channel.value(file("${params.fixture}/genomes/GRCh38/genes.gtf", checkIfExists: true)),
        Channel.value(file("${params.fixture}/databases/polya_atlas/polyasite2_GRCh38.tsv", checkIfExists: true)),
        Channel.value(file("${projectDir}/../../../assets/NO_FILE_PRIMING_BLACKLIST", checkIfExists: true)),
        Channel.value(file("${params.fixture}/indexes/GRCh38/star_index", type: 'dir', checkIfExists: true)),
        Channel.value(file("${params.fixture}/genomes/GRCh38/chrom.sizes", checkIfExists: true))
    )
}
