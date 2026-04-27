process RENDER_APA_REPORT {
    tag "render_apa_report"
    label 'process_single'
    label 'process_python'

    conda "${projectDir}/envs/python.yml"
    container params.apptainer_cache_dir ? "${params.apptainer_cache_dir}/scpolaseq-python.sif" : null
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    input:
    path site_catalog
    path apa_usage
    path apa_stats
    path track_paths, stageAs: 'tracks/?/*'
    path qc_paths,    stageAs: 'qc/?/*'

    output:
    path "scpolaseq_report.html", emit: report
    path "report_plots", emit: plots
    path "report_summary.tsv", emit: summary

    script:
    def trackArgs = track_paths.collect { it.toString() }.join(' ')
    def qcArgs = qc_paths.collect { it.toString() }.join(' ')
    """
    python3 ${projectDir}/bin/render_apa_report.py \\
        --site-catalog ${site_catalog} \\
        --apa-usage ${apa_usage} \\
        --apa-stats ${apa_stats} \\
        --track-dir tracks \\
        --qc-dir qc \\
        --out-html scpolaseq_report.html \\
        --out-plots report_plots \\
        --out-summary report_summary.tsv
    """

    stub:
    """
    mkdir -p report_plots
    echo "<html><body><h1>scPolASeq report</h1></body></html>" > scpolaseq_report.html
    printf "artifact_type\\tcount\\ntrack_paths\\t0\\nqc_paths\\t0\\n" > report_summary.tsv
    """
}
