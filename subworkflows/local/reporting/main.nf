include { RENDER_APA_REPORT } from '../../../modules/local/report/render_apa_report/main'

workflow REPORTING_SC {
    take:
    site_catalog
    apa_usage
    apa_stats
    track_bundle
    qc_bundle

    main:
    RENDER_APA_REPORT(site_catalog, apa_usage, apa_stats, track_bundle.collect(), qc_bundle.collect())

    emit:
    report  = RENDER_APA_REPORT.out.report
    plots   = RENDER_APA_REPORT.out.plots
    summary = RENDER_APA_REPORT.out.summary
}
