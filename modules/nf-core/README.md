# nf-core vendoring target

This scaffold keeps custom local modules as the executable default and tracks preferred nf-core equivalents in [`../../modules.json`](../../modules.json).

Recommended vendoring targets for later hardening:

- `star/genomegenerate`
- `star/align`
- `star/starsolo`
- `samtools/*`
- `bedtools/*`
- `ucsc/bedclip`
- `ucsc/bedgraphtobigwig`
- `multiqc`
- `gffread`
