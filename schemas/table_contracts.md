# Table Contracts

## samplesheet.csv

`sample_id,library_id,protocol,chemistry,condition,replicate_id,fastq_r1,fastq_r2,bam,matrix_path,barcode_whitelist`

## cell_metadata.tsv

`sample_id,barcode_raw,barcode_corrected,cell_id,cluster_id,cell_type,condition,batch,label_source`

## group_map.tsv

`sample_id,barcode_corrected,group_level,group_id`

## site_catalog.tsv

`site_id,gene_id,chrom,start,end,strand,site_class,site_source,priming_flag`

## apa_usage.tsv

`group_level,group_id,gene_id,site_id,read_count,umi_count,usage_fraction,pdui_like`

## apa_stats.tsv

`contrast_id,gene_id,site_id,group_level,logFC,delta_usage,pvalue,fdr,test_class`
