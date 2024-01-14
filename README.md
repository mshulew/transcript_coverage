# Gene Coverage LSG0038

This pipeline was developed as a companion to gene_coverage_lsg0038 to support SEQuoia Express. It calculates the percent of bases of the top 1,000 expressed transcripts at different read coverages (0x - 20X). This pipeline was used in verification and validation and is being archived for reference. See improvements for bugs fixes for future use.

Author: Mark Shulewitz (mark_shulewitz@bio-rad.com) 

## Description of pipeline
### assignTranscripts
Run featureCounts using deduplicated BAM (if deduplication was performed) or STAR BAM (if no deduplication was performed) to assign reads to transcripts (note: reads can be assigned to more than 1 transcript)  

### transcriptList
Generates list of transcripts by gene from gtf annotation file  

### topTranscripts
Generates a list of most highly expressed transcript for each gene with read counts for that transcript 

### makeBed  
1. Identifies top 1,000 expressed transcripts (top transcript per gene)
2. Filters gtf annotation file for exons from these top 1,000 transcripts and generates a bed file

### filterBam
1.  Filters gtf annotation file for top 1,000 transcripts
2.  Rerun featureCounts using deduplicated BAM (if deduplication was performed) or STAR BAM (if no deduplication was performed) to assign reads to top 1,000 transcripts
3.  Filters for assigned alignments smaller than 1,000 bp (SEQuoia Express inserts should be ~300 bp; set to 1,000 bp to accomodate spliced introns while removing nonsense)

## transcriptCoverage  
Calculates base pair coverage using bed file from makeBed process and BAM file from filterBam step

## foldCoverage
1. Uses output from transcriptCoverage step to calculate fold coverage of each base (0X - 20X)
2. Calculates overall coverage across all bases of all top 1,000 transcripts


## Preflight

Install Docker and Nextflow
Build Docker Container

```bash
cd transcript_coverage_lsg0038
docker build -t transcripts_coverage:latest .
```

## Running the pipeline

```bash
nextflow run transcript_coverage_lsg0038/main.nf --inDir /path/to/SEQuoia-Express/output/directory --outDir /path/to/output/directory
```

## Description of output

{sample id}_foldcoverage.tsv - fold coverage across each of the top 1,000 genes 
  - percent of bases with 0X coverage (no reads), 1X (at least 1 read), 5X (at least 5 reads), 10X (at least 10 reads) and 20X (at least 20 reads)
  - also includes coverage across all genes ("all genes") - this is the stat used for the DITM for LSG0038

{sample id)_stats.tsv 
  - total bases in top 1,000 genes  
  - min, max, mean and median reads per base  
  - lower, middle and upper quartile  

{sample_id}_NS16_10M_coverage_by_base.cov - raw data, number of reads mapping to each base of top 1,000 genes 

### Intermediate files:  

assignTranscripts > 
Aligned.sortedByCoord.[deduplicated.]out.bam.featureCounts.bam - featureCounts output bam  
transcript_counts.tsv - read count file from featureCounts  

bed > 
toptranscripts - list of top 1,000 expressed transcripts (highest expressed transcript per gene)
top.bed - top transcript and read counts for each gene  

filterBam > filtered.bam - filtered featureCounts output Bam file from filterBam process

topTranscripts > top_transcripts.tsv List of most highly expressed transcript for each gene with read counts for that transcript from topTranscripts process
transcriptList > transcriptlist.tsv - list of transcripts and corresponding gene from transcriptList process  



## Running pipeline without Docker

It is possible to run this pipeline without Docker. Dockerfile can be used as a guide to install all of the required software in a Linux instance (python3.9, Subread, Bedtools, Sambamba). ADditional scripts (ie: bash and python scripts) are in the scr directory. main.nf can be used as a guide for the command lines that need to be run sequentially (command lines are in the script block for each process). 

## Additional resources

- The bash script topgene.sh can be used to compare the top 1,000 genes between processed data sets. 
- The bash script batchprocess.sh can be used to sequentially process multiple SEQuoia Express Toolkit outputs.

## Improvements
1. In the filterBam process, the featureCounts output BAM file is not correctly filtered for Assiged reads. In addition, reads with unmapped mates should be removed. samtools view *.featureCounts.bam | grep 'Assigned' | awk 'sqrt(\$9^2) < 1000' > body.sam should changed to samtools view -F 8 *.featureCounts.bam | grep 'XS:Z:Assigned' | awk 'sqrt(\$9^2) < 1000' > body.sam  
2. The transcriptCoverage process is not strand specific. Add strand = params.reverseStrand ? "-S" : "-s" as the first line in the script block
3. Update bedtools coverage -a $top_bed -b $top_bam -s -d > bedtools_coverage to bedtools coverage -a $top_bed -b $top_bam $strand -d > bedtools_coverage
