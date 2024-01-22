# Transcript Coverage

This toolkit calculates the percent of bases of the top 1,000 expressed transcripts at different read coverages (0x - 20X). This toolkit is a companion for Bio-Rad's SEQuoia Express Toolkit, using the SEQuoia Express Toolkit output as its input. With slight modification to main.nf, it can be used for any BAM file with paired-end read alignments.

Author: Mark Shulewitz (markshulewitz@gmail.com)

## Description of main workflow (main.nf)
### assignGened
Runs featureCounts using deduplicated BAM (if deduplication was performed) or STAR BAM (if no deduplication was performed) to assign reads to genes. Output is filtered for alignments assigned to reads and is not piped to any downstream process in the toolkit - the BAM and indexed BAM files are for viewing on IGV or similar genome viewer.

### transcriptList
Generates list of transcripts by gene from gtf annotation file

### assignTranscripts
Runs featureCounts using deduplicated BAM (if deduplication was performed) or STAR BAM (if no deduplication was performed) to assign reads to transcripts (note: reads can be assigned to more than 1 transcript)   

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
nextflow run transcript_coveragemain.nf --inDir /path/to/SEQuoia-Express/output/directory --outDir /path/to/output/directory
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

- igv.nf is used to generate a file of percent GC content for IGV or a similar genome viewer. Run using command line nextflow transcript_coverage/igv.nf -profile docker --help for list of options.
- The bash script topgene.sh can be used to compare the top 1,000 genes between processed data sets. 
- The bash script batchprocess.sh can be used to sequentially process multiple SEQuoia Express Toolkit outputs.
