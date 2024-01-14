#!/usr/bin/env nextflow

"""
Coverage of top 1,000 transcripts
2021 Mark Shulewitz Bio-Rad Laboratories, Inc.

"""

def helpMessage() {
  log.info Header()
  log.info """
  Usage:

  The typical command for running the pipeline is as follows:

  Required:
    --inDir                 Path to directory containing SEQuoia Express Toolkit output for sample
    --outDir                Path to output directory
    -profile docker         Requires docker container
 
    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

longRNAgtfFile = file(params.genomes[params.genome][params.spikeType].longRNAgtfFile)

inputDir = file("${params.inDir}", checkIfExists: true)

// check for deduplication and get input bam file
dedup = true
    
Channel
  .fromPath("${params.inDir}/dedup/*rumi_dedup.sort.bam")
  .ifEmpty { dedup = false }
  .into { bam1_ch; bam2_ch }
    
if (!dedup){
  Channel
    .fromPath("${params.inDir}/star/*sortedByCoord.out.bam")
    .ifEmpty { exit 1, "STAR BAM file not found" }
    .into { bam1_ch; bam2_ch }  
}



if (params.inDir.toString().substring(params.inDir.toString().length() - 1) == "/"){
  namingprefix = params.inDir.toString().substring(0,params.inDir.toString().length() - 1).substring(params.inDir.toString().substring(0,params.inDir.toString().length() - 1).lastIndexOf("/") + 1)
  
} else { namingprefix = params.inDir.toString().substring(params.inDir.toString().lastIndexOf("/") + 1)
}

def summary = [:]
summary['Run Name'] = workflow.runName
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Config Profile'] = workflow.profile
summary['Input directory'] = params.inDir
summary['Name Prefix'] = namingprefix
summary['Output directory'] = params.outDir
summary['Deduplication'] = dedup
summary['Reverse strand'] = params.reverseStrand
log.info Header()
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"
 
process assignTranscripts {
  tag "assignTranscripts on ${namingprefix}"
  label 'mid_cpu'
  publishDir "${params.outDir}/assignTranscripts", mode: 'copy'

  input:
  file bam from bam1_ch
  file longRNAgtfFile

  output:
  file 'transcript_counts.tsv' into transcript_counts_ch
  file '*featureCounts.bam'

  script:
  strand = params.reverseStrand ? "-s 2" : "-s 1"
  """
  featureCounts -p -T $task.cpus --primary -M -O -t exon -g transcript_id $strand -Q 1 --fracOverlap 0.8 \
   -a $longRNAgtfFile \
   -o ./transcript_counts \
   -R BAM ${bam} 
  cut -f 1,7 transcript_counts | tail -n +3 > transcript_counts.tsv
  """
}
  
process transcriptList {
  tag "transcriptList on ${namingprefix}"
  publishDir "${params.outDir}/transcriptList", mode: 'copy'
    
  input:
  file longRNAgtfFile
    
  output:
  file 'transcriptlist.tsv' into transcriptlist_ch
    
  script:
  """
  python3.6 /opt/biorad/src/gtfToTranscript.py ${longRNAgtfFile} 'transcriptlist.tsv'
  """
}
  
process topTranscripts {
  tag "topTranscripts on ${namingprefix}"
  publishDir "${params.outDir}/topTranscripts", mode: 'copy'
  
  input:
  file transcriptlist from transcriptlist_ch
  file transcriptcounts from transcript_counts_ch
  
  output:
  file 'top_transcripts.tsv' into makebed_ch
  
  script:
  """
  python3.6 /opt/biorad/src/TopTranscripts.py ${transcriptlist} ${transcriptcounts} 'top_transcripts.tsv'
  """
}

process makeBed {
  tag "makeBed on ${namingprefix}"
  publishDir "${params.outDir}/bed", mode: 'copy'
  
  input:
  file longRNAgtfFile
  file transcriptcounts from makebed_ch
  
  output:
  file 'top.bed' into bed_ch
  file 'toptranscripts' into toptranscripts_ch
  
  script:
  """
  bash /opt/biorad/src/top1000.sh $transcriptcounts toptranscripts
  python3.6 /opt/biorad/src/gtftobed_filter.py ${longRNAgtfFile} toptranscripts top_bed
  awk 'BEGIN{FS=" ";OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,\$6}' top_bed > top.bed
  """ 
}

process filterBam {
  tag "filterBam on ${namingprefix}"
  label 'mid_cpu'
  publishDir "${params.outDir}/filterBam", mode: 'copy'

  input:
  file toptranscripts from toptranscripts_ch
  file longRNAgtfFile
  file bam from bam2_ch

  output:
  file 'filtered.bam' into coverage_bam_ch

  script:
  strand = params.reverseStrand ? "-s 2" : "-s 1"
  """
  grep -w -F -f $toptranscripts $longRNAgtfFile > toptranscripts.gtf
  featureCounts -p -T $task.cpus --primary -M -t exon -g transcript_id $strand -Q 1 --fracOverlap 0.8 \
   -a toptranscripts.gtf \
   -o ./transcript_counts \
   -R BAM ${bam} 
  samtools view -H *featureCounts.bam > header.sam
  samtools view *featureCounts.bam | grep 'Assigned' | awk 'sqrt(\$9^2) < 1000' > body.sam
  cat header.sam body.sam | samtools view -Sb | samtools sort > filtered.bam
  sambamba index -t ${task.cpus} filtered.bam
  """
}


process transcriptCoverage {
  tag "transcriptCoverage on ${namingprefix}"
  publishDir "${params.outDir}", mode: 'copy'
  label 'mid_memory'
  
  input:
  file top_bed from bed_ch
  file top_bam from coverage_bam_ch
    
  output:
  file '*.tsv' into fullcoverage_ch
    
  script:
  """
  bedtools coverage -a $top_bed -b $top_bam -s -d > bedtools_coverage
  awk 'BEGIN {FS="\t";OFS="\t"} { print \$1, \$2 + \$7 - 1, \$4, \$5, \$8 }' bedtools_coverage | sort -k1,1 -k2n,2 > ${namingprefix}_coverage_by_base.tsv
  """
}

process foldCoverage {
  tag "foldCoverage on ${namingprefix}"
  publishDir "${params.outDir}", mode: 'copy'
    
  input:
  file fullcoverage from fullcoverage_ch
    
  output:
  file '*.tsv'
    
  script:
  """
 python3.9 /opt/biorad/src/foldCoverage.py $fullcoverage ${namingprefix}_foldcoverage.tsv 
 mv stats ${namingprefix}_stats.tsv
  """
}

 
def Header() {
    return """
    Transcript coverage
    version 1.0
    August 23, 2021
    """.stripIndent()
    }
