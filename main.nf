#!/usr/bin/env nextflow

"""
Coverage of top 1,000 transcripts
2024 Mark Shulewitz (markshulewitz@gmail.com)

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

//include modules containing processes
include { assignTranscripts }   from './modules/assignTranscripts'
include { transcriptList }      from './modules/transcriptList'
include { topTranscripts }      from './modules/topTranscripts'
include { makeBed }             from './modules/makeBed'
include { filterBam }           from './modules/filterBam'
include { transcriptCoverage }  from './modules/transcriptCoverage'
include { foldCoverage }        from './modules/foldCoverage'

longRNAgtfFile = file(params.genomes[params.genome][params.spikeType].longRNAgtfFile)

inputDir = file("${params.inDir}", checkIfExists: true)

// check for deduplication and get input bam file
dedup = true

bam_file = Channel
  .fromPath("${params.inDir}/dedup/*rumi_dedup.sort.bam")
  .ifEmpty { dedup = false }

if (!dedup){
  bam_file = Channel
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

workflow {
        transcriptList(longRNAgtfFile)
        assignTranscripts(namingprefix,bam_file,longRNAgtfFile)
        topTranscripts(assignTranscripts.out.transcript_read_counts,transcriptList.out)
        makeBed(topTranscripts.out,longRNAgtfFile)
        filterBam(makeBed.out.top_transcripts,longRNAgtfFile,bam_file)
        transcriptCoverage(makeBed.out.top_bed,filterBam.out)
        foldCoverage(transcriptCoverage.out)
}



def Header() {
    return """
    Transcript coverage
    version 2.0
    January 21, 2024
    """.stripIndent()
    }