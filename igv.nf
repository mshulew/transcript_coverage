#!/usr/bin/env nextflow

"""
Creates percent gc for igv (integrated genome viewer) 
"""

def helpMessage() {
  log.info Header()
  log.info """
  Usage:

  The typical command for running the pipeline is as follows:

  Options:
    --fasta                 Path to genome reference fasta file; default is ./ref_data/genome_annotations/hg38/fasta_ncbi/hg38_no_alt_analysis_set.fna' 
    --outDir                Path to output directory; default is ./results
    --width		    Window size in bases; default is 50 bp
    --chromosome	    limit coverage to a single chromosome by number; default is 'all'

    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

//include modules containing processes
include { igv }      from './modules/igv.nf'

fasta_file = file("${params.fasta}", checkIfExists: true)

if(fasta_file.isEmpty()){
 print "genome fasta file not found"
 exit 1
}


def summary = [:]
summary['Run Name'] = workflow.runName
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Config Profile'] = workflow.profile
summary['Window'] = params.width
summary['Chromosome'] = params.chromosome
summary['Output Dir']=params.outDir
log.info Header()
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"

workflow {
	igv(fasta_file)
}



def Header() {
    return """
    IGV reference file generator
    version 1.0
    January 22, 2024
    """.stripIndent()
    }
