process assignTranscripts {
  tag "assignTranscripts on ${name}"
  label 'mid_cpu'
  publishDir "${params.outDir}/assignTranscripts", mode: 'copy'

  input:
  val name
  path bam
  path longRNAgtfFile

  output:
  tuple val(name), path('transcript_counts.tsv'), emit:transcript_read_counts
  path '*featureCounts.bam'

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