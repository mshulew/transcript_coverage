process assignGenes {
  tag "assignGenes on ${name}"
  label 'mid_cpu'
  publishDir "${params.outDir}/assignGenes", mode: 'copy'

  input:
  val name
  path bam
  path longRNAgtfFile

  output:
  path '*featureCounts.bam*'
  path 'gene_counts.tsv'

  script:
  strand = params.reverseStrand ? "-s 2" : "-s 1"
  """
  featureCounts -p -T $task.cpus --primary -M -t exon -g gene_id $strand -Q 1 --fracOverlap 0.8 \
   -a $longRNAgtfFile \
   -o ./gene_counts \
   -R BAM ${bam}
  cut -f 1,7 gene_counts | tail -n +3 > gene_counts.tsv
  """
}
