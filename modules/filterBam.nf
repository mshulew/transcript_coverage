process filterBam {
  tag "filterBam on ${name}"
  label 'mid_cpu'
  publishDir "${params.outDir}/filterBam", mode: 'copy'

  input:
  tuple val(name),path(toptranscripts)
  path longRNAgtfFile
  path bam

  output:
  tuple val(name),path('filtered.bam')

  script:
  strand = params.reverseStrand ? "-s 2" : "-s 1"
  """
  grep -w -F -f $toptranscripts $longRNAgtfFile > toptranscripts.gtf
  featureCounts -p -T $task.cpus --primary -M -t exon -g transcript_id $strand -Q 1 --fracOverlap 0.8 \
   -a toptranscripts.gtf \
   -o ./transcript_counts \
   -R BAM ${bam}
  samtools view -H *featureCounts.bam > header.sam
  samtools view -F 8 *.featureCounts.bam | grep 'XS:Z:Assigned' | awk 'sqrt(\$9^2) < 1000' > body.sam
  cat header.sam body.sam | samtools view -Sb | samtools sort > filtered.bam
  sambamba index -t ${task.cpus} filtered.bam
  """
}
