process transcriptCoverage {
  tag "transcriptCoverage on ${name}"
  publishDir "${params.outDir}", mode: 'copy'
  label 'mid_memory'

  input:
  path(top_bed)
  tuple val(name),path(top_bam)

  output:
  tuple val(name),path('*.tsv')

  script:
  strand = params.reverseStrand ? "-S" : "-s"
  """
  bedtools coverage -a $top_bed -b $top_bam $strand -d > bedtools_coverage
  awk 'BEGIN {FS="\t";OFS="\t"} { print \$1, \$2 + \$7 - 1, \$4, \$5, \$8 }' bedtools_coverage | sort -k1,1 -k2n,2 > ${name}_coverage_by_base.tsv
  """
}
