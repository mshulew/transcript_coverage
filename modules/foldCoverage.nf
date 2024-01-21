process foldCoverage {
  tag "foldCoverage on ${name}"
  publishDir "${params.outDir}", mode: 'copy'

  input:
  tuple val(name),path(fullcoverage)

  output:
  path '*.tsv'

  script:
  """
  python3.9 /opt/biorad/src/foldCoverage.py $fullcoverage ${name}_foldcoverage.tsv
  mv stats ${name}_stats.tsv
  """
}
