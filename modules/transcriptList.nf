process transcriptList {
  tag "transcriptList"
  publishDir "${params.outDir}/transcriptList", mode: 'copy'

  input:
  path longRNAgtfFile

  output:
  path 'transcriptlist.tsv'

  script:
  """
  python3.6 /opt/biorad/src/gtfToTranscript.py ${longRNAgtfFile} 'transcriptlist.tsv'
  """
}
