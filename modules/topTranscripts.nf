process topTranscripts {
  tag "topTranscripts on ${name}"
  publishDir "${params.outDir}/topTranscripts", mode: 'copy'

  input:
  tuple val(name),path(transcriptcounts)
  path transcriptlist

  output:
  tuple val(name),path('top_transcripts.tsv')

  script:
  """
  python3.6 /opt/biorad/src/TopTranscripts.py ${transcriptlist} ${transcriptcounts} 'top_transcripts.tsv'
  """
}
