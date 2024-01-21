process makeBed {
  tag "makeBed on ${name}"
  publishDir "${params.outDir}/bed", mode: 'copy'

  input:
  tuple val(name),path(transcriptcounts)
  path longRNAgtfFile

  output:
  tuple val(name),path('toptranscripts'), emit:top_transcripts
  path 'top.bed', emit:top_bed


  script:
  """
  bash /opt/biorad/src/top1000.sh $transcriptcounts toptranscripts
  python3.6 /opt/biorad/src/gtftobed_filter.py ${longRNAgtfFile} toptranscripts top_bed
  awk 'BEGIN{FS=" ";OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,\$6}' top_bed > top.bed
  """
}
