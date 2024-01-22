process igv {
  tag "create igv ref files"
  label 'mid_cpu'
  publishDir "${params.outDir}", mode: 'copy'

  input:
  path fasta

  output:
  path 'genomeref.fa*'
  path '*.igv'

  script:
  window = params.width
  if(params.chromosome == 'all'){
  """
  cp $fasta genomeref.fa
  samtools faidx genomeref.fa
  cut -f 1,2 genomeref.fa.fai > genomeref.sizes
  bedtools makewindows -g genomeref.sizes -w $window > genome_${window}_bps.bed
  bedtools nuc -fi genomeref.fa -bed genome_${window}_bps.bed > genome_nuc_${window}.txt
  gawk -v w=$window 'BEGIN{FS="\t"; OFS="\t"} { if (FNR>1) {print \$1,\$2,\$3, "GCpc_"w"bps",\$5} }' genome_nuc_${window}.txt > genome_nuc_${window}.igv
  """
  } else {
  """
  cp $fasta genomeref.fa
  samtools faidx genomeref.fa
  cut -f 1,2 genomeref.fa.fai > genomeref.sizes
  bedtools makewindows -g genomeref.sizes -w $window > genome_${window}bps.bed
  mv genome_{window}bps.bed "temp.bed
  awk -v chromosome=$chromosome '\$1 == chromosome' temp.bed > genome_${window}bps.bed
  rm temp.bed
  bedtools nuc -fi genomeref.fa -bed genome_${window}_bps.bed > genome_nuc_${window}.txt
  gawk -v w=$window 'BEGIN{FS="\t"; OFS="\t"} { if (FNR>1) {print \$1,\$2,\$3, "GCpc_"w"bps",\$5} }' genome_nuc_${window}.txt > genome_nuc_${window}.igv
  mv genome_nuc_${window}.igv > chromosome_${chromosome}_nuc_${window}.igv
  """
  }
}
