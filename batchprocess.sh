#!/bin/bash
# batch process gc content

inDir=$1
outDir=$2

for dir in $inDir/*; do
  prefix=$(echo $dir | awk '{print $1}' FS=\_ RS=\/ | tail -1)
  echo "processing ${dir} (${prefix}) ..."
  nextflow run transcript_coverage/main.nf -profile docker --inDir $dir --outDir ${outDir}/${prefix} 
  sudo rm -r work
done



