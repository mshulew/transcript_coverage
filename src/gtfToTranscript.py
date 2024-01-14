#!/usr/bin/env python3
# coding: utf-8

"""
creates list of all transcripts per gene from gtf
"""

import sys

def getValue(annotation,query_metric):
    extractedValue = ''
    for element in annotation.split(';'):
        if len(element.split()) > 0:
            metric = element.split()[0]
            if metric == query_metric:
                extractedValue = element.split()[1].replace('"','')
                break
        
    return extractedValue

if __name__ == "__main__":

    gtf_filename = sys.argv[1]
    output_filename = sys.argv[2]
        
# create list of all transcripts per gene
    genelist = {}
    for line in open(gtf_filename, 'r'):
        line = line.strip('\n').split('\t')   
        if line[2] == 'exon':
            annotation = line[8]
            gene_name = getValue(annotation,'gene_id')
            transcript_id = getValue(annotation,'transcript_id')
                                     
            genelist.setdefault(gene_name,[])
            if transcript_id not in genelist[gene_name]:
                genelist[gene_name].append(transcript_id)
                
# write to file
    with open(output_filename, 'w') as outputfile:
        for genename,transcripts in genelist.items():
            for transcript in transcripts:
                outputfile.write('{}\t{}\n'.format(transcript,genename))
