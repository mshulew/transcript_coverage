#!/usr/bin/env python3
# coding: utf-8

"""
creates bed file from gtf file; filters for genes from list
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

    input_filename = sys.argv[1]
    transcript_list_filename = sys.argv[2]
    output_filename = sys.argv[3]
    transcript_list = []
    bed_data = []
    
# import gene_list
    for line in open(transcript_list_filename, 'r'):
        transcript_list.append(line.strip('\n'))
    
# iterate through gtf file
    for line in open(input_filename, 'r'):
        line = line.strip('\n').split('\t')
        if line[2] == 'exon':
            chromosome = line[0]
            startpos = int(line[3]) - 1
            stoppos = line[4]
            annotation = line[8]
            strand = line[6]
            source = line[1]
            
            gene_name = getValue(annotation,'gene_id')
            transcript_id = getValue(annotation,'transcript_id')
            exon_number = getValue(annotation, 'exon_number')
            db_xref = getValue(annotation,'db_xref')
            gbkey = getValue(annotation,'gbkey')
            gene = getValue(annotation,'gene')
                                     
            new_annotation = 'gene_id "{}"; transcript_id "{}"; exon_number "{}"; db_xref "{}"; gbkey "{}"; gene "{}";'.format(gene_name,\
                                                                                                                               transcript_id,\
                                                                                                                               exon_number,\
                                                                                                                               db_xref,\
                                                                                                                               gbkey,\
                                                                                                                               gene)
            
            for element in annotation.split(';'):
                if 'gene_id' in element:
                    gene_name = element.split('gene_id ')[1].replace('"','')
                    break
                     
            if transcript_id in transcript_list:
                bed_data.append([chromosome,str(startpos),stoppos,gene_name,transcript_id,strand,source,'exon',line[5],new_annotation])
                
# write to file
    with open(output_filename, 'w') as outputfile:
        for entry in bed_data:
            outputfile.write('{}\n'.format('\t'.join(entry)))
