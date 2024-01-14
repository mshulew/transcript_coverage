#!/usr/bin/env python3
# coding: utf-8

"""
creates list of top transcript for each gene with read counts for that transcript
"""

import sys

if __name__ == "__main__":

    transcripts_filename = sys.argv[1]
    counts_filename = sys.argv[2]
    output_filename = sys.argv[3]
        
# create dictionary of all transcripts with gene name
    transcriptgene = {}
    for line in open(transcripts_filename, 'r'):
        line = line.strip('\n').split('\t') 
        transcriptgene.setdefault(line[0],line[1])
        
# create dictionary of all genes with transcripts and read count for transcript
    readcounts = {}
    for line in open(counts_filename, 'r'):
        line = line.strip('\n').split('\t')
        transcriptname = line[0]
        genename = transcriptgene[transcriptname]
        transcriptreadcount = int(line[1])
        readcounts.setdefault(genename,[])
        readcounts[genename].append([transcriptname,transcriptreadcount])
        
# identify transcript per gene with highest read count
    toptranscript = []
    for gene,transcripts in readcounts.items():
        max_transcript = ['TBD',-1]
        for transcript in transcripts:
            transcriptname = transcript[0]
            transcriptreadcounts = transcript[1]
            if transcriptreadcounts > max_transcript[1]:
                max_transcript = [transcriptname, transcriptreadcounts]
        toptranscript.append([gene,max_transcript[0],str(max_transcript[1])])
        
                
# write to file
    with open(output_filename, 'w') as outputfile:
        for entry in toptranscript:
                outputfile.write('{}\n'.format('\t'.join(entry)))
