"""
calculates fold coverage from bedtools coverage -d output
"""

import sys
import statistics

if __name__ == "__main__":

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
# iterate through bedtools coverage output and count coverage at each base for all genes (all_cov) and each individual gene (each_cov)
# entries are as follows: total bases, bases with 0 reads (0X), bases with >0 reads (1X), bases with >4 reads (5X), bases with >9 reads (10X), bases with >19 reads (>20X)

    all_cov = [0,0,0,0,0,0]
    each_cov = {}
    fold_cov = []
    transcript_key = {}
   
    for line in open(input_filename, 'r'):

        
        line = line.strip('\n').split('\t')
        gene = line[2]
        transcript = line[3]
        readperbp = int(line[4])
        
        transcript_key.setdefault(transcript,gene)
        
        fold_cov.append(readperbp)
        
# create entry for gene if does not exist
        each_cov.setdefault(transcript,[0,0,0,0,0,0])
        
# total bases +1
        all_cov[0] += 1
        each_cov[transcript][0] += 1
        
# 0X
        if readperbp == 0:
            all_cov[1] += 1
            each_cov[transcript][1] += 1
           
# 1X
        else:
            all_cov[2] += 1
            each_cov[transcript][2] += 1
# 5X            
            if readperbp > 4:
                all_cov[3] += 1
                each_cov[transcript][3] += 1
                
# 10X 
                if readperbp > 9:
                    all_cov[4] += 1
                    each_cov[transcript][4] += 1
                
# 20X
                    if readperbp > 19:
                        all_cov[5] += 1
                        each_cov[transcript][5] += 1
                
# convert each_cov to list
    each_coverage = []
    for transcript,stats in each_cov.items():
        gene = transcript_key[transcript]
        each_coverage.append([gene,transcript,stats[0],stats[1],stats[2],stats[3],stats[4],stats[5]])
        
# write to file
    with open(output_filename, 'w') as outputfile:
        outputfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('gene','transcript','0X','1X','5X','10X','20X'))
        zeroX = round((all_cov[1]/all_cov[0])*100)
        oneX = round((all_cov[2]/all_cov[0])*100)
        fiveX = round((all_cov[3]/all_cov[0])*100)
        tenX = round((all_cov[4]/all_cov[0])*100)
        twentyX = round((all_cov[5]/all_cov[0])*100)
        outputfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('all genes','all transcripts',zeroX,oneX,fiveX,tenX,twentyX))
        for entry in sorted(each_coverage):
            gene = entry[0]
            transcript = entry[1]
            totalreads = entry[2]
            zeroX = round((entry[3]/totalreads)*100)
            oneX = round((entry[4]/totalreads)*100)
            fiveX = round((entry[5]/totalreads)*100)
            tenX = round((entry[6]/totalreads)*100)
            twentyX = round((entry[7]/totalreads)*100)
            
            outputfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene,transcript,zeroX,oneX,fiveX,tenX,twentyX))
            
# stats

    with open('stats', 'w') as outputfile:
        outputfile.write('{}\t{}\n'.format('total bases',len(fold_cov)))
        outputfile.write('{}\t{}\n'.format('min reads per base',min(fold_cov)))
        outputfile.write('{}\t{}\n'.format('max reads per base',max(fold_cov)))
        outputfile.write('{}\t{}\n'.format('mean reads per base',round(statistics.mean(fold_cov),1)))
        outputfile.write('{}\t{}\n'.format('median reads per base',statistics.median(fold_cov)))
        outputfile.write('{}\t{}\n'.format('lower quartile',statistics.quantiles(fold_cov)[0]))
        outputfile.write('{}\t{}\n'.format('middle quartile',statistics.quantiles(fold_cov)[1]))
        outputfile.write('{}\t{}\n'.format('upper quartile',statistics.quantiles(fold_cov)[2]))

