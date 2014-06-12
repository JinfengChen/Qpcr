#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import math
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Qpcr.py --input date.txt

input: 
Well	Fluor	Target	Content	Sample	Cq
B02	SYBR	ubq	Unkn-01	nb	20.19
B03	SYBR	ubq	Unkn-01	nb	20.24
B04	SYBR	ubq	Unkn-01	nb	20.30
B05	SYBR	ubq	Unkn-01	nb	20.29
B06	SYBR	ubq	Unkn-02	heg4	18.75
B08	SYBR	ubq	Unkn-02	heg4	18.67
B09	SYBR	ubq	Unkn-02	heg4	18.85
C02	SYBR	ef1	Unkn-03	nb	18.14
C03	SYBR	ef1	Unkn-03	nb	18.23
C05	SYBR	ef1	Unkn-03	nb	18.16
C06	SYBR	ef1	Unkn-04	heg4	16.62
C07	SYBR	ef1	Unkn-04	heg4	16.60
C09	SYBR	ef1	Unkn-04	heg4	16.77


    '''
    print message



def sqr(x):
    return x*x

def readinput(infile, outfile, plot):
    ctl_gene = 'ubq'
    ctl_sample = 'nb'
    '''Raw data input'''
    data = defaultdict(lambda : defaultdict(lambda: list()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('Well'): 
                unit = re.split(r'\t',line)
                data[unit[2]][unit[4]].append(float(unit[5]))
    '''Mean and std'''
    data_sum = defaultdict(lambda : defaultdict(lambda: list())) 
    for gene in data.keys():
        for sample in data[gene].keys():
            #print gene, sample, np.mean(data[gene][sample]), np.std(data[gene][sample]), data[gene][sample]
            data_sum[gene][sample] = [np.mean(data[gene][sample]), np.std(data[gene][sample])]
            #print gene, sample, data_sum[gene][sample]
           
    
    #DeltaCt
    genes = ['Sample']
    data_final = defaultdict(lambda : defaultdict(lambda: list()))
    data_expr = defaultdict(lambda : defaultdict(lambda: list()))
    for gene in sorted(data_sum.keys()):
        if gene not in ctl_gene:
            genes.append(gene)
            genes.append('std')
            for sample in data_sum[gene].keys():
                DeltaCt = float(data_sum[gene][sample][0]) - float(data_sum[ctl_gene][sample][0])
                S_std   = math.sqrt(sqr(float(data_sum[gene][sample][1])) + sqr(float(data_sum[ctl_gene][sample][1])))
                data_expr[gene][sample] = [DeltaCt, S_std]
                expr_r  = 2**(-DeltaCt)
                Range_botton = 2**(-DeltaCt-data_expr[gene][sample][1])
                Range_top    = 2**(-DeltaCt+data_expr[gene][sample][1])
                STDEV        = ((expr_r-Range_botton)+(Range_top-expr_r))/2 
                data_final[sample][gene] = [expr_r, STDEV]

 
    ofile = open(outfile, 'w')
    #header = '\t'.join(headers)     
    header = '\t'.join(genes)
    print >> ofile, header
    for sample in sorted(data_final.keys()):
        expr = [sample]
        for gene in sorted(data_final[sample].keys()):
            expr.append(str(data_final[sample][gene][0]))
            expr.append(str(data_final[sample][gene][1]))
        print >> ofile, '\t'.join(expr)
    ofile.close()
    if plot == 'bar': 
        R_cmd_bar(outfile)
    #elif plot == 'point':
    #    R_cmd_point(outfile)

def R_cmd_bar(exprsumfile):
    R = '''
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

pdf("Qpcr.pdf",width=2, height=3) # width=2 if two sample, 3 if three sample.
data =read.table("''' + exprsumfile  + '''", header = T)
par(mar=c(3,4,2,2))
for (i in seq(2,length(data),by=2)){
    index = 1
    if ( nchar(floor(1/max(data[,i]))) == 4 ){
        index = 1000
    } else if ( nchar(floor(1/max(data[,i]))) == 3 ){
        index = 100
    } else if ( nchar(floor(1/max(data[,i]))) == 2 ){
        index = 10
    } else if ( nchar(floor(1/max(data[,i]))) == 1 ){
        if (floor(1/max(data[,i])) > 0){
           index = 1
        }
    }

    barx=barplot(data[,i]*index, ylim=c(0,max(data[,i]*index*1.3)), space = 0.6, border=F,axis.lty=1, ylab=paste(names(data)[i],"Ubq",index,sep='/'))
    error.bar(barx, data[,i]*index, data[,i+1]*index)
    axis(1,c(0.6,max(barx)+0.6),line=0,labels=c("",""))
    text(barx,rep(-0.07,length(data[,i])),offset=2,labels=data$Sample,srt=0,xpd=TRUE)
}
dev.off()

''' 
    infile = open ('Qpcr.R', 'w')
    print >> infile, R
    infile.close() 
    cmd = 'cat Qpcr.R | R --slave'
    #os.system(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-p', '--plot')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if args.output is None:
        args.output = 'Qpcr.sum'
    if args.plot is None:
        args.plot = 'bar'
    readinput(args.input, args.output, args.plot)

if __name__ == '__main__':
    main()

