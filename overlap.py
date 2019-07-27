# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 22:37:36 2019

@author: Maosong Pei
"""

from bx.intervals.intersection import Interval, IntervalTree

import sys

def parse_gene_coordinate(infile):
    header=1
    for line in open(infile):
        if header:
            header-=1
        else:
            cols = line.strip().split('\t')
            geneid = cols[2]
            chr= cols[0]
            start=cols[1]
            end=cols[6]
            start = int(start)
            end = int(end)
            yield chr, Interval(start, end, value={'geneid':geneid})


def find_overlap(genes, qfile):
    # CHROM  POS     ID      REF     ALT
    
    r=open('out','w')
    for line in open (qfile):
        cols=line.strip().split("\t")
        chrom=cols[0]
        start=int(cols[3])
        end=int(cols[4])
        if cols[6]=="+":
            up3k_start=start-3001
            up3k_end=start-1
            down1k_end=end+1001
            down1k_start=end+1
        elif cols[6]=="-":
            up3k_end=end+3001
            up3k_start=end+1
            down1k_start=start-1001
            down1k_end=start-1
        if chrom in genes:
            if genes[chrom].find(start,end):
                for hit in genes[chrom].find(start,end):
                    
                    r.write(hit.value['geneid']+"\t"+str(hit.start)+"\t"+str(hit.end)+"\t"+"gene_body"+"\t"+line)
            elif genes[chrom].find(up3k_start,up3k_end):
                for hit in genes[chrom].find(up3k_start,up3k_end):
                    
                    r.write(hit.value['geneid']+"\t"+str(hit.start)+"\t"+str(hit.end)+"\t"+"up3k"+"\t"+line)
            elif genes[chrom].find(down1k_start,down1k_end):
                for hit in genes[chrom].find(down1k_start,down1k_end):
                    
                    r.write(hit.value['geneid']+"\t"+str(hit.start)+"\t"+str(hit.end)+"\t"+"down1k"+"\t"+line)
    r.close()

def main():
    infile = sys.argv[1]
    qfile = sys.argv[2]

    genes = {}
    for chr, gene in parse_gene_coordinate(infile):
        if chr not in genes:
            genes[chr] = IntervalTree()

        genes[chr].insert_interval(gene)

    find_overlap(genes, qfile)

if __name__=="__main__":
    main()