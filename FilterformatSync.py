#!/usr/bin/python
#########################################
# File Name: FilterformatSync.py
#
# Purpose: Take as intput a .sync file as generated by popoolation2 (hence
# filtered to base quality and around indel positions). Outputs a list
# of bi-polymorphic sites filtered:
#   for mincount of minor allele in allpop
#   for mincov in each pop
#   for maxcov in each pop
#
#
# Created by: Etienne Loire
#
###########################################



# Define some variables (to pass as arguments later)
MINCOV = 10 # Minimum global coverage
MAXCOV = 100 # Maximum global coverage
MAF = 0.05 # Minor allele frequency
MCTS = 3 # Minimum counts of minor allele

# Define some output types:
KEEPMONO = True # whether or not to output monomorphic positions
OUTPUTMAXCOV= True # whether or not to write SNP at positions with coverage > MAXCOV in a separate file

import sys

if OUTPUTMAXCOV:
    maxcovfile=open(sys.argv[1]+".bigcov_file.sync","w")

# open .sync file
print "Treating "+sys.argv[1]+" file ..."
infile = open(sys.argv[1],"r")
# open out file
outfile = open(sys.argv[1]+".filtered_"+str(MINCOV)+"_"+str(MAXCOV)+"_"+str(MCTS)+".sync", "w")

# compute coverage for a position and a pop
def popcov(counts):
    c = counts.split(":")
    cov=0
    for i in c[:4]:
        cov += int(i)
    return cov

# compute coverage for each base (ATGC) in a position (all pop included)

def basecounts(allcounts):
    # Split all columns and add them to data
    data = []
    for i in allcounts:
        data.append(i.split(":"))
    # Now data countains [[0,1,0,2,0,0],[2,5,4,20,0,0], ...]
    # basecounts will store counts of each base A,T,G,C observed in ALL pop
    basecounts= [0,0,0,0]
    # for each base:
    for i in range(0,4):
    # add for  each pop:
        for pop in data:
            basecounts[i]+=int(pop[i])
    return basecounts


# Main loop: Treat all lines in file:

# Just to keep order of the bases in the sync files
dicobase = {"A":0,"T":1,"C":2,"G":3}
dicopos =  {0:"A",1:"T",2:"C",3:"G"}


for line in infile:
    c = line.split() # split line, sep = tabulation
    # store some infos
    chrom = c[0]
    pos = int(c[1])
    refbase = c[2]
    # compute each base coverage over all populations
    basecov=basecounts(c[3:])
    # filter all parameters:

    flag = True # OK by default for this line but
    # If one pop is not covered enough, pass
    for population in range(3,len(c)):
        if popcov(c[population])<MINCOV:
            flag=False
            continue

    # If it's not a biallelic loci (more than 2 alleles with cov >= MCTS)
    # (or MAF)
    alleles=0
    for ma in basecov: # global coverage of each base which is not ref:
        if ma>=MCTS:
            alleles +=1
    if alleles>2:
        # more than 2 alleles: filterout
        flag=False
    if alleles == 1:
        # Monomorphic (keep ot remove)
        if not KEEPMONO:
            flag=False    # stay in the loop
    if flag:
        #last test: if one pop has a coverage higher than MAXCOV, we need to output the line in a separate file:
        maxcovflag=False
        for population in range(3,len(c)):
            if popcov(c[population])>MAXCOV:
                maxcovflag=True
                continue
        if not maxcovflag:
            outfile.write(line) # keep line because it's okay !
        elif OUTPUTMAXCOV:
            maxcovfile.write(line)

# Close file handler
if OUTPUTMAXCOV:
    maxcovfile.close()
outfile.close()
