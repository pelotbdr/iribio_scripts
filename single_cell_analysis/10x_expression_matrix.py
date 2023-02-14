#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: beaudier
"""

# =============================================================================
# Import des packages et des arguments de la ligne de commande
# =============================================================================

import sys
import re
import csv

maxInt=sys.maxsize
while True:
	try:
		csv.field_size_limit(maxInt)
		break
	except OverflowError:
		maxInt=int(maxInt/10)
        
        
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--data',nargs='+')
parser.add_argument('--reference',type=str)
parser.add_argument('--out',type=str)
parser.add_argument('--minsize',type=int)  
args=parser.parse_args()

data=args.data
ref=args.reference
outname=args.out
sizelim=args.minsize


# =============================================================================
# Functions
# =============================================================================
#Convert isoform to gene name
def iso_transf(gene):
    no=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','x','y','z']
    points=([m.start() for m in re.finditer('\.',gene)])
    if(len(points)>1):
        b=gene[:points[-1]]
    else:
        b=gene
    if(b[-1] in no):
        b=b[:-1]
    return b

#Split string by symbol
def cutter(predec,symbol):
    sortie={}
    pos=([m.start() for m in re.finditer(symbol,predec)])
    adj=0
    for i in range(len(pos)):
        if(predec[adj:pos[i]] not in sortie):
            sortie[predec[adj:pos[i]]]=1
        else:
            sortie[predec[adj:pos[i]]]+=1
        adj=pos[i]+1
    if(len(pos)>0):
        if(predec[adj:] not in sortie):
            sortie[predec[adj:]]=1
        else:
            sortie[predec[adj:]]+=1
    else:
        if(predec not in sortie):
            sortie[str(predec)]=1
    sortie2=[]
    for elem in sortie:
        for i in range(sortie[elem]):
            if(len(elem)>1):
                sortie2.append(elem)
    return(sortie2)


# =============================================================================
# Write expression matrix
# =============================================================================
genes_ids={}
with open(ref,mode='r') as handle:
    for line in handle:
        if('>' in line):
            pos=([m.start() for m in re.finditer('gene=',line)])
            if(len(pos)>0):
                isoform=line[1:pos[0]-1]
                gene=iso_transf(isoform)
                genes_ids[gene]=0
handle.close()


###Instead of a pandas dataframe, a new line is directly written for each barcode

#Write header
with open(outname+'.csv',mode='w') as out:
    for elem in genes_ids:
        out.write(','+elem)
    out.write('\n')
    
    #Go through each file to get the count per gene and write it
    for elem in data:
        print(elem)
        with open(elem,mode='r') as handle:
            handle2=csv.reader(handle,delimiter='\t')
            for row in handle2:
                #Clean the str and split it by comma
                erf=row[1].replace('[','').replace(']','').replace(' ','').replace("'","")
                genes=cutter(erf,',')
                
                #Apply the minimum size limit
                if(len(genes)>sizelim):
                    genes2={}
                    #Convert isoform to gene name
                    for gen in genes_ids:
                        genes2[gen]=0
                    for gen in genes:
                        genes2[iso_transf(gen)]+=1
                    
                    #Write barcode name with file name
                    out.write(row[0]+'_'+elem)
                    #Write gene counts 
                    for gen in genes_ids:
                        out.write(','+str(genes2[gen]))
                    out.write('\n')
        handle.close()
out.close()