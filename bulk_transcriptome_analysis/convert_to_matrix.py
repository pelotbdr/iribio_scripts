#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 15:18:03 2022

@author: beaudier
"""

# =============================================================================
# Import des packages et des arguments de la ligne de commande
# =============================================================================

import sys
import re
import pandas as pd
import csv
import statistics
import os

maxInt=sys.maxsize
while True:
	try:
		csv.field_size_limit(maxInt)
		break
	except OverflowError:
		maxInt=int(maxInt/10)
        
        
import argparse
print('Parsing arguments...')
parser = argparse.ArgumentParser()
parser.add_argument('--model',type=str)
parser.add_argument('--reference',type=str)
parser.add_argument('--data',nargs='+')
parser.add_argument('--out',type=str)  
args=parser.parse_args()

model=args.model
ref=args.reference
names=args.data
outname=args.out


# =============================================================================
# Fonctions pour le script
# =============================================================================
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

def col_index(corr,name):
    return(corr.columns.get_loc(name))

def row_index(corr,name):
    return(corr.index.get_loc(name)) 


# =============================================================================
# Processing du fichier d'alignement pour créer la matrice d'expression
# =============================================================================
print('Reading reference...')
genes={}
ids={}
if(model=='human'):
    with open(ref,mode='r') as handle:
        for line in handle:
            if('>' in line):
                pos=([m.start() for m in re.finditer('gene:ENSG',line)])
                if(len(pos)>0):
                    pos2=([m.start() for m in re.finditer('\.',line[pos[0]:])])
                    gene=iso_transf(line[pos[0]+5:pos[0]+pos2[0]])
                    pos3=([m.start() for m in re.finditer(r'[^\S\n\t]+',line)])
                    transcript=line[1:pos3[0]]
                    ids[transcript]=gene
                    genes[gene]=0
    handle.close()
if(model=='celegans'):
    with open(ref,mode='r') as handle:
        for line in handle:
            if('>' in line):
                pos=([m.start() for m in re.finditer('gene=',line)])
                if(len(pos)>0):
                    isoform=line[1:pos[0]-1]
                    gene=iso_transf(isoform)
                    ids[isoform]=gene
                    genes[gene]=0
    handle.close()
print('Creating matrix...')
all_reads={}
assigned_reads={}

d = pd.DataFrame(0, index=genes.keys(),columns=names)
for elem in names:
    all_r=0
    assigned_r=0
    done={}
    temp={}
    for gen in genes: 
        temp[gen]=0
    for file in os.listdir(elem):
        if('.sam' in file):
            with open(elem+'/'+file,mode='r',encoding="ISO-8859-1") as handle:
                handle2=csv.reader(handle,delimiter='\t')
                for row in handle2:
                    if('@' not in row[0]):
                        if(row[0] not in done):
                            all_r+=1
                            if('*' not in row[2]):
                                assigned_r+=1
                                if('REFERENCE' not in row[2]):
                                    if(row[2] in ids):
                                        gene=ids[row[2]]
                                        a=iso_transf(gene)
                                        temp[a]+=1
                                done[row[0]]=0
            handle.close()
            for key in temp:
                #d.iloc[row_index(d,key),col_index(d,elem)]=int(temp[key])
                d.loc[key,elem]=int(temp[key])
            all_reads[elem]=all_r
            assigned_reads[elem]=assigned_r
print('Writing matrix...')
d.to_csv(outname+'_matrix.csv',sep=',')

# =============================================================================
# Stats sur les fichiers
# =============================================================================
print('Generating statisitcs...')
with open(outname+'_stats.txt',mode='w') as out:
    for elem in names:
        tailles={}
        for file in os.listdir(elem):
            if('.fastq' in file):
                with open(elem+'/'+file,mode='r') as handle:
                    inc=1
                    for line in handle:
                        if(inc==2):
                            tailles[len(line.replace('\n',''))]=0
                        inc+=1
                        if(inc==5):
                            inc=1
                handle.close()
        percent=(assigned_reads[elem]/all_reads[elem])*100
        out.write(elem+'\n')
        out.write(str(all_reads[elem])+' reads with an average length of '+ str(statistics.mean(tailles.keys()))+' bases\n')
        out.write(str(percent)+'% of reads have been mapped on the reference\n')
        out.write('\n\n')
out.close()

print('Done')