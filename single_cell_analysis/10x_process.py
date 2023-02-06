#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: beaudier
"""

### 10x single cell data processing
#Input: Read file n°1 containing 26 nucleotides: Barcode (0-15) and UMI (16-25)
#Output: 2 TSV file with column 1 containing a barcode name and column 2 
#the barcode filtered content (either reads ids or gene names)

# =============================================================================
# Initialisation
# =============================================================================

import sys
import argparse
import csv
import re

maxInt=sys.maxsize
while True:
	try:
		csv.field_size_limit(maxInt)
		break
	except OverflowError:
		maxInt=int(maxInt/10)

print('Parsing arguments...')
parser = argparse.ArgumentParser()
parser.add_argument('--read1',type=str)
parser.add_argument('--sam',type=str)
parser.add_argument('--out',type=str)
args=parser.parse_args()

read1=args.ids
sam=args.sam
out=args.out

# =============================================================================
# Filter Barcodes and UMIs
# =============================================================================

uniques=[] 		#Sequences list
ids=[]			#Reads ids list
with open(read1,mode='r') as handle:
    inc=1
    for line in handle :
        if(inc==1):
            ids.append(str(line)[str(line).find('@') + 1:str(line).rfind(' N')])
        if(inc==2):
            uniques.append(str(line))
        inc+=1
        if(inc==5):
            inc=1
handle.close()


umis={}								#List of UMIs in a barcode
ids_assigned={}						#List of reads ids in a barcode
for i in range(len(uniques)):
	if(uniques[i][0:16] not in umis):
		umis[uniques[i][0:16]]=[]
		ids_assigned[uniques[i][0:16]]=[]
		umis[uniques[i][0:16]].append(uniques[i][16:26])
		ids_assigned[uniques[i][0:16]].append(ids[i])
	else:
		if(uniques[i][16:26] not in umis[uniques[i][0:16]]):
			umis[uniques[i][0:16]].append(uniques[i][16:26])
			ids_assigned[uniques[i][0:16]].append(ids[i])

with open(out+"_filtered.tsv",mode='w') as handle:
	for key in ids_assigned:
		handle.write(str(key)+'\t'+str(ids_assigned[key])+'\n')
handle.close()


# =============================================================================
# Connect read id to mapped gene
# =============================================================================

mapped={}   #Genes per mapped ids
with open(sam,mode='r') as handle:
    done={}
    handle2=csv.reader(handle,delimiter='\t', quotechar='"')
    for row in handle2:
        if('@' not in row[0]):
            if('*' not in row[2]):
                if(row[0] not in done):
                    mapped[row[0]]=row[2]
                    done[row[0]]=0
handle.close()

ids={}
with open(out+'_filtered.tsv',mode='r') as handle:
    handle2=csv.reader(handle,delimiter='\t', quotechar='"')
    for row in handle2:
        ids[row[0]]=[]
        clean=row[1].replace('[','').replace(']','').replace(' ','').replace("'","")
        #Ids are separated by commas so we differentiate them this way
        #The classic split function was not used because it was too slow 
        pos=([m.start() for m in re.finditer(',',clean)])
        
        adj=0       #Keep count of position in string
        
        #Go on each id one by one and check if it corresponds to a mapped gene
        for i in range(len(pos)):
            if(clean[adj:pos[i]] in mapped):
                ids[row[0]].append(mapped[clean[adj:pos[i]]])
            adj=pos[i]+1
            
        #If only one gene in the barcode
        if(len(pos)==0):
            if(clean in mapped):
                ids[row[0]].append(mapped[clean])
        #If not, there is still the last id behind the last comma
        else:
            if(clean[adj:] in mapped):
                ids[row[0]].append(mapped[clean[adj:]])
handle.close()

del mapped

with open(out+'_mapped.tsv',mode='w') as handle:
	for key in ids:
		handle.write(str(key)+'\t'+str(ids[key])+'\n')
handle.close()