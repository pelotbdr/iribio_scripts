#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: beaudier
"""
# =============================================================================
# Packages import and argument parsing
# =============================================================================

import csv
import sys
import matplotlib.pyplot as plt
import statistics
import argparse
import os

maxInt=sys.maxsize
while True:
	try:
		csv.field_size_limit(maxInt)
		break
	except OverflowError:
		maxInt=int(maxInt/10)


parser = argparse.ArgumentParser()
parser.add_argument('--sam',nargs='+')
parser.add_argument('--names',nargs='+')
parser.add_argument('--out',type=str)
args=parser.parse_args()

sam_files=args.sam
names={}
for i in range(len(args.names)):
    names[sam_files[i]]=args.names[i]
outname=args.out


# =============================================================================
# Functions
# =============================================================================

### N50 stat: value at which half of the data is contained within longer reads
def n50(numbers):
    total=sum(numbers)
    numbers.sort()
    count=0
    for elem in numbers:
        count+=elem
        if(count>total/2):
            return elem

### Processing of total amount of Insertions, Deletions and Soft clips in sequence
def cigar_process(seq):
    nbrs=['1','2','3','4','5','6','7','8','9','0']
    point=0
    soft=[]
    insert=[]
    delet=[]
    for i in range(len(seq)):
        if(seq[i] not in nbrs):
            if(seq[i]=='S'):
                soft.append(seq[point:i])
            if(seq[i]=='I'):
                insert.append(seq[point:i])
            if(seq[i]=='D'):
                delet.append(seq[point:i])
            point=i+1
    total=0
    for elem in soft:
        total+=int(elem)
    for elem in insert:
        total+=int(elem)
    for elem in delet:
        total-=int(elem)
    
    return total
               



def sam_process(filename):
    #Reads sizes
    sizes=[]
    #Read 5' end
    start_sites=[]
    #Read 3' end
    end_sites=[]
    #Number of complete reads (with a 3% margin of error allowed)
    complete_reads=0
    #Number of total reads
    read_count=0
    #Number of mapped reads
    mapped_reads=0
    with open(filename,mode='r',encoding="ISO-8859-1") as handle:
        handle2=csv.reader(handle,delimiter='\t')
        length=0
        for row in handle2:
            if('@' in row[0]):
                if('DNACS' not in row[1]):
                    if('LN' in row[2]):
                        #Reference sequence length
                        length=int(row[2].replace('LN:','').replace('\n',''))    
            else:
                read_count+=1
                if(row[2] not in ['*','DNACS']):
                    mapped_reads+=1
                    sizes.append(int(len(row[9])))
                    start_sites.append(int(row[3]))
                    
                    if(length*0.97<len(row[9])<length*1.03):
                        complete_reads+=1
                    
                    correct=cigar_process(row[9])
                    end=int(row[3])+correct
                    if(end<length):
                        end_sites.append(end)
                    elif(end<length*1.03):
                        end_sites.append(length)           
    handle.close()
    
    return sizes,start_sites,end_sites,complete_reads,read_count,mapped_reads,length


# =============================================================================
# Results printing
# =============================================================================

sizes={}
start_sites={}
end_sites={}

os.system('mkdir '+outname)
with open(outname+'/stats.txt',mode='w') as handle:
    for elem in sam_files:
        print('Reading '+names[elem]+'.......')
        handle.write('File '+names[elem]+'\n')
        
        a,b,c,complete_reads,read_count,mapped_reads,length=sam_process(elem)
        sizes[elem]=a
        start_sites[elem]=b
        end_sites[elem]=c
        
        handle.write(str(mapped_reads)+' reads mapped out of '+str(read_count)+'\n')
        handle.write('Mean size: '+str(statistics.mean(a))+'\nMedian size: '+ 
                     str(statistics.median(a))+'\nN50:'+str(n50(a))+'\n\n\n')
handle.close()

nfiles=len(sam_files)

### Log and non-log distribution of read sizes
fig,ax=plt.subplots(nfiles,figsize=(6*nfiles,4*nfiles))
i=0
for elem in sizes:
    if(nfiles==1):
        ax.hist(sizes[elem],bins=60,range=(0,length*1.1))
        ax.title.set_text(names[elem])
    else:
        ax[i].hist(sizes[elem],bins=60,range=(0,length*1.1))
        ax[i].title.set_text(names[elem])
    i+=1
plt.xlabel('Read size')
plt.yscale('log',base=10)
plt.tight_layout()
plt.savefig(outname+'/log_distribution.png',dpi=150)

plt.clf()

fig,ax=plt.subplots(nfiles,figsize=(6*nfiles,4*nfiles))
i=0
for elem in sizes:
    if(nfiles==1):
        ax.hist(sizes[elem],bins=60,range=(0,length*1.1))
        ax.title.set_text(names[elem])
    else:
        ax[i].hist(sizes[elem],bins=60,range=(0,length*1.1))
        ax[i].title.set_text(names[elem])
    i+=1

plt.xlabel('Read size')
plt.tight_layout()
plt.savefig(outname+'/non-log_distribution.png',dpi=150)

plt.clf()


### Reads end sites plot
fig,ax = plt.subplots(len(sam_files),figsize=(5*nfiles,3*nfiles))

start_points={}
end_points={}
i=0
for elem in sam_files:
    for tip in start_sites[elem]:
        if(tip not in start_points):
            start_points[tip]=1
        else:
            start_points[tip]+=1
    for tip in end_sites[elem]:
        if(tip not in end_points):
            end_points[tip]=1
        else:
            end_points[tip]+=1
    
    if(nfiles==1):
        ax.scatter(start_points.keys(),start_points.values(),s=5,color='blue')
        ax.scatter(end_points.keys(),start_points.values(),s=5,color='red')
        ax.set_yscale('log',base=10)
        ax.title.set_text(names[elem])
    else:
        ax[i].scatter(start_points.keys(),start_points.values(),s=5,color='blue')
        ax[i].scatter(end_points.keys(),start_points.values(),s=5,color='red')
        ax[i].set_yscale('log',base=10)
        ax[i].title.set_text(names[elem])
    i+=1
plt.tight_layout()
plt.xlabel('Read length',fontsize=20)      
plt.savefig(outname+'/reads_ends.png',dpi=150)