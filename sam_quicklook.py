#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: beaudier
"""

### Produce global informations on the results of the mapping contained in the sam file
#Gives mapping % and read lengths

import csv
import sys
import matplotlib.pyplot as plt

name=sys.argv[1]

maxInt=sys.maxsize
while True:
	try:
		csv.field_size_limit(maxInt)
		break
	except OverflowError:
		maxInt=int(maxInt/10)


read_total=0
mapped_subtotal=0
control_subtotal=0
tailles=[]
tailles_control=[]
done={}
with open(name, mode='r',encoding="ISO-8859-1") as handle:
    handle2=csv.reader(handle,delimiter='\t', quotechar='"')
    for row in handle2:
        if('@' not in row[0]):
            if(row[0] not in done):
                read_total+=1
                done[row[0]]=0
                if('*' not in row[2]):
                    if(row[2]!='REFERENCE'):
                        mapped_subtotal+=1
                        tailles.append(len(row[9]))
                    else:
                        control_subtotal+=1
                        tailles_control.append(len(row[9]))

handle.close()
perc_mapped=round((mapped_subtotal/read_total)*100,2)
perc_control=round((control_subtotal/read_total)*100,2)
print(str(perc_mapped)+'% mapped on the Reference')
print(str(perc_control)+'% mapped on the Internal DNA/RNA control')

labels='Non-mapped reads','Mapped reads','Internal DNA/RNA control'
sizes=[100-perc_mapped-perc_control,perc_mapped,perc_control]
colors='Red','lightgreen','lightblue'
explode=(0,0.1,0)
figs, axs = plt.subplots(2)
axs[0].pie(sizes, explode=explode,labels=labels, autopct='%1.1f%%', colors=colors, shadow=False, startangle=190)
axs[0].axis('equal')
axs[1].hist(tailles,bins=50,color='lightgreen',label='Reference',alpha=0.8)
axs[1].hist(tailles_control,bins=50,color='lightblue',label='Internal DNA/RNA control',alpha=0.8)
plt.legend()
plt.title('Mapped reads lengths')
plt.tight_layout()
plt.show()
