# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 09:52:18 2019

@author: Junyu

Record Mismatch Positions
"""

from Bio import SeqIO
for rna in SeqIO.parse("vibrio.fasta", "fasta"):
    #print(rna.seq)
    #print(rna.id)
    print(len(rna))
    
Seq = "AGACTTGGAGGTTGTGGCCT" #mutate 1st C to A 2d T to G

    
def CountMatch(Probe):
    print(Probe)
    Count = dict()
    #Mismatch = dict()
    mismatch = []
    for i in range(len(rna) - len(Probe)):
        match = 0
        for k in range(len(Probe)):
            if rna[i + k] == Probe[k]:
                match += 1
            #else: #mismatch 
                #mismatch.append(k) no need to record all position
        
        Count[i] = match
        #Mismatch[i] = mismatch
        
        #print(Mismatch)
    MaxCount = max(Count, key=Count.get)
    for i in range(len(Probe)):
        if rna[MaxCount+i] != Probe[i]:
            mismatch.append(i+1)
        
    
    print(MaxCount)
    print(mismatch)
    return Count
            
CountMatch(Seq)
