import numpy as np
import pandas as pd
#from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

#Parse 16S Seq
for rna in SeqIO.parse("vibrio.fasta", "fasta"):
    #print(rna.seq)
    #print(rna.id)
    print(len(rna))

#Probe Seq
Probe = Seq("AGGCCACAACCTCCAAGTAG") 

#convert probe string to mutable object
mutable_probe = Probe.tomutable()
global MutNum

#generate the random mutation position
def RandomPosition(Probe, MutateNum):
    rp = np.random.choice(range(len(Probe)), MutateNum, replace=False)
    RP = rp.tolist()
    RP.sort() #sort the random position in convenient of checking
    return RP

#generate the random mutation base
def RandomBase(ExistBase):
    Base = ["A","T","C","G"]
    Base.remove(ExistBase)
    rb = np.random.choice(Base)
    return rb

#generate the finall random mutaion Probe
def RandomProbe(probe, MutateNum):
    #Probe Seq
    Probe = Seq("AGGCCACAACCTCCAAGTAG")
    #convert probe string to mutable object
    mutable_probe = Probe.tomutable()
    global Pos
    Pos = RandomPosition(probe, MutateNum)
    #print(Pos)
    for p in range(len(Pos)): #single position in mutaition positions 
        RB = mutable_probe[Pos[p]]
        mutable_probe[Pos[p]] = RandomBase(RB)
    return mutable_probe
    
#MutateNum = 4
#RandomProbe(mutable_probe, MutateNum)

#Define output frame
f = pd.DataFrame()

#Match Probe to rRNA
AbsPos = 778
def CountMatch(probe):
    #print(Probe)
    Count = dict()
    
    for i in range(len(rna) - len(Probe)):
        match = 0
        for k in range(len(Probe)):
            if rna[i + k] == probe[k]:
                match += 1
        Count[i] = match
        #print(Count)
    MaxCount = max(Count, key=Count.get)
    if MaxCount == AbsPos:
        
        global f 
        f = f.append({'MutNum':MutNum, "MatchNum":Count[MaxCount], 'MutPosition':Pos, 'MutProbe':str(MutProbe)}, ignore_index=True)
    #print(MaxCount)
    #print(Count[MaxCount])
    #print(f)
    return MaxCount
            

for MutNum in range (1,15):
    print(MutNum)
    count = 0
    while MutNum<11:
        MutProbe = RandomProbe(mutable_probe, MutNum).toseq()#not 1 or 20!
        r = MutProbe.reverse_complement()
        if AbsPos == CountMatch(r):
            count += 1
        if count >= 10:
            break 

f.to_csv("Probe.csv", encoding = "utf-8")

    
        
    