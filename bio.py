
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

infile=open("input.txt","r")
sequence=infile.readlines()
count=0
clss=[]
asize=[]
for i in range(0,len(sequence)):
    current_seq=sequence[i]
    part=current_seq[0:19]
    count+=1
    print("Sequence: "+ part+ " Length:"+str(len(current_seq)))
    
    try:
        infile=open('dna_lab6_'+str(count)+'.xml', "r")
        print("Using saved file.")
    except FileNotFoundError:
        print("Perfoming online BLAST search")
        handle = NCBIWWW.qblast("blastn", "nt", current_seq)
        with open('dna_lab6_'+str(count)+'.xml', "w") as out_handle:    
            out_handle.write(handle.read())
            handle.close()

    blast_infile= open('dna_lab6_'+str(count)+'.xml', "r")
    records = NCBIXML.parse(blast_infile)
    blast_record = next(records)
    asize.append(blast_record.alignments[0].length)
    if "Vitis" in blast_record.alignments[0].title:
        clss.append(1)
    else:
        clss.append(0)
    blast_infile.close()
print(clss)
print(asize)
plt.hist(asize, normed=False, bins=50)
plt.ylabel('Count')
plt.xlabel('Sequence Length')
plt.title('Histogram of Matching Sequence Lengths')

a = np.array(clss)
unique, counts = np.unique(a, return_counts=True)
print(dict(zip(unique, counts)))

x=[1,0]
y=[3,7]
fig=plt.figure()
ax=fig.add_subplot(111)
ax.bar(x,y,width=1,color=['green','blue'])
ax.set_xticks([0,1])
plt.show()

fig3=plt.figure()
plt.boxplot(asize)    
plt.ylabel('Count')
plt.xlabel('Sequence Length')
plt.title('Histogram of Matching Sequence Lengths')
        
  
   
