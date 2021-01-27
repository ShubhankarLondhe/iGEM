from Bio import Entrez
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
import random as rnd
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Alphabet.IUPAC import ambiguous_dna
import numpy as np

table = CodonTable.ambiguous_dna_by_id[1]

#Entrez.email = "nishant.baruah@students.iiserpune.ac.in" 
#handle = Entrez.efetch(db="nucleotide",id="U00096", rettype="fasta", retmode="text")

s=('A','T','C','G')
handle=open(r"D:\\iGEM\\Sequence.txt")
seq1 = handle.read()
seq2 = seq1.replace('\n','')
my_seq = MutableSeq(seq2, IUPAC.ambiguous_dna)  #making the sequence Mutable

new_seq=my_seq

#print(str(new_seq))
"""
for i in range(1):
        
    #For frameshift mutations
    b=rnd.randint(1,len(my_seq)+1)
    g=rnd.choice(s)                              
    for l in range(len(my_seq)-b):
        my_seq[b+l]=my_seq[b+l-1]
        my_seq[b-1]=g
        
        
    #For point mutations
    a=rnd.randint(0,len(my_seq))
    f=rnd.choice(s)
        
    #For ensuring that mutations happen to different base pairs  
    while(f==my_seq[a]):
        f=rnd.choice(s)
    print(b)
    print(g)
    print(my_seq[a])
    my_seq[a]=f
    print(f)
    print(a)
           """
#seq = my_seq.toseq()
aa_seq = Seq.translate(my_seq)
aa=str(aa_seq)

#print(aa)

handle.close()

proteins = aa.split('*')
print(proteins)
