from Bio import Entrez
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
import random as rnd
from Bio.Alphabet import IUPAC
def translate(a):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  #this table dictionary is pre-created
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein=""
    if len(a)%3==0:
        for i in range(0,len(a),3):
            codon=a[i:i+3]
            if(codon in ('ATA','ATC','ATT','ATG','ACA','ACC','ACG','ACT','AAC','AAT','AAA','AAG','AGC','AGT','AGA','AGG','CTA','CTC','CTG','CTT','CCA','CCC','CCG','CCT','CAC','CAT','CAA','CAG','CGA','CGC','CGG','CGT','GTA','GTC','GTG','GTT', 'GCA','GCC','GCG','GCT','GAC','GAT','GAA','GAG','GGA','GGC','GGG','GGT', 'TCA','TCC','TCG','TCT','TTC','TTT','TTA','TTG', 'TAC','TAT','TAA','TAG','TGC','TGT','TGA','TGG','ATA')):
                protein+=table[codon]
            else:protein+='*'
    return(protein)



# Downloading the Escherichia coli str. K-12 substr. MG1655, complete genome from NCBI

Entrez.email = "nishant.baruah@students.iiserpune.ac.in" 
handle = Entrez.efetch(db="nucleotide",id="U00096", rettype="fasta", retmode="text")

s=('A','T','C','G')
a=handle.read()
my_seq = MutableSeq(a, IUPAC.unambiguous_dna)  #making the sequence Mutable
#for i in range(1):
        
    #For frameshift mutations
      #  b=rnd.randint(1,len(my_seq)+1)
      #  g=rnd.choice(s)                              
      #  for l in range(len(my_seq)-b):
      #      my_seq[b+l]=my_seq[b+l-1]
      #      my_seq[b-1]=g
        
        
    #For point mutations
        #a=rnd.randint(0,len(my_seq))
        #f=rnd.choice(s)
        
    #For ensuring that mutations happen to different base pairs  
        #while(f==my_seq[a]):
        #    f=rnd.choice(s)
        #print(b)
        #print(g)
        #print(my_seq[a])
        #my_seq[a]=f
        #print(f)
        #print(a)
myn_seq=my_seq.toseq()
mynm_seq=str(myn_seq)                                     


print(translate(mynm_seq))   
