import numpy as np
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import Align
from Bio.SubsMat import MatrixInfo
import time
import random

#Uses mutate() to create array of mutants, n is the variants per generation
def variants(seq, v, med):
    dna_arr = []
    if med=='r':
        for j in range(v):
            mut_dna = mutate_rich(seq)
            dna_arr.append(mut_dna)
    
    elif med=='p':
        for j in range(v):
            mut_dna = mutate_poor(seq)
            dna_arr.append(mut_dna)
            
    return dna_arr

#Takes the best sequence and mutates it many times to create variants, dna_arr is a list
#Nested inside variants()
def mutate_rich(seq):           #Mutations in Rich medium       
    random.seed()
    dna = list(seq)
    #print(dna)
    #print(seq)
    i=3
    #print(len(dna))
    while i<len(dna):            #Don't want to mutate start codon therefore start with 2
        if dna[i] == 'N':
            dna.pop(i)
        else :
            z = random.random()
            if z<8e-4:                         #Actual mutation rate in rich medium is 8 x 10^-6 /nt
                x = random.gauss(0,1)          #For Checking mutation type     
                if x<=1.64:                    #Transitions
                    if dna[i]=='A':
                        dna[i]='G'
                    elif dna[i]=='T':
                        dna[i]='C'
                    elif dna[i]=='G':
                        dna[i]='A'
                    else:
                        dna[i]='T' 
                                    
                elif 1.64<x<=1.96:              #Transversion
                    y = random.uniform(0,1)     #Choose which transversion
                    if y<=0.5:                                   
                        if dna[i]=='A':
                            dna[i]='T'
                        elif dna[i]=='T':
                            dna[i]='A'
                        elif dna[i]=='G':
                            dna[i]='C'
                        else:
                            dna[i]='G'      
                        
                    else:
                        if dna[i]=='A':
                            dna[i]='C'
                        elif dna[i]=='T':
                            dna[i]='G'
                        elif dna[i]=='G':
                            dna[i]='T'
                        else:
                            dna[i]='A'  
                                
                elif 1.96<x<=3.9:               #Deletion 
                    dna.pop(i)                  
                else:                           #Addition
                    s = random.choice('ATGC')  
                    dna.insert(i+1,s)
        i+=1    
        
    mut_dna = str(''.join(dna))
    return mut_dna

def mutate_poor(seq):          #Mutations in poor medium
    random.seed()
    dna = list(seq)
    #print(dna)
    #print(seq)
    i=3
    #print(len(dna))
    while i<len(dna):            #Don't want to mutate start codon therefore start with 2
        if dna[i] == 'N':
            dna.pop(i)
        else :
            z = random.random()
            if z<1e-5:                            #Actual mutation rate in poor medium is 10^-7 /nt
                x = random.gauss(0,1)          #For Checking mutation type     
                if x<=0.88:                    #Transversion
                    y = random.uniform(0,1)    #Choose which transversion
                    if y<=0.5:                                   
                        if dna[i]=='A':
                            dna[i]='T'
                        elif dna[i]=='T':
                            dna[i]='A'
                        elif dna[i]=='G':
                            dna[i]='C'
                        else:
                            dna[i]='G'      
                        
                    else:
                        if dna[i]=='A':
                            dna[i]='C'
                        elif dna[i]=='T':
                            dna[i]='G'
                        elif dna[i]=='G':
                            dna[i]='T'
                        else:
                            dna[i]='A'      
                                    
                elif 0.88<x<=1.96:          #Transitions
                    if dna[i]=='A':
                        dna[i]='G'
                    elif dna[i]=='T':
                        dna[i]='C'
                    elif dna[i]=='G':
                        dna[i]='A'
                    else:
                        dna[i]='T'      
                                
                elif 1.96<x<=3.9:              #Deletion 
                    dna.pop(i)                  
                else:                          #Addition
                    s = random.choice('ATGC')  
                    dna.insert(i+1,s)
        i+=1    
        
    mut_dna = str(''.join(dna))
    return mut_dna

#Takes Many sequences and chooses the best one
def choose(dna_arr):
    prot_arr = to_prot(dna_arr)   #Call function to_prot(), may give since prot_arr is not defined as an list
    scr = []
    #seq = []              To return multiple sequences with same score
    #indx = []             Find python function to return all max value indices in a list 
    
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = MatrixInfo.blosum62
    for p in range(len(prot_arr)):
        alignments = aligner.align(prot_best, prot_arr[p])
        scr.append(alignments.score)
    
    indx = np.argmax(scr)
    seq = dna_arr[indx]
    return seq           

#Takes an array of many DNAs, gives protein sequences
#Nested inside choose()
def to_prot(dna_arr):
    prot_arr = []                    
    for m in range(len(dna_arr)):
        if len(dna_arr[m])%3!=0:      #To make the sequence length divisible by three for translation
            seq = list(dna_arr[m])
            q = len(dna_arr[m])
            for i in range(q%3):
                seq.append('N')
            dna_arr[m] = str(''.join(seq))    
        
        my_seq = MutableSeq(dna_arr[m], IUPAC.ambiguous_dna)
        prot = Seq.translate(my_seq,stop_symbol='')
        prot_arr.append(str(prot))
    return prot_arr

#Main Body of the Code------------------------------------------------------------------

#dna_ori = 'ATGAATATCCAGATCGGCGAGCTTGCCAAGCGCACCGCATGCCCGGTGGTGACCATTCGCTTCTACGAACAAGAAGGGCTGTTGCCGCCGCCGGGCCGCAGCCGGGGGAATTTTCGCCTGTATGGCGAGGAGCACGTGGAGCGCTTGCAGTTCATTCGTCACTGCCGGTCTCTGGATATGCCGTTGAGCGACGTACGGACCTTATTGAGTTACCGGAAGCGGCCCGACCAGGATTGCGGTGAAGTCAATATGCTCTTGGATGAGCACATCCGTCAGGTCGAATCTCGGATCGGAGCTTTGCTCGAACTGAAGCACCATTTGGTGGAACTGCGCGAAGCCTGTTCTGGTGCCAGGCCCGCCCAATCGTGCGGGATTCTGCAGGGACTGTCGGACTGCGTGTGTGATACGCGGGGGACCACCGCCCATCCAAGCGACTAG'
dna_ori = 'ATGATATATATATAT' #input('Enter original gene :',dna_ori)
z = len(dna_ori)
#print(z)

#Find More Efficient Protein
dna_best = 'ATGGCGCGCGCGCGC' #input('Enter gene of better efficiency :',dna_best)
my_seq = MutableSeq(dna_best, IUPAC.unambiguous_dna)
prot_b = Seq.translate(my_seq,stop_symbol='')
prot_best = str(prot_b)


gens = int(input('How many generations do you want to run it for ?'))
v = int(input('How many variants per generation ?'))
n = int(input('Selection after how many generations ?'))

while True:
    med = input('Which medium do you want the bacteria in? Rich(r) or Poor(p) ')
    a = med.lower()
    if a=='r' or a=='p':
        break    
    else:
        print("Enter 'r' or 'p'")

dna = dna_ori
dna_arr = []

for k in range(1,gens+1):        #We want generation to start from 1 and not 0
    if dna == dna_best:
        break;
    elif k%n == 0:
        dna_arr = variants(dna,v, med)
        dna = choose(dna_arr)
    elif k%n == 1:
        dna_arr = variants(dna,v, med)
    else:
        if med=='r':
            for k in range(v):
                dna_arr[k] = mutate_rich(dna)
        elif med=='p':
            for k in range(v):
                dna_arr[k] = mutate_poor(dna)
        
if k%n != 0:
    dna = choose(dna_arr)

if k<gens:
    print('It took',k,'generations for original gene',dna_ori,'to become the more efficient gene',dna_best)
else :
    print('Gene surviving after',gens,'generations is',dna,'\n','Original gene was',dna_ori)
