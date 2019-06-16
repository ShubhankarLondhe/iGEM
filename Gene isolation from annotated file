import re
f = open(r"D:\\iGEM\\Genes.txt")
g = f.read()
h = g.replace('\n','')
raw_seq = re.findall(r'[A-Z]+', h) #Keeps only uppercase characters
print(len(raw_seq))
#print(raw_seq)

seq = [x for x in raw_seq if len(x) >7]  #Removes all sequences whose length is less than 7 (Id Tags till length 7)
        
print(len(seq)) 
print(seq)

f.close()
