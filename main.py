import random

def HammingDistance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError('Undefined for sequences of unequal length.')
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def generateProfile(Motifs):
    arr = [[0]*4]*10
    print(arr,"before")
    for line in Motifs:
        for line2 in Motifs[k]:
            print()

def profile(motifs):
    prof = []
    for i in range(len(motifs[0])):
        col = ''.join([motifs[j][i] for j in range(len(motifs))])
        prof.append([float(col.count(nuc))/float(len(col)) for nuc in 'ACGT'])
    return prof

# Using readlines() to read from txt dna file.-> returns an array which has every line as element
file1 = open('dna.txt', 'r')
Lines = file1.readlines()
k = 10

count = 0
# Strips the newline character
temp_List = []
print("motif list is :",temp_List)
for line in Lines:
    count += 1
    line = line.upper(); # Mutationları upper yaptık
    line = line.strip(); # /n karakterleri ile boşlukları kırptık
    motifStart=random.randint(0,(500-k)); # k ile değişecek. işlemler girilecek k üzerinden olmalı.
    temp_List.append(line[motifStart:(motifStart+k)])
    print(line)
    # select random motifs from lines

print("motif list is :",temp_List)

profileArr=profile(temp_List)
print(profileArr);
# select random motifs from lines

