import random


def score(motifs):
    score = 0
    for i in range(len(motifs[0])):
        motif = ''.join(motifs[j][i] for j in range(len(motifs)))
        score += min([HammingDistance(motif,homogenous*len(motif)) for homogenous in 'ACGT'])
    return score

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
def find_most_probable_kmer_in_line(line, k, prof):
    score = 0.0
    index = 0
    kmerScore=1.0
    for x in range (500-k):
        curKmer=line[x:x+10]
        for y in range(10):
            kmerIndex=curKmer[y]
            print("KMER INDEX IS :"+kmerIndex+" ITS enumerateion is ="+(kmerIndex in enumerate('ACGT')))
            kmerScore=kmerScore*prof[y][kmerIndex in enumerate('ACGT')]
        print(kmerScore)
        if(kmerScore>score):
            index=x
    return line[x:x+10]
def profile_most_probable_kmer(dna, k, prof):
    dna = dna.splitlines()
    nuc_loc = {nucleotide:index for index,nucleotide in enumerate('ACGT')}
    motif_matrix = []
    max_prob = [-1, None]
    for i in range(len(dna)):
        motif_matrix.append(max_prob)
    for i in range(len(dna)):
        for chunk in window(dna[i],K):
            current_prob = 1
            for j, nuc in enumerate(chunk):
                current_prob*=prof[j][nuc_loc[nuc]]
            if current_prob>motif_matrix[i][0]:
                motif_matrix[i] = [current_prob,chunk]
    return list(list(zip(*motif_matrix))[1])

def profile_most_probable_kmer(dna, k, prof):
    newMotifs= []


# Using readlines() to read from txt dna file.-> returns an array which has every line as element
file1 = open('dna.txt', 'r')
Lines = file1.readlines()
k = 10

count = 0
# Strips the newline character
motif_List = []
print("motif list is :",motif_List)
for line in Lines:
    count += 1
    line = line.upper(); # Mutationları upper yaptık
    line = line.strip(); # /n karakterleri ile boşlukları kırptık
    motifStart=random.randint(0,(500-k)); # k ile değişecek. işlemler girilecek k üzerinden olmalı.
    # select random motifs from lines
    motif_List.append(line[motifStart:(motifStart+k)])
    print(line)
    # select random motifs from lines
profileArr=profile(motif_List)
print(profileArr);

for line in Lines:
    print("MOST PROBAB KMER IN LINE" + find_most_probable_kmer_in_line(line, k,profileArr));
print("motif list is :",motif_List)

#calculate profile from motifs







