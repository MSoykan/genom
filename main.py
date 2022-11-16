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
    nucName = "ACGT"
    kmerScore=1.0
    for x in range (500-k):
        curKmer=line[x:x+10]
        for y in range(10):
            kmerIndex=curKmer[y]
            nuc = nucName.find(kmerIndex)
            #print(type(kmerIndex))
            #print("KMER INDEX IS :",kmerIndex," ITS enumerateion is =" , nuc)
            #print("ZEYN ZEYN ZEYN",prof[y][nuc])
            kmerScore=kmerScore*prof[y][nuc]
        #print(kmerScore)
        if(kmerScore>score):
            index=x
    return (index,score);

def iterate_randomized_search(dna,k,prof,motif_list):
    new_motif_List = []
    for i in range(len(dna)):
        newLine=dna[i]
        most_probab_kmer_stats= find_most_probable_kmer_in_line(newLine, k, prof)
        if(most_probab_kmer_stats[1]==0):
            most_probab_kmer=motif_list[i]
        if(most_probab_kmer_stats[0]!=0 && most_probab_kmer_stats[])
        #print("MOST PROBAB KMER IN LINE :" +most_probab_kmer );
        new_motif_List.append(most_probab_kmer);
    print("NEW MOTIFS====", new_motif_List, "NEW SCORE===", score(new_motif_List))
    return new_motif_List


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
#profileArr=profile(motif_List)

profileArr = profile(motif_List)
for iter in range(10):
    print("PROFİLE aRRaY İS ", profileArr)
    new_motifffs=iterate_randomized_search(Lines, k, profileArr)
    motif_List=new_motifffs
    print("OLD MOTIFS===", profileArr, "OLD SCORE====", score(motif_List))
    prof=profile(new_motifffs)
    profileArr=prof
    #print("NEW MOTIFS====", prof, "NEW SCORE===", score(new_motifffs))
    print("IS PROF CHaNGING???", prof[0])
print(profileArr)


#calculate profile from motifs







