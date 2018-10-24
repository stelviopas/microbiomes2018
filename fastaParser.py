#simple example of parsing FASTA
from Bio import SeqIO
import numpy as np
from numpy import median

counter = 0
#set min,max Read to the first one to start
minReadName = maxRead = ""
minRead = maxRead = len(next(SeqIO.parse("/Users/study/Downloads/GreenShake.fasta", "fasta")))
gc_content = 0
lengthSum = 0
np_length = []
np_Ids = []
#counters for distribution
#since the shortest sequence was about 90 bp
lessThan100 = 0
lessThan200 = 0
lessThan500 = 0
lessThan700 = 0
lessThan1000 = 0
lessThan2000 = 0
lessThan4000 = 0
lessThan10000 = 0
#since the longest sequence was about 13000 bp
lessThan14000 = 0

for seq_record in SeqIO.parse("/Users/study/Downloads/GreenShake.fasta", "fasta"):
    #we use counter to track the number of reads, if indexing was broken somewhere
    counter = counter + 1
    #sum up all length of reads
    lengthSum = lengthSum + len(seq_record)
    #longest Read
    #shortest read
    if len(seq_record) < minRead:
        minRead = len(seq_record)
        minReadName = seq_record.id
    if len(seq_record) > maxRead:
        maxRead = len(seq_record)
        maxReadName = seq_record.id
    #print(seq_record.id)
    #print(seq_record.seq)
    #print(len(seq_record))
    np_length.append(len(seq_record))
    np_Ids.append(seq_record.id)
    #counters for distribution
    ##########################
    seq = len(seq_record)
    if seq < 100:
        lessThan100 += 1
    elif seq < 200:
        lessThan200 += 1
    elif seq < 500:
        lessThan500 += 1
    elif seq < 700:
        lessThan700 += 1
    elif seq < 1000:
        lessThan1000 += 1
    elif seq < 2000:
        lessThan2000 += 1
    elif seq < 4000:
        lessThan4000 += 1
    elif seq < 10000:
        lessThan10000 += 1
    elif seq < 14000:
        lessThan14000 += 1
    ##########################
    read = seq_record.seq
    for i in range(len(read) - 1):
        if read[i] == "G":
            gc_content += 1
        elif read[i] == "C":
            gc_content += 1

#Total number of reads
print ("Total number of reads", counter)
#shortest Read
print(minReadName + " has length", minRead)
#longestRead
print(maxReadName + " has length", maxRead)
#mean or average Length
print ("Mean length is ", lengthSum/counter)
#print (np)
#find a median automatically using numpy
medianReadValue = median(np_length)
index = np_length.index(medianReadValue)
medianReadID = np_Ids[index]
#print(index)
print("A Median read " + medianReadID + " has length", medianReadValue)
print ("Total length is ", lengthSum)
#distribution stuff
list = [lessThan100,
lessThan200,
lessThan500,
lessThan700,
lessThan1000,
lessThan2000,
lessThan4000,
lessThan10000,
lessThan14000]
print (list)
sum = 0
for a in list:
    print(100 * (a / counter), " % ")
    print((round(100*(a/counter))), " % ")

print ("GC content is ", (100*gc_content/lengthSum), " %")

