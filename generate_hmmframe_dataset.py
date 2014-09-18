#!/usr/bin/python
# this program accepts a pseudogene file and generate a fasta file
# to run hmmframe.
import sys
# read the sequences in fasta file
def ReadFastaFile(inFileName):
  inFile = open(inFileName, 'Ur')
  seqName = ''
  seq = ''
  seqDict = {}
  for line in inFile:
    strippedLine = line.rstrip()
    if strippedLine:
      if strippedLine[0] == '>':
        if seqName:
          seqDict[seqName] = seq
        seqName = strippedLine[1:]
        seq = ''
      else:
        seq += strippedLine
  # add the last sequence.
  seqDict[seqName] = seq
  return seqDict      

# read pseudogene file and output fasta file.
def ReadPseudogeneFile(inFileName, seqDict, outFileName):
  inFile = open(inFileName, 'Ur')
  outFile = open(outFileName, 'w')
  extension = 30 # extend the original subsequence.
  for line in inFile:
    if line.rstrip():
      items = line.rstrip().split()
      seqName = items[0]
      domainAccs = items[1]
      seqBegin = int(items[2]) 
      seqEnd = int(items[3]) 
      assert seqName in seqDict
      outSeq = seqDict[seqName][seqBegin-1-extension:seqEnd+extension]
      domains = domainAccs.split(',')
      domainSet = set(domains)
      for domain in domainSet:
        if domains.count(domain) > 1:
          outFile.write('>%s_%d_%d %s\n' % (seqName, seqBegin, seqEnd, domain)) 
          outFile.write(outSeq + '\n')
  inFile.close()
  outFile.close()

# main function.
if len(sys.argv) < 4:
  sys.stderr.write('Usage: %s <fasta file> <pseudogene file> <output file>\n' % sys.argv[0])
  sys.exit()
seqDict = ReadFastaFile(sys.argv[1])
ReadPseudogeneFile(sys.argv[2], seqDict, sys.argv[3])

