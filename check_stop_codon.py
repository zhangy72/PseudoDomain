#!/usr/bin/python
# this script is to check whether there are stop codons 
# in the frame which generate hmmer hits.
import sys
def ReadFaaFile(inFileName):
  faaFile = open(inFileName, 'Ur')
  seqDict = {}
  seqName = ''
  seq = ''
  for line in faaFile:
    if not line.strip():
      continue
    if line[0] == '>':
      if seqName:
        seqDict[seqName] = seq
        seq = ''
      seqName = line.rstrip()[1:]
    else:
      seq += line.rstrip()
  seqDict[seqName] = seq
  faaFile.close()
  return seqDict

def ReadHmmerFile(inFileName, seqDict, outFileName):
  hmmerFile = open(inFileName, 'UR')
  outFile = open(outFileName, 'w')
  for line in hmmerFile:
    items = line.rstrip().split()
    seqName = items[0]
    seqBegin = int(items[6])
    seqEnd = int(items[7])
    # there is stop codon.
    stopCodonNum = seqDict[seqName][seqBegin-1:seqEnd].count('X')
    #print seqName, seqDict[seqName]
    outFile.write('%s %d\n' % (line.rstrip(), stopCodonNum))
  hmmerFile.close()
  outFile.close()
  
if len(sys.argv) < 4:
  sys.stderr.write('Usage: %s <faa file> <hmmer result file> <output file>\n' % sys.argv[0])
  sys.exit()
seqDict = ReadFaaFile(sys.argv[1])
ReadHmmerFile(sys.argv[2], seqDict, sys.argv[3])

      

