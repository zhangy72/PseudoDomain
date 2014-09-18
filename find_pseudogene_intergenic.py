#!/usr/bin/python
import sys
def ReadHits(inFileName, domainLengthDict):
  hits = []
  inFile = open(inFileName, 'Ur')
  for line in inFile:
    items = line.rstrip().split()
    seqName = items[0]
    domainAcc = items[1][0:7]
    if domainLengthDict[domainAcc] < 50:
      continue
    score = float(items[2])
    evalue = float(items[3])
    domainBegin = int(items[4]) 
    domainEnd = int(items[5])
    hitBegin = int(items[6]) # relative begin of hit in the region.
    hitEnd = int(items[7]) # relative end of hit in the region.
    absoluteBegin = 3 * (hitBegin - 1) + 1
    absoluteEnd = 3 * hitEnd
    stopCodonNum = int(items[8])
    hits.append([seqName, [domainAcc], domainBegin, domainEnd, absoluteBegin, absoluteEnd, False, stopCodonNum, score, evalue])
  inFile.close()  
  # sort hits according to sequence name and beginning position in sequence.
  hits.sort(key=lambda hit: (hit[0], hit[4]))
  return hits
# show ordered hits. for debug purpose.
def ShowHits(hits, outFileName):
  # write original hits in order to files.
  hitFile = open(outFileName, 'w')
  for hit in hits:
    hitFile.write('%s\n' % hit)
  hitFile.close()

def RemoveRedundantHits(hits, overlapRate):
  if len(hits) < 2:
    return hits
  nonRedundantHits = []
  currentHit = []
  sameDomain = 0
  diffDomain = 0
  for hit in hits:
    if currentHit:
      if currentHit[0] == hit[0]:
        minLength = min(currentHit[5]-currentHit[4]+1, hit[5]-hit[4]+1)
        # two neighboring hits have significant overlap.
        overlap = Overlap(currentHit[4], currentHit[5], hit[4], hit[5])
        if overlap > overlapRate * minLength:
          #print currentHit, hit, overlap
          if currentHit[1][0] == hit[1][0]: sameDomain += 1
          else: diffDomain += 1
          currentHitEvalue = currentHit[-1]
          hitEvalue = hit[-1]
          if currentHitEvalue > hitEvalue:
            currentHit = hit       
        # two neighboring hits do not have significant overlap.
        else:
          nonRedundantHits.append(currentHit)
          currentHit = hit     
      else:
        nonRedundantHits.append(currentHit)
        currentHit = hit
    else:
      currentHit = hit
  # handle the last hit.
  nonRedundantHits.append(currentHit)
  #print 'Samedomain', sameDomain, 'diffDomain', diffDomain
  return nonRedundantHits 

# this function removes singleton hits with poor domain coverage.
def DomainFiltration(outHits, domainLengthDict, domainCoverage):
  filteredOutHits = []
  for hit in outHits:
    if hit[6]:
      filteredOutHits.append(hit)
    else:
      # for singleton hit, we require domain coverage.
      domainAcc = hit[1][0]
      domainLength = domainLengthDict[domainAcc]
      hitDomainLength = hit[3] - hit[2] + 1
      if hitDomainLength >= domainLength * domainCoverage:
        filteredOutHits.append(hit)
  return filteredOutHits

# this function merges neighboring hits.    
def MergeHits(hits, genomicGap):
  outHits = []
  if not hits:
    return outHits 
  elif len(hits) == 1:
    return hits
  else:
    currentHit = [] 
    for hit in hits:
      if currentHit:
        # hits from the same genomic sequence.
        if currentHit[0] == hit[0]:
          if Distance(currentHit[4], currentHit[5], hit[4], hit[5]) <= genomicGap:
            # merge two neighboring hits.
            hit[1] += currentHit[1] # combine all the domains.
            hit[4] = currentHit[4]
            hit[6] = True # this hit has been merged.
            hit[7] += currentHit[7] # add the stop codons
          else:
            outHits.append(currentHit)
        # now change the genomic sequence.  
        else:
          outHits.append(currentHit)
      currentHit = hit
  # handle the last hit.  
  outHits.append(currentHit)
  return outHits       

# overlap is defined as 0.
def Distance(begin1, end1, begin2, end2):
  if max(begin1, begin2) - min(end1, end2) <= 1:
    return 0
  else:
    return max(begin1, begin2) - min(end1, end2) - 1    

def Overlap(begin1, end1, begin2, end2): 
  if max(begin1, begin2) - min(end1, end2) >= 1:
    return 0
  else:
    return min(end1, end2) - max(begin1, begin2) + 1   

def OutputHits(outHits, outFileName):
  outFile = open(outFileName, 'w')
  if outHits:
    for hit in outHits:
      domainAccs = ''
      for domainAcc in hit[1]:
        domainAccs = domainAccs + domainAcc + ','
      outFile.write('%s %s %d %d %d\n' % (hit[0], domainAccs, hit[4], hit[5], hit[7]))
  outFile.close()

def GetDomainLength(inFileName):
  domainLengthFile = open(inFileName, 'Ur')
  domainLengthDict = {}
  for line in domainLengthFile:
    domainAcc = line.rstrip().split()[0]
    domainLength = int(line.rstrip().split()[1])
    domainLengthDict[domainAcc] = domainLength
  return domainLengthDict

# main function.
if len(sys.argv) != 6:
  sys.stderr.write('Usage: %s <hmmer file> <genomic gap in bp> <output file>\
 <domain length file> <domain coverage>\n' % sys.argv[0]) 
  sys.exit()
hitFileName = sys.argv[1]
genomicGap = int(sys.argv[2])
outFileName = sys.argv[3]
domainLengthDict = GetDomainLength(sys.argv[4])
domainCoverage = float(sys.argv[5])
overlapRate = 0.5

hits = ReadHits(hitFileName, domainLengthDict)
#ShowHits(hits, 'original_hits.list')
nonRedundantHits = RemoveRedundantHits(hits, overlapRate)
#ShowHits(nonRedundantHits, 'nonredundant_hits.list')
outHits = MergeHits(nonRedundantHits, genomicGap)
#ShowHits(outHits, 'out_hits.list')
filteredOutHits = DomainFiltration(outHits, domainLengthDict, domainCoverage)
OutputHits(filteredOutHits, outFileName)
