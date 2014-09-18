#!/bin/bash
# this file is a pipeline to run selected pfam domains against the compressed reads to see which family is transcribed
if [ $# -ne 2 ];then
  echo "Usage: $0 <hmm file> <fasta file>"
  exit
fi
hmmsearch --cut_ga $1 $2 | grep -E '^Accession:|^>>|^\ *[1-9]+ !|^\ *[1-9]+ \?' | awk '
BEGIN {
  detail = ""
}
{
  if ($1 == ">>") {
    seq=$2
  } else if ($1 == "Accession:") {
    hmm_name = $2
  } else {
    evalue = $5
    score = $3
    model_begin = $7
    model_end = $8
    align_begin = $13
    align_end = $14
    print seq, substr(hmm_name, 1, 7), score, evalue, model_begin, model_end, align_begin, align_end 
  }
}'
