#!/bin/bash
# Copyright (c) 2012 Yuan Zhang, Yanni Sun.
# You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
# pipeline of PseudoDomain. 

# input: 
# -m: hmm file;
# -f: genomic sequence file in fasta format;
# -c: domain coverage threshold, default: 0.5.

# output:
# -o: identified processed pseudogene in the following format:
# [genomic sequence name][a list of domains separated by comma]
# [beginning position of pseudogene][ending position of pseudogene] [stop codon number]
# [frameshift number]

usage() {
  echo "./PseudoDomain.sh -m <Pfam HMM file> -f <genomic sequence fasta file> -o <output file> [options] 
  Options:
    -h:  show this message
    -c:  domain coverage threshold (default: 0.5)" >&2
    #-r:  consider reverse complement of the input sequence" >&2
}

hmm=
fasta=
coverage=0.5
consider_reverse_complement=0

while getopts "hm:f:c:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    m)
      hmm=$OPTARG
      ;;
    f)
      fasta=$OPTARG
      ;;
    c)
      coverage=$OPTARG
      ;;
    o)
      out=$OPTARG     
      ;;
    esac
done

if [ "$hmm" == "" ];then
  echo "Please specify the hmm file!"
  usage
  exit
fi

if [ "$fasta" == "" ];then
  echo "Please specify the input fasta file!"
  usage
  exit
fi

if [ "$out" == "" ];then
  echo "Please specify the output file!"
  usage
  exit
fi

# generate an empty output file.
cat /dev/null >$out

# keeps the length of hmms.
if [ ! -f hmm_acc_length.list ];then 
cat $hmm | awk '
/^ACC/ {
  acc = $2 
}
/^LENG/ {
  len = $2
  print substr(acc, 1, 7), len
}' >hmm_acc_length.list
fi

# run hmmer and generate all raw hits.
# hmmer results from all domains. stop codon are annotated.
cat /dev/null >${fasta}_stopcodon.hmmer
if [ ${consider_reverse_complement} -eq 0 ];then
  ./DNA2Protein 1-3 $fasta $fasta
  for i in {1..3}
  do
    ./hmmer3_pipeline_ga.sh $hmm $fasta.frame${i} \
      >${fasta}_frame${i}.hmmer
    python ./check_stop_codon.py $fasta.frame${i} \
      ${fasta}_frame${i}.hmmer \
      ${fasta}_frame${i}_stopcodon.hmmer
    cat ${fasta}_frame${i}_stopcodon.hmmer\
      >>${fasta}_stopcodon.hmmer
  done
else
  ./DNA2Protein 1-6 $fasta $fasta
  for i in {1..6}
  do
    ./hmmer3_pipeline_ga.sh $hmm $fasta.frame${i} \
      >${fasta}_frame${i}.hmmer
    python ./check_stop_codon.py $fasta.frame${i} \
      ${fasta}_frame${i}.hmmer \
      ${fasta}_frame${i}_stopcodon.hmmer
    cat ${fasta}_frame${i}_stopcodon.hmmer \
      >>${fasta}_stopcodon.hmmer
  done
fi
# clean some tmp files.
rm $fasta.frame?
rm ${fasta}_frame?.hmmer

# find pseudogenes from raw hits.
python ./find_pseudogene_intergenic.py ${fasta}_stopcodon.hmmer \
  60 ${fasta}_stopcodon.pseudogene \
  hmm_acc_length.list $coverage

rm ${fasta}_frame?_stopcodon.hmmer
rm ${fasta}_stopcodon.hmmer
rm hmm_acc_length.list

# find frameshifts.
# generate data set of hmmframe.
python ./generate_hmmframe_dataset.py $fasta ${fasta}_stopcodon.pseudogene \
  ${fasta}_pseudogene_for_hmmframe.fasta
# run hmmframe on the data set.
./run_hmmframe_on_pseudogene.sh ${fasta}_pseudogene_for_hmmframe.fasta \
  ${fasta}_pseudogene_for_hmmframe.hmmframe
# combine pseudogene result and hmm result
python ./output_pseudogene.py ${fasta}_stopcodon.pseudogene \
  ${fasta}_pseudogene_for_hmmframe.hmmframe $out  
  
rm ${fasta}_pseudogene_for_hmmframe.fasta
rm ${fasta}_pseudogene_for_hmmframe.hmmframe
rm ${fasta}_stopcodon.pseudogene
