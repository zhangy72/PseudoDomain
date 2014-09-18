#!/bin/bash
# this script runs fasta file against its parent's family
# using hmmer3 and hmmframe 
if [ $# -ne 2 ];then
  echo "Usage: <fasta file> <hmmframe output file>"
  exit
fi
cat /dev/null >$2
cat $1 | while read line
do
  seq_name=`echo $line | awk '{print substr($1,2)}'`
  domain=`echo $line | awk '{print $2}'`
  read line
  seq=$line
  echo '>'$seq_name >${seq_name}_${hmm}.fna
  echo $seq >>${seq_name}_${hmm}.fna
  if [ ! -e ${domain}.hmm ];then
    ./get_hmm.sh -a $domain >$domain.hmm
  fi
  ./run_hmmframe.sh -m $domain.hmm -s ${seq_name}_${hmm}.fna -o ${seq_name}_${hmm}_hmmframe.fna
  cat ${seq_name}_${hmm}_hmmframe.fna >>$2
  rm ${seq_name}_${hmm}.fna
  rm ${seq_name}_${hmm}_hmmframe.fna
  rm $domain.hmm 
done 
