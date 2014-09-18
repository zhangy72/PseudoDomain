#!/bin/bash
# Copyright (c) 2011 Yuan Zhang, Yanni Sun.
# You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
# run hmmframe. 

# input: 
# -m: hmm file;
# -f: fasta file;
# -e: error model option. 

# output:
# -o: output file.

usage() {
  echo "./run_hmmframe.sh -m <Pfam HMM file> -s <fasta file> -o <output file> [options] 
  Options:
    -h:  show this message
    -e:  use self-trained model" >&2
}

hmm=
fasta=
out=
error_model=0
while getopts "hem:s:o:" OPTION
do
  case $OPTION in
    h)
      usage
      exit 1
      ;;
    m)
      hmm=$OPTARG
      ;;
    s)
      fasta=$OPTARG
      ;;
    o)
      out=$OPTARG     
      ;;
    e)
      error_model=1
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

if [ $error_model -eq 0 ];then
  ./hmmframe $hmm $fasta $out 0
else
  ./hmmframe $hmm $fasta $out 1
fi 
