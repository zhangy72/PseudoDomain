Last updated: 9/21/2013

Installation
===================
1. Untar the source file called "pseudodomain.tar.gz".
2. g++ compiler is required in your Unix/Linux system. To install component bin files of PseudoDomain, run the Makeme file using the following command:
make
This will generate the executable file hmmframe for frameshift detection. 
3. Make sure all .sh files and hmmframe are executable in your environment. If not you can use the command:
chmod 755 *.sh
chmod 755 hmmframe
4. Python is required in your Unix/Linux system and in the path. 

Run PseudoDomain
===================
To run PseudoDomain pipeline, use the following command:

./PseudoDomain.sh -m <Pfam HMM file> -f <fasta file> -o <output file> [other options]

Options:
  -h:  show this message
  -c:  domain coverage threshold (default: 0.5)

The hmm file should be in HMMER3.0's hmm file format. It is recommended you download the Pfam-A.hmm file from Pfam ftp. The nucleotide sequence file should be in fasta format.
 
Output
===================
The output file specifies the following fields of the identified processed pseudogene:
[input sequence name] [protein domains separated by comma] [beginning position] [ending position] [stop codon number] [frameshift number]

License
============
Copyright (C) 2012 Yuan Zhang, Yanni Sun.
You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.

