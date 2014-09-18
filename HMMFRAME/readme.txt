Version
===================
This software was last updated 3/16/2012. To report bugs, please send emails to zhangy72@msu.edu.

Installation
===================
1. Untar the source file called "HMMFRAME.tar.gz".
2. g++ compiler is required in your Unix system. To install HMM-FRAME, run the Makeme file using the following command:
make

Run HMM-FRAME
===================
The installation will generate the bin file: hmmframe. To run hmmframe, use the following command:
./run_hmmframe.sh -m <hmm file> -s <nucleotide fasta file> -o <output file> [-e]

The hmm file can contain multiple hmm models and should be in HMMER3.0's hmm file format. These files can be downloaded from Pfam ftp. The nucleotide sequence file should in be fasta format. The default reference model uses 0.0007 for non-homopolymer regions and 0.0044 for homopolymer regions. When the option "-e" is set, then the error model will shift to our self-trained model. For example, to run the fasta file test.fa, against nifH family model file nifH.hmm using the reference error model:

./run_hmmframe -m nifH.hmm -s test.fa -o test.hmmframe
 
Output
===================
The output of HMMFRAME is in fasta format. It keeps the sequence name but adds additional information such as Pfam family, Viterbi score, number of errors and error positions in both the original sequence and the model. Also the sequences are corrected based on the predicted errors. 
For example, for the following input sequence against nifH family:

>GB4XUSJ05FQ71N
CACCCGTCTGATGCTTCACTCCAAAGCTCAAACCACCGTACTACACTTAGCTGCTGAACGCGGTGCAGTAGAAGACTTAGAACTCCACGAAGTAATGTTGACCGGTTTCCGTGGCGTTAAGTGCGTAGAATCTGGTGGTCCAGAACCCGGTGTAGGTTGCGCCGGTCGTGGTATCATCACCGCCATTAACTTCTTAGAAGAAAACGGCGCTTTACCAAGACCTAGACTTCGTATCCTACGACGTATTGGGTGACGTTGTATGTGGTGGTTTCGCTATGCCTATCCGTGAAGGTAAAGCACAAGAAATCTACATCGTTACC

The output of HMMFRAME is:

>GB4XUSJ05FQ71N hmm_name=nifH score=216.586 error_num=1 seq_insertion_pos=211, seq_deletion_pos= state_insertion_pos=116, state_deletion_pos=
ACCCGTCTGATGCTTCACTCCAAAGCTCAAACCACCGTACTACACTTAGCTGCTGAACGCGGTGCAGTAGAAGACTTAGAACTCCACGAAGTAATGTTGACCGGTTTCCGTGGCGTTAAGTGCGTAGAATCTGGTGGTCCAGAACCCGGTGTAGGTTGCGCCGGTCGTGGTATCATCACCGCCATTAACTTCTTAGAAGAAAACGGCGCiTTACCAAGACCTAGACTTCGTATCCTACGACGTATTGGGTGACGTTGTATGTGGTGGTTTCGCTATGCCTATCCGTGAAGGTAAAGCACAAGAAATCTACATCGTTACC

The sequence name is GB4XUSJ05FQ71N. The query hmm family is nifH family. The alignment score is 216.586 bits. There is one error.

License
===================
Copyright (C) 2011 Yuan Zhang, Yanni Sun.
You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.
If you find HMM-FRAME useful, please cite the following paper:
Y. Zhang and Y. Sun, "HMM-FRAME: accurate protein domain classification for metagenomic sequences containing frameshift errors", BMC Bioinformatics, vol. 12, no. 1, p. 198, 2011.

