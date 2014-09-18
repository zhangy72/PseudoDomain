all:
	g++ -O3 -w -o hmmframe HMM-FRAME/hmmframe.cpp HMM-FRAME/viterbi_model.cpp HMM-FRAME/DNA_seq.cpp HMM-FRAME/hmm_model.cpp
	g++ -O3 -w -o DNA2Protein HMM-FRAME/DNA2Protein.cpp
