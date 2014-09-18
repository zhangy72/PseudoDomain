all:
	g++ -O3 -w -o hmmframe HMM-FRAME/hmmframe.cpp HMMFRAME/viterbi_model.cpp HMMFRAME/DNA_seq.cpp HMMFRAME/hmm_model.cpp
	g++ -O3 -w -o DNA2Protein HMM-FRAME/DNA2Protein.cpp
