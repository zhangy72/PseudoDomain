#ifndef DNASEQ_H_
#define DNASEQ_H_
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>

using namespace std;

class DNASeq {
public:
  DNASeq(const string& seq_name, const string&seq, 
         float non_homo_error_rate, float homo_error_rate);
	string seq_name() const {return seq_name_;}
	string seq() const {return seq_;}
	int length() const {return length_;}
	void set_corrected_seq(string out_seq) {corrected_seq_ = out_seq;}
	string corrected_seq() const {return corrected_seq_;} 
	vector<string> ProteinSeqs(); // Return 6-frame translations of seq_
	bool error_model() const {return error_model_;}
	void set_error_model(bool error_model) {error_model_ = error_model;}
  float non_homo_error_rate() const {return non_homo_error_rate_;}
  float homo_error_rate() const {return homo_error_rate_;}
  string ReverseComplement() const;
	friend int Codon2AminoAcidIndex(const char* codon);
	friend char Codon2AminoAcid(const char* codon);
private:
	string seq_name_; 
	string seq_; 
	int length_;
	string corrected_seq_; // The output sequence after error correction
	bool error_model_;
  float non_homo_error_rate_;
  float homo_error_rate_;
	static const int kCodonNum_;
	static const char* kCodonTable_; // Amino acids for 64 possible codons plus 'X' for codons wtih unexpected characters other than A, C, G and T
	static const int kCodonIndexTable_[65]; // Indices for each amino acid according to their position in hmm files
	static const int kMaxLineLength_;
};
// The index begins with 1 to be consistent with Viterbi.
vector<vector<float> > CalculateErrorScore(string seq, const bool mode,
                                           float non_homo_error_rate, 
                                           float homo_error_rate); 
int Codon2AminoAcidIndex(const char* codon);
char Codon2AminoAcid(const char* codon);

#endif
  