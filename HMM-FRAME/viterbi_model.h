#ifndef VITERBI_MODEL_H_
#define VITERBI_MODEL_H_
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <stack>
#include "hmm_model.h"
#include "DNA_seq.h"

using namespace std;
int Codon2AminoAcidIndex(const char* codon);
void Viterbi(HmmModel& hmm, DNASeq& dna ,ofstream& out_file);
int p7_ReconfigLength(HmmModel& hmm, int seq_length);
float HiddenCodon2AminoAcidIndex(char* codon, int missing_position, vector<float> background, 
                                 vector<float> match_emission, char& missing_base);
void TraceBack(const vector<vector<vector<int> > >& flag,
                           const vector<vector<int> >& E_state_flag, 
                           const vector<vector<vector<char> > >& deleted_base,
                           int L, int M, const char* seq, string& out_seq,
                           int& error_num, vector<vector<int> >& seq_error_pos, 
                           vector<vector<int> >& state_error_pos);
string StripOutputSeq(const string& seq);
void OutputPath(stack<string>& stack);
#endif