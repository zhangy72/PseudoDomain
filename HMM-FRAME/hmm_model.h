#ifndef HMM_MODEL_H_
#define HMM_MODEL_H_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>


using namespace std;
class HmmModel { // All the numbers are scores. No probabilities.
public:
	void set_length(const int length) {length_ = length;}
	void set_accession(const char* accession) {accession_ = string(accession);}
	void set_acession(string accession) {accession_ = accession;}
	void set_name(const char* name) {name_ = string(name);}
	void set_name(const string name) {name_ = name;}
	void set_nj(const float nj) {nj_ = nj;}
	const float nj() const {return nj_;}
	const int length() const {return length_;}
	const string accession() const {return accession_;}
	const string name() const {return name_;}
	// Static const variables
	static const float kInf_; // The possibly lowest score
	static const int kAlphaSize_; // The number of alphabet, i.e. the number of possible valid amino acid characters
	static const int p7P_N, p7P_E, p7P_C, p7P_B, p7P_J;
	static const int p7P_LOOP, p7P_MOVE;
	static const int M_STATE, I_STATE, D_STATE, G_STATE, B_STATE, E_STATE, J_STATE, N_STATE, C_STATE;
	static const int kTransitionNum_;
	static const int kSpecialStateNum_; // There are 4 special states
	static const int M2M, M2I, M2D, I2M, I2I, D2M, D2D, G2M, M2G, G2G, B2M, M2E;
	static const int kBeginTransitionNum_; // Transitions at the beginning stage of the HMM model
	static const int B2M1, B2I0,B2D, I02M1, I02I0;

	vector<vector<float> > match_emission_; // Emission scores for all match states. Match states are from 0 to length_. M0 is for background. 
											//Amino acids are from 0 to 20
										   // The last emssion score is for 'X'
	vector<vector<float> > insert_emission_; // Emission scores of 20 amino acids for all insertion states.
	vector<vector<float> > transition_; // Transition index begins with 1. 0 in the vector is vacant. Transition scores in the order of M->M,M->I,M->D,I->M,I->I
	vector<float> begin_transition_;    // B->M1,B->I0,B->D,I0->M1,I0->I0
	vector<vector<float> > special_transition_;			// transition scores for 5 special states, p7P_B, p7P_C, p7P_E, p7P_N and p7P_J.
														// Each state includes p7P_LOOP and p7P_MOVE

private:
	int length_;       
	string accession_; // Accession number of the Pfam domain
	string name_; // Name of the Pfam domain
	float nj_;	// Expected number of uses of J states
};
bool ConstructHmmModel(ifstream& in_file, HmmModel& hmm, float G_state_prob); // Set the parameters which have not been set in constructor
void TokenizeString(const string& line, vector<string>& tokens, const string& delimiters);// Tokenize the input line and output a vector of all separate strings
float FLogSum(float *vec, int n);
float FMax(float *vec, int n);
int p7_hmm_CalculateOccupancy(const HmmModel& hmm, float *mocc, float *iocc); // Written by HMMER3.0 group. Thanks!
#endif