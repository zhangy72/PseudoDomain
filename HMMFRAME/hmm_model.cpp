#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include "hmm_model.h"
using namespace std;
const float HmmModel::kInf_ = -32768;
const int HmmModel::kAlphaSize_ = 20; // The number of alphabet, i.e. the number of possible valid amino acid characters
const int HmmModel::kSpecialStateNum_ = 5;// B, C, E, N and J are 4 special states
const int HmmModel::p7P_B = 0;
const int HmmModel::p7P_C = 1;
const int HmmModel::p7P_E = 2;
const int HmmModel::p7P_N = 3; 
const int HmmModel::p7P_J = 4;
const int HmmModel::p7P_LOOP = 0;
const int HmmModel::p7P_MOVE = 1;
const int HmmModel::M_STATE = 0;
const int HmmModel::I_STATE = 1;
const int HmmModel::D_STATE = 2;
const int HmmModel::G_STATE = 3;
const int HmmModel::B_STATE = 4;
const int HmmModel::E_STATE = 5;
const int HmmModel::J_STATE = 6;
const int HmmModel::N_STATE = 7;
const int HmmModel::C_STATE = 8;

const int HmmModel::kTransitionNum_ = 12;
const int HmmModel::M2M = 0;
const int HmmModel::M2I = 1;
const int HmmModel::M2D = 2;
const int HmmModel::I2M = 3;
const int HmmModel::I2I = 4;
const int HmmModel::D2M = 5;
const int HmmModel::D2D = 6;
const int HmmModel::G2M = 7;
const int HmmModel::M2G = 8;
const int HmmModel::G2G = 9;
const int HmmModel::B2M = 10;
const int HmmModel::M2E = 11;

const int HmmModel::kBeginTransitionNum_ = 5; 
const int HmmModel::B2M1 = 0;
const int HmmModel::B2I0 = 1;
const int HmmModel::B2D = 2;
const int HmmModel::I02M1 = 3;
const int HmmModel::I02I0 = 4;

bool ConstructHmmModel(ifstream& hmm_file, HmmModel& hmm, float G_state_prob) {
  assert(hmm_file);
  const double kXProb = 0.05;
  string line;
  stringstream sstr;
    while (getline(hmm_file, line)) {
      if (line.find_first_not_of(" \t") == string::npos) {
        continue;
      } else {
      sstr.clear();
      sstr.str(line);
      string field;
      sstr >> field;
      if (field == "//") {
        return true; // this file has not been done.
      } else if (field == "LENG") {
        sstr >> field;
        hmm.set_length(atoi(field.c_str()));
        // Note that we have I0 state and 'X'    
        hmm.insert_emission_ = vector<vector<float> >(hmm.length() + 1, vector<float>(HmmModel::kAlphaSize_ + 1));
        // M0 is background emissions        
        hmm.match_emission_ = vector<vector<float> >(hmm.length() + 1, vector<float>(HmmModel::kAlphaSize_ + 1));
        // Transitions start with 1.Transition_[0] is vacant            
        hmm.transition_ = vector<vector<float> >(hmm.length() + 1, vector<float>(HmmModel::kTransitionNum_));         
      } else if (field == "NAME") {
        sstr >> field;
        hmm.set_name(field);
      } else if (field == "ACC") {
        hmm.set_acession(field.substr(0, 7));
      } else if ( field == "COMPO") { // This block reads the three lines after the COMPO line inclusively.
        vector<string> tokens;
        TokenizeString(line, tokens, " \t");
        // Give the average probability for 20 amino acids to residue 'X'
        hmm.match_emission_[0][HmmModel::kAlphaSize_] = log(kXProb);
        for (int i = 0; i < HmmModel::kAlphaSize_; i++) {
          if (tokens[i + 1] == "*") {
            hmm.match_emission_[0][i] = HmmModel::kInf_;
          } else {
        // Background probability
          hmm.match_emission_[0][i] = 0.0f - atof(tokens[i + 1].c_str()); // Convert to log(probability)
          }
        }
        getline(hmm_file, line); // I0 insertion scores
        TokenizeString(line, tokens, " ");
        // Give the average probability for 20 amino acids to residue 'X'
        hmm.insert_emission_[0][HmmModel::kAlphaSize_] = log(kXProb);
        for (int i = 0; i < HmmModel::kAlphaSize_; i++) {
          if (tokens[i] == "*") {
            hmm.insert_emission_[0][i] = HmmModel::kInf_;
          } else {
            hmm.insert_emission_[0][i] = 0.0f - atof(tokens[i].c_str());
          }
        }
        getline(hmm_file, line); // Begin transition line
        hmm.begin_transition_ = vector<float>(HmmModel::kBeginTransitionNum_);
        for (int i = 0; i < HmmModel::kBeginTransitionNum_; i++) { 
        // Only consider the first 5 fields. The last two are for deletion state,
        // which does not exist in the beginning stage.
          hmm.begin_transition_[i] = 0.0f - atof(tokens[i].c_str());
        }

        // assign the values of begin transition to transition[0].
        for (int i = 0; i < HmmModel::kBeginTransitionNum_; ++i) {
          hmm.transition_[0][i] = hmm.begin_transition_[i];
        }
        for (int i = 1; i < hmm.length(); i ++) {
          getline(hmm_file, line); // Match state emissions
          TokenizeString(line, tokens," ");
          hmm.match_emission_[i][HmmModel::kAlphaSize_] = log(kXProb);
          for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
            hmm.match_emission_[i][j] = 0.0f - atof(tokens[j + 1].c_str()); // The first column is the state index                    
          }
          getline(hmm_file, line); // Insertion state emissions.
          TokenizeString(line, tokens," ");
          hmm.insert_emission_[i][HmmModel::kAlphaSize_] = log(kXProb);
          for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
            hmm.insert_emission_[i][j] = 0.0f - atof(tokens[j].c_str());
          }
          getline(hmm_file, line); // Transition scores.
          TokenizeString(line, tokens," ");
          for (int j = 0; j < 7; j++) {
            hmm.transition_[i][j] = 0.0f - atof(tokens[j].c_str());
          }
        }
        // Consider the last block, which has the transition from the last state to END
        getline(hmm_file, line); // Match state emissions
        TokenizeString(line, tokens, " ");
        // Emission score for 'X' should be the minimum value
        hmm.match_emission_[hmm.length()][HmmModel::kAlphaSize_] = log(kXProb);
        for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
          hmm.match_emission_[hmm.length()][j] = 0.0f - atof(tokens[j + 1].c_str());
        }
        getline(hmm_file, line); // Insertion state emissions
        TokenizeString(line, tokens," ");
        hmm.insert_emission_[hmm.length()][HmmModel::kAlphaSize_] = log(kXProb);
        for (int j = 0; j < HmmModel::kAlphaSize_; j++) {
          hmm.insert_emission_[hmm.length()][j] = 0.0f - atof(tokens[j].c_str());
        }
        getline(hmm_file, line); // Transition scores
        TokenizeString(line, tokens," ");
        for (int j = 0; j < 7; j++) {
          if (tokens[j] == "*") {
            hmm.transition_[hmm.length()][j] = HmmModel::kInf_;
          } else {
            hmm.transition_[hmm.length()][j] = 0.0f - atof(tokens[j].c_str());
          }
        }
        // Calculate the scores for G state and END state
        hmm.special_transition_ = vector<vector<float> >(HmmModel::kSpecialStateNum_, vector<float>(2));
        for (int i = 1; i < hmm.length(); i++) {
          hmm.transition_[i][HmmModel::M2G] = log(G_state_prob);
          float vec[2] = {0.0f, hmm.transition_[i][HmmModel::M2G]}; //log(1) = 0.
          float log_sum = FLogSum(vec, 2); // log(1 + G_state_prob)
          hmm.transition_[i][HmmModel::M2M] -= log_sum;
          hmm.transition_[i][HmmModel::M2I] -= log_sum;
          hmm.transition_[i][HmmModel::M2D] -= log_sum;
          hmm.transition_[i][HmmModel::M2G] -= log_sum;
          hmm.transition_[i][HmmModel::G2G] = hmm.transition_[i][HmmModel::M2G];
          hmm.transition_[i][HmmModel::G2M] = log(1 - G_state_prob);
          hmm.transition_[i][HmmModel::M2E] = 0.0f; // No penalty for any M state to END state
        }
        hmm.special_transition_[HmmModel::p7P_E][HmmModel::p7P_MOVE] = 0.0f;   
        hmm.special_transition_[HmmModel::p7P_E][HmmModel::p7P_LOOP] = HmmModel::kInf_;  
        hmm.set_nj(0.0f);
        // Transition between BEGIN state and match states
        for (int i = 1; i <= hmm.length(); ++i) {
          hmm.transition_[i - 1][HmmModel::B2M] = 0.0f;
        }
      } else {
	    // useless lines, skip.
        continue;
      }
    } 
  } 
  // the file is done.
  return false;
  
    //float Z = 0.0f;
    //float *occ = (float *)malloc(sizeof(float ) * hmm.length() + 1);
 //   if ((p7_hmm_CalculateOccupancy(hmm, occ, NULL)) != 0) {
    //    cout << "B2M score calculation error!" << endl;
    //    hmm_file.close();
    //    exit(1);
    //}
 //   for (int k = 1; k <= hmm.length(); k++) {
    //  Z += occ[k] * float(hmm.length() - k + 1);
 //   }
 //   for (int k = 1; k <= hmm.length(); k++) {
    //    if (occ[k] <= 0) {
    //        hmm.transition_[k - 1][HmmModel::B2M] = HmmModel::kInf_;
    //    } else {
    //        hmm.transition_[k-1][HmmModel::B2M] = log(occ[k] / Z);
    //    }
    //  //cout << hmm.transition_[k - 1][HmmModel::B2M] << endl;
 //   } 
}
//Given vector of log p's; return log of summed p's
float FLogSum(float *vec, int n) 
{
  int x;
  float max, sum;
  max = FMax(vec, n);
  sum = 0.0;
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      sum += exp(vec[x] - max);
  sum = log(sum) + max;
  return sum;
}
float FMax(float *vec, int n) {
  int i;
  float best;
  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] > best) best = vec[i];
  return best;
}
void TokenizeString(const string& str, vector<string>& tokens, const string& delimiters) {
  tokens.clear(); // Initialize tokens because tokens.push_back() will be used.
  string::size_type last_pos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, last_pos);
  while (string::npos != pos || string::npos != last_pos)  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(last_pos, pos - last_pos));
    // Skip delimiters.  Note the "not_of"
    last_pos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, last_pos);
  }
}
int p7_hmm_CalculateOccupancy(const HmmModel& hmm, float *mocc, float *iocc) {
  int k;
  const int length = hmm.length();
  mocc[0] = 0.;                                /* no M_0 state */
  mocc[1] = exp(hmm.begin_transition_[HmmModel::M2I]) + exp(hmm.begin_transition_[HmmModel::M2M]);   /* initialize w/ 1 - B->D_1 */
  for (k = 2; k <= length; k++){
      mocc[k] = mocc[k-1] * (exp(hmm.transition_[k-1][HmmModel::M2M]) + exp(hmm.transition_[k-1][HmmModel::M2I])) +
          (1.0-mocc[k-1]) * exp(hmm.transition_[k-1][HmmModel::D2M]);  
  }
  if (iocc != NULL) {
      iocc[0] = exp(hmm.begin_transition_[HmmModel::M2I] - hmm.begin_transition_[HmmModel::I2I]);
    for (k = 1; k <= length; k++)
        iocc[k] = mocc[k] * exp(hmm.transition_[k][HmmModel::M2I] - hmm.transition_[k][HmmModel::I2I]);
  }
  return 0;
}
