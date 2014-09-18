#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <string>
#include <stack>
#include "viterbi_model.h"
using namespace std;
// In this function all indices should begin with 0 to be consistent
void Viterbi(HmmModel& hmm, DNASeq& dna, ofstream& out_file) {
  const int M = hmm.length(); // M is the number of states
  const int L = dna.length();
  char* seq = new char[L+2];
  seq[0] = '$';
  strncpy(seq+1, dna.seq().c_str(), dna.length());
  // Use reference error model. Seq error scores begin with 1
  vector<vector<float> > seq_error_score = CalculateErrorScore(dna.seq(), dna.error_model(), 
    dna.non_homo_error_rate(), 
    dna.homo_error_rate()); 
  // Store deleted bases in DP. Note that different may lead to position conflicts. 
  // So 0, 1 and 2 of the first dimension represented a
  // deleted base after seq[i], seq[i - 1] and seq[i - 2].
  vector<vector<vector<char> > > deleted_base(3, vector<vector<char > >(M+1, vector<char>(L+1, '*')) );
  p7_ReconfigLength(hmm, L / 3); // Configure some special parameters for hmm model.
  //The following matrix for DP begins with 1. The first position is unused.
  vector<vector<float> > mat_M(M + 1, vector<float>((L + 1), HmmModel::kInf_)); // Score matrix for match staes
  vector<vector<float> > mat_I(M + 1, vector<float>((L + 1), HmmModel::kInf_)); // Score matrix for insertion states
  vector<vector<float> > mat_D(M + 1, vector<float>((L + 1), HmmModel::kInf_)); // Score matrix for deletion states
  vector<vector<float> > mat_G(M + 1, vector<float>((L + 1), HmmModel::kInf_)); // Score matrix for G states
  vector<vector<float> > mat_X(4, vector<float>(L + 1, HmmModel::kInf_)); // B, C, E and N states
  vector<vector<vector<int> > > flag(4,vector<vector<int> >(M + 1, vector<int>(L + 1, -1))); // Flags for M, I, D and G states
  //  E_state_flag[i][0] stores the value of j which gives largest END score; 
  // E_state_flag[i][1] indicates whether E state leads to largest score. 
  vector<vector<int> >E_state_flag(L + 1, vector<int>(2, -1)); 
  vector<float> background = hmm.match_emission_[0]; // Background emission scores
  // Do some initialization
  for (int i = 0; i <= L; ++i) {
    mat_X[HmmModel::p7P_B][i] = 0.0f;
  }
  // codon[0]: no error inside the codon; codon[1]: seq[i - 1] is an insertion; codon[2]: seq[i - 2] is an insertion.
  // codon[3]: a deletion after seq[i]; codon[4]: a deletion after seq[i - 1]; codon[5]: a deletion after seq[i - 2]
  char codon[6][4];
  for (int i = 0; i < 6; i ++) {
    codon[i][3] = '\0';
  }
  //amino_acid_index is used to store the translated amino acid indices in different cases.
  // index 0 means there is no error in the codon; index 1 means seq[i - 1] is an insertion; 
  // index 2 means seq[i - 2] is an insertion. we cannot get deletion case for now.
  int amino_acid_index[3]; 
  float max_end_score = HmmModel::kInf_;  // currently best end score.

  // Begin the main DP.
  for (int i = 3; i < L + 1; i++) {
    // case 1-4. no insertion.
    codon[0][0] = seq[i - 2];
    codon[0][1] = seq[i - 1];
    codon[0][2] = seq[i];
    amino_acid_index[0] = Codon2AminoAcidIndex(codon[0]);
    // case 5-8. seq[i - 1] is an insertion.
    codon[1][0] = seq[i - 3];
    codon[1][1] = seq[i - 2];
    codon[1][2] = seq[i];        
    amino_acid_index[1] = Codon2AminoAcidIndex(codon[1]);
    // case 9-12. seq[i - 2] is an insertion.
    codon[2][0] = seq[i - 3];
    codon[2][1] = seq[i - 1];
    codon[2][2] = seq[i];        
    amino_acid_index[2] = Codon2AminoAcidIndex(codon[2]);

    // For deleted base cases, we cannot calculate the index for now 
    // because it is dependent on states.
    // A deletion after seq[i].
    // case 13-16. there is a deletion after seq[i].
    codon[3][0] = seq[i - 1];
    codon[3][1] = seq[i];
    codon[3][2] = '*';    

    // A deletion after seq[i - 1].
    // case 17-20. there is a deletion after seq[i - 1].
    codon[4][0] = seq[i - 1];
    codon[4][1] = '*';
    codon[4][2] = seq[i];

    // A deletion after seq[i - 2].
    // case 20-23. there is a deletion after seq[i - 2].
    codon[5][0] = '*';
    codon[5][1] = seq[i - 1];
    codon[5][2] = seq[i];

    float max, temp_score;  // current best score ending with i,j.
    for (int j = 1; j < M + 1; j++) {
      // Match state
      // The first block contains cases without sequencing errors in the current codon.

      // Case 0: BEGIN state.  
      max = mat_X[HmmModel::p7P_B][i - 1] + hmm.transition_[j - 1][HmmModel::B2M] 
      + hmm.match_emission_[j][amino_acid_index[0]]
      - background[amino_acid_index[0]];  
      flag[HmmModel::M_STATE][j][i] = 0;
      // Case 1: previous match state without error.	
      temp_score = mat_M[j - 1][i - 3] + hmm.transition_[j - 1][HmmModel::M2M] 
      + hmm.match_emission_[j][amino_acid_index[0]]
      - background[amino_acid_index[0]]; 
      if (temp_score > max) {
        max = temp_score;
        flag[HmmModel::M_STATE][j][i] = 1; 
      }
      // Case 2: previous insertion state without error.	
      temp_score = mat_I[j - 1][i - 3] + hmm.transition_[j - 1][HmmModel::I2M] 
      + hmm.match_emission_[j][amino_acid_index[0]]
      - background[amino_acid_index[0]];
      if (temp_score > max) {
        max = temp_score;
        flag[HmmModel::M_STATE][j][i] = 2; 
      }
      // Case 3: previous deletion state without error	
      temp_score = mat_D[j - 1][i - 3] + hmm.transition_[j - 1][HmmModel::D2M] 
      + hmm.match_emission_[j][amino_acid_index[0]]
      - background[amino_acid_index[0]];
      if (temp_score > max) {
        max = temp_score;
        flag[HmmModel::M_STATE][j][i] = 3; 
      }
      // Case 4: previous G state without error
      temp_score = mat_G[j - 1][i - 3] + hmm.transition_[j - 1][HmmModel::G2M] 
      + hmm.match_emission_[j][amino_acid_index[0]]
      - background[amino_acid_index[0]];
      if (temp_score > max) {
        max = temp_score;
        flag[HmmModel::M_STATE][j][i] = 4; 
      }
      if (i > 3) {
        // The second block contains cases where seq[i - 1] is an inserted base in the current codon
        // Case 5: previous match state with seq[i - 1] being an inserted base.
        temp_score = mat_M[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::M2M] 
        + hmm.match_emission_[j][amino_acid_index[1]]
        + seq_error_score[0][i - 1] 
        - background[amino_acid_index[1]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 5;
        }
        // Case 6: previous insertion state with seq[i - 1] being an inserted base.
        temp_score = mat_I[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::I2M] 
        + hmm.match_emission_[j][amino_acid_index[1]]
        + seq_error_score[0][i - 1] 
        - background[amino_acid_index[1]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 6; 
        }
        // Case 7: previous deletion state with seq[i - 1] being an inserted base.
        temp_score = mat_D[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::D2M] 
        + hmm.match_emission_[j][amino_acid_index[1]]
        + seq_error_score[0][i - 1] 
        - background[amino_acid_index[1]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 7;
        }
        // Case 8: previous G state with seq[i - 1] being an inserted base.
        temp_score = mat_G[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::G2M] 
        + hmm.match_emission_[j][amino_acid_index[1]]
        + seq_error_score[0][i - 1] 
        - background[amino_acid_index[1]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 8; 
        }

        // The third block contains cases where seq[i - 2] is an inserted base in the current codon
        // Case 9: previous match state with seq[i - 2] being an inserted base.
        temp_score = mat_M[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::M2M] 
        + hmm.match_emission_[j][amino_acid_index[2]]
        + seq_error_score[0][i - 2] 
        - background[amino_acid_index[2]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 9; 
        }
        // Case 10: previous insertion state with seq[i - 2] being an inserted base.
        temp_score = mat_I[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::I2M] 
        + hmm.match_emission_[j][amino_acid_index[2]]
        + seq_error_score[0][i - 2] 
        - background[amino_acid_index[2]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 10; 
        }
        // Case 11: previous deletion state with seq[i - 2] being an inserted base.
        temp_score = mat_D[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::D2M] 
        + hmm.match_emission_[j][amino_acid_index[2]]
        + seq_error_score[0][i - 2] 
        - background[amino_acid_index[2]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 11; 
        }
        // Case 12: previous G state with seq[i - 2] being an inserted base.
        temp_score = mat_G[j - 1][i - 4] + hmm.transition_[j - 1][HmmModel::G2M] 
        + hmm.match_emission_[j][amino_acid_index[2]]
        + seq_error_score[0][i - 2] 
        - background[amino_acid_index[2]];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 12; 
        }

        char hidden_base = '\0';
        float max_missing_base_score;
        // The forth block contains cases where there is a deleted base after seq[i] in the current codon.
        // Calculate the maximum score by guessing missing base.
        // position begins with 0 
        max_missing_base_score = 
          HiddenCodon2AminoAcidIndex(codon[3], 2, background, hmm.match_emission_[j], hidden_base); 
        deleted_base[0][j][i] = hidden_base; // Keep a record of the deleted base.
        // Case 13: previous match state with a deleted base after seq[i].
        temp_score = mat_M[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::M2M] 
        + max_missing_base_score
          + seq_error_score[1][i];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 13; 
        }
        // Case 14: previous insertion state with a deleted base after seq[i].
        temp_score = mat_I[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::I2M]
        + max_missing_base_score
          + seq_error_score[1][i];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 14; 
        }
        // Case 15: previous deletion state with a deleted base after seq[i]. 
        temp_score = mat_D[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::D2M] 
        + max_missing_base_score
          + seq_error_score[1][i];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 15; 
        }
        // Case 16: previous G state with a deleted base after seq[i].
        temp_score = mat_G[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::G2M]
        + max_missing_base_score
          + seq_error_score[1][i];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 16; 
        }

        // The fifth block contains cases where there is a deleted base after seq[i - 1] in the current codon.
        max_missing_base_score = 
          HiddenCodon2AminoAcidIndex(codon[4], 1, background, hmm.match_emission_[j], hidden_base); // position begins with 0
        deleted_base[1][j][i] = hidden_base; // Keep a record of the deleted base. 
        // Case 17: previous match state with a deleted base after seq[i - 1].
        temp_score = mat_M[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::M2M] 
        + max_missing_base_score
          + seq_error_score[1][i - 1];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 17; 
        }
        // Case 18: previous match state with a deleted base after seq[i - 1].
        temp_score = mat_I[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::I2M] 
        + max_missing_base_score
          + seq_error_score[1][i - 1];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 18; 
        }
        temp_score = mat_D[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::D2M]
        + max_missing_base_score
          + seq_error_score[1][i - 1];
        // Case 19: previous deletion state with a deleted base after seq[i -1]. 
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 19; 
        }
         // Case 20: previous G state with a deleted base after seq[i -1]. 
        temp_score = mat_G[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::G2M] 
        + max_missing_base_score
          + seq_error_score[1][i - 1];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 20;
        }

        // The sixth block contains cases where there is a deleted base after seq[i - 2] in the current codon.
        // position begins with 0. 
        max_missing_base_score = 
          HiddenCodon2AminoAcidIndex(codon[5], 0, background, hmm.match_emission_[j], hidden_base); 
        deleted_base[2][j][i] = hidden_base; // Keep a record of the deleted base 
        // Case 21: previous match state with a deleted base after seq[i - 2].
        temp_score = mat_M[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::M2M] 
        + max_missing_base_score
          + seq_error_score[1][i - 2];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 21; 
        }
        // Case 22: previous insertion state with a deleted base after seq[i-2].
        temp_score = mat_I[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::I2M] 
        + max_missing_base_score
          + seq_error_score[1][i - 2];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 22; 
        }
        // Case 23: previous deletion state with a deleted base after seq[i-2]
        temp_score = mat_D[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::D2M]
        + max_missing_base_score
          + seq_error_score[1][i - 2];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 23;  
        }
        // Case 24: previous G state with a deleted base after seq[i-2].
        temp_score = mat_G[j - 1][i - 2] + hmm.transition_[j - 1][HmmModel::G2M]
        + max_missing_base_score
          + seq_error_score[1][i - 2];
        if (temp_score > max) {
          max = temp_score;
          flag[HmmModel::M_STATE][j][i] = 24;  
        }
      } // if (i>3).
      mat_M[j][i] = max;
      //Insertion state
      max = mat_M[j][i - 3] + hmm.transition_[j][HmmModel::M2I] 
      + hmm.insert_emission_[j][amino_acid_index[0]] 
      - background[amino_acid_index[0]];
      flag[HmmModel::I_STATE][j][i] = 0; // Case 0: from previous match state.
      temp_score = mat_I[j][i - 3] + hmm.transition_[j][HmmModel::I2I] 
      + hmm.insert_emission_[j][amino_acid_index[0]] 
      - background[amino_acid_index[0]];
      if (temp_score > max) {
        max = temp_score;
        flag[HmmModel::I_STATE][j][i] = 1; // Case 1: a insertion state loop.
      }
      mat_I[j][i] = max;
      // Deletion state
      max = mat_M[j - 1][i] + hmm.transition_[j - 1][HmmModel::M2D];
      flag[HmmModel::D_STATE][j][i] = 0; // Case 0: from previous match state.
      temp_score = mat_D[j - 1][i] + hmm.transition_[j - 1][HmmModel::D2D];
      if (temp_score > max) {
        max = temp_score;
        flag[HmmModel::D_STATE][j][i] = 1; // Case 1: from previous deletion state.
      }
      mat_D[j][i] = max;

      // G state
      max = mat_M[j][i - 1] + hmm.transition_[j][HmmModel::M2G];
      flag[HmmModel::G_STATE][j][i] = 0; // Case 0: from previous match state.
      temp_score = mat_G[j][i - 1] + hmm.transition_[j][HmmModel::G2G];
      if (temp_score > max) {
        max = temp_score;
        flag[HmmModel::G_STATE][j][i] = 1; // Case 1: a G state loop.
      }
      mat_G[j][i] = max;

      // E state. Here do not consider the last state to be a deletion state.
      temp_score = mat_M[j][i] + hmm.transition_[j][HmmModel::M2E];
      if (temp_score > mat_X[HmmModel::p7P_E][i]) {
        E_state_flag[i][0] = j;
        mat_X[HmmModel::p7P_E][i] = temp_score;
      }
    } // for (j).
    if (mat_X[HmmModel::p7P_E][i] > max_end_score) {
      max_end_score = mat_X[HmmModel::p7P_E][i];
      E_state_flag[i][1] = 1;
    } else {
      E_state_flag[i][1] = 0;
    }
  } // for (i).
  float score = max_end_score / log(2.0f);
  string out_seq;
  int error_num = 0;
  vector< vector<int> > seq_error_pos(2);
  vector< vector<int> > state_error_pos(2);
  TraceBack(flag, E_state_flag, deleted_base, L, M, seq, out_seq, error_num, seq_error_pos, state_error_pos);
  dna.set_corrected_seq(StripOutputSeq(out_seq));
  out_file << ">" << dna.seq_name() << " hmm_name=" 
    << hmm.name() << " score=" << score 
    << " error_num=" << error_num;
  out_file << " seq_insertion_pos=";
  for (int i = 0; i < seq_error_pos[0].size(); ++i) {
    out_file << seq_error_pos[0][i] << ",";
  }
  out_file << " seq_deletion_pos=";
  for (int i = 0; i < seq_error_pos[1].size(); ++i) {
    out_file << seq_error_pos[1][i] << ",";
  }
  out_file << " state_insertion_pos=";
  for (int i = 0;i < state_error_pos[0].size(); ++i) {
    out_file << state_error_pos[0][i] << ",";
  }
  out_file << " state_deletion_pos=";
  for (int i = 0;i < state_error_pos[1].size(); ++i) {
    out_file << state_error_pos[1][i] << ",";
  }
  out_file << endl;
  out_file << dna.corrected_seq() << endl;
}

void TraceBack(const vector<vector<vector<int> > >& flag, 
               const vector<vector<int> >& E_state_flag, 
               const vector<vector<vector<char> > >& deleted_base,
               int L, int M, const char* seq, string& out_seq, 
               int& error_num, vector<vector<int> >& seq_error_pos, 
               vector<vector<int> >& state_error_pos) {
  int begin_base, end_base, begin_state, end_state; 
  end_base = L;
  for (int i = L; i > 0; i--) {
    if(E_state_flag[i][1] == 1)
    {
      end_base = i;
      break;
    }
  }
  end_state = E_state_flag[end_base][0];
  out_seq.clear();  // The output sequence. Begin base not known yet. Index begins with 1. out_seq[0] is '$' symbo
  int choice = 0; // choice[0]: match state; choice[1]: insertion state; choice[2]: deletion state; choice[3]: G state; choice[4]:BEGIN
  stack<string> path; // The traceback path
  int index1 = end_base;
  int index2 = end_state;
  const int kLabelLength = 10; // The length of each label of the DNA base
  char label[kLabelLength];
  char buff[kLabelLength];
  char del_base;
  char temp[2];
  // insertion position in the original sequence.
  vector<int> seq_insertion_pos;
  // deletion position after some base in the original sequence. 
  // 2 means there is a deletion after the second base.
  vector<int> seq_deletion_pos;  // insertion and deletion position
  vector<int> state_insertion_pos;
  vector<int> state_deletion_pos;
  while (index1 > 1 && index2 > 0) {
    if(choice == 0) {
      switch(flag[0][index2][index1]) {
      case 0: // B->M.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(seq[index1 - 2]);
        strcpy(label,"M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push(label);
        path.push(label);
        index1 -= 2;
        begin_base = index1;
        begin_state = index2;
        sprintf(buff,"%d", begin_base);
        strcpy(label, "B");
        strcat(label, buff); 
        path.push(label);
        choice = 4;
        continue;
      case 1: // M->M.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(seq[index1 - 2]);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push(label);
        path.push(label);
        index1 -= 3;
        index2--;
        choice = 0;
        continue;
      case 2: //I->M.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(seq[index1 - 2]);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push(label);
        path.push(label);
        index1 -= 3;
        index2--;
        choice = 1;
        continue;
      case 3: //D->M.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(seq[index1 - 2]);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push(label);
        path.push(label);
        index1 -= 3;
        index2--;
        choice = 2;
        continue;
      case 4: //G->M.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(seq[index1 - 2]);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label,buff);
        path.push(label);
        path.push(label);
        path.push(label);
        index1 -= 3;
        index2--;
        choice = 3;
        continue;
      case 5:  // M->M, index1 - 1 is an insertion. 
        out_seq.push_back(seq[index1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 2]);
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 1);
        state_insertion_pos.push_back(index2);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push("i*");
        path.push(label);
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 0;
        continue;
      case 6:  // I->M, index - 1 is an insertion.
        out_seq.push_back(seq[index1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 2]);
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 1);
        state_insertion_pos.push_back(index2);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push("i*");
        path.push(label);
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 1;
        continue;
      case 7:  // D->M, index - 1 is an insertion.
        out_seq.push_back(seq[index1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 2]);
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 1);
        state_insertion_pos.push_back(index2);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push("i*");
        path.push(label);
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 2;
        continue;
      case 8:  // G->M, index - 1 is an insertion.
        out_seq.push_back(seq[index1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 2]);
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 1);
        state_insertion_pos.push_back(index2);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push("i*");
        path.push(label);
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 3;
        continue;
      case 9:  // M->M, index1 - 2 is an insertion.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 2);
        state_insertion_pos.push_back(index2);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label,buff);
        path.push(label);
        path.push(label);
        path.push("i*");
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 0;
        continue;
      case 10: // I->M, index1 - 2 is an insertion.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 2);
        state_insertion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff, "%d", index2);
        strcat(label, buff);
        path.push(label);
        path.push(label);
        path.push("i*");
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 1;
        continue;
      case 11:  // D->M, index1 - 2 is an insertion.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 2);
        state_insertion_pos.push_back(index2);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label,buff);
        path.push(label);
        path.push(label);
        path.push("i*");
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 2;
        continue;
      case 12:  // G->M, index1 - 2 is an insertion.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back('i');
        out_seq.push_back(seq[index1 - 3]);
        error_num++;
        seq_insertion_pos.push_back(index1 - 2);
        state_insertion_pos.push_back(index2);
        strcpy(label, "M");
        sprintf(buff, "%d", index2);
        strcat(label,buff);
        path.push(label);
        path.push(label);
        path.push("i*");
        path.push(label);
        index1 -= 4;
        index2--;
        choice = 3;
        continue;
      case 13:  // M->M, there is a deletion after index1.
        del_base = tolower(deleted_base[0][index2][index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(temp);
        path.push(label);
        path.push(label);
        choice = 0;
        index1 -= 2;
        index2--;
        continue;
      case 14:  // I->M, there is a deletion after index1.
        del_base = tolower(deleted_base[0][index2][index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(temp);
        path.push(label);
        path.push(label);
        choice = 1;
        index1 -= 2;
        index2--;
        continue;
      case 15:  // D->M, there is a deletion after index1.
        del_base = tolower(deleted_base[0][index2][index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(temp);
        path.push(label);
        path.push(label);
        choice = 2;
        index1 -= 2;
        index2--;
        continue;
      case 16:  // G->M, there is a deletion after index1.
        del_base = tolower(deleted_base[0][index2][index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(temp);
        path.push(label);
        path.push(label);
        choice = 3;
        index1 -= 2;
        index2--;
        continue;
      case 17:  // M->M, there is a deletion after index1 - 1.
        del_base = tolower(deleted_base[1][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1 - 1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(temp);
        path.push(label);
        choice = 0;
        index1 -= 2;
        index2--;
        continue;    
      case 18:  // I->M, there is a deletion after index1 - 1
        del_base = tolower(deleted_base[1][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1 - 1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(temp);
        path.push(label);
        choice = 1;
        index1 -= 2;
        index2--;
        continue;
      case 19:  // D->M, there is a deletion after index1 - 1.
        del_base = tolower(deleted_base[1][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1 - 1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(temp);
        path.push(label);
        choice = 2;
        index1 -= 2;
        index2--;
        continue;
      case 20:  // G->M, there is a deletion after index1 - 1.
        del_base = tolower(deleted_base[1][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(del_base);
        out_seq.push_back(seq[index1 - 1]);
        error_num++;
        seq_deletion_pos.push_back(index1 - 1);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(temp);
        path.push(label);
        choice = 3;
        index1 -= 2;
        index2--;
        continue;
      case 21:  // M->M, there is a deletion after index - 2.
        del_base = tolower(deleted_base[2][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(del_base);
        error_num++;
        seq_deletion_pos.push_back(index1 - 2);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(label);
        path.push(temp);
        choice = 0;
        index1 -= 2;
        index2--;
        continue;
      case 22:  // I->M, there is a deletion after index1 - 2.
        del_base = tolower(deleted_base[2][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(del_base);
        error_num++;
        seq_deletion_pos.push_back(index1 - 2);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(label);
        path.push(temp);
        choice = 1;
        index1 -= 2;
        index2--;
        continue;
      case 23:  // D->M, there is a deletion after index1 - 2.
        del_base = tolower(deleted_base[2][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(del_base);
        error_num++;
        seq_deletion_pos.push_back(index1 - 2);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(label);
        path.push(temp);
        choice = 2;
        index1 -= 2;
        index2--;
        continue;
      case 24:  // G->M, there is a deletion after index1 - 2.
        del_base = tolower(deleted_base[2][index2][index1]);
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(del_base);
        error_num++;
        seq_deletion_pos.push_back(index1 - 2);
        state_deletion_pos.push_back(index2);
        strcpy(label,"M");
        sprintf(buff,"%d",index2);
        strcat(label,buff);
        temp[0] = del_base;
        temp[1] = '\0';
        path.push(label);
        path.push(label);
        path.push(temp);
        choice = 3;
        index1 -= 2;
        index2--;
        continue;
      }
    }
    if (choice == 1)
    {    
      switch (flag[1][index2][index1])
      {
      case 0:  // M->I.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(seq[index1 - 2]);
        strcpy(label, "I");
        sprintf(buff, "%d", index2);
        strcat(label, buff);    
        path.push(label);
        path.push(label);    
        path.push(label);
        index1 -= 3;
        choice = 0;
        continue;
      case 1:  // I->I.
        out_seq.push_back(seq[index1]);
        out_seq.push_back(seq[index1 - 1]);
        out_seq.push_back(seq[index1 - 2]);
        strcpy(label, "I");
        sprintf(buff, "%d", index2);
        strcat(label, buff);    
        path.push(label);
        path.push(label);
        path.push(label);
        index1 -= 3;
        choice = 1;
        continue;
      }
    }
    if (choice == 2)
    {
      switch (flag[2][index2][index1]) {
      case 0:  // M->D.
        strcpy(label, "D");
        sprintf(buff, "%d", index2);
        strcat(label, buff);    
        path.push(label);
        index2--;
        choice = 0;
        continue;
      case 1:  // D->D.
        strcpy(label, "D");
        sprintf(buff, "%d", index2);
        strcat(label, buff);    
        path.push(label);
        index2--;
        choice = 2;
        continue;
      }
    }
    if (choice == 3)
    {
      switch (flag[3][index2][index1])
      {
      case 0:  // M->G.
        out_seq.push_back('i');
        error_num++;
        seq_insertion_pos.push_back(index1);
        state_insertion_pos.push_back(index2);
        strcpy(label, "G");
        sprintf(buff, "%d", index2);
        strcat(label,buff);
        path.push(label);
        index1--;
        choice = 0;
        continue;
      case 1:  // G->G.
        out_seq.push_back('i');
        error_num++;
        seq_insertion_pos.push_back(index1);
        state_insertion_pos.push_back(index2);
        strcpy(label, "G");
        sprintf(buff,"%d", index2);
        strcat(label, buff);
        path.push(label);
        index1--;
        choice = 3;
        continue;
      }
    }
    if (choice == 4)
    {
      break;
    }        
  }
  seq_error_pos[0] = seq_insertion_pos;
  seq_error_pos[1] = seq_deletion_pos;
  state_error_pos[0] = state_insertion_pos;
  state_error_pos[1] = state_deletion_pos;
  reverse(out_seq.begin(), out_seq.end());  // output the original sequence order.
}
// This function is called in main function due to L
int p7_ReconfigLength(HmmModel& hmm, int seq_length)
{
  float ploop, pmove;  
  /* Configure N,J,C transitions so they bear L/(2+nj) of the total
  * unannotated sequence length L. 
  */
  pmove = (2.0f + hmm.nj()) / ((float) seq_length + 2.0f + hmm.nj()); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;
  hmm.special_transition_[HmmModel::p7P_N][HmmModel::p7P_LOOP] =  hmm.special_transition_[HmmModel::p7P_C][HmmModel::p7P_LOOP]
  = hmm.special_transition_[HmmModel::p7P_J][HmmModel::p7P_LOOP] = log(ploop);
  hmm.special_transition_[HmmModel::p7P_N][HmmModel::p7P_MOVE] =  hmm.special_transition_[HmmModel::p7P_C][HmmModel::p7P_MOVE] 
  = hmm.special_transition_[HmmModel::p7P_J][HmmModel::p7P_MOVE] = log(pmove);
  return 0;
}

float HiddenCodon2AminoAcidIndex(char* codon, int missing_position, vector<float> background, vector<float> match_emission, char& missing_base) {

  float max = HmmModel::kInf_; // max_missing_base_score
  float temp; 
  int guess_amino_acid_index; // Guess a base
  int max_amino_acid_index; // The base which maximizes the emission score minus background score
  // Begin to guess the missing base
  codon[missing_position] = 'A';
  guess_amino_acid_index = Codon2AminoAcidIndex(codon);
  temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
  if (temp > max) {
    max = temp;
    max_amino_acid_index = guess_amino_acid_index;
    missing_base = 'A';
  }
  codon[missing_position] = 'C';
  guess_amino_acid_index = Codon2AminoAcidIndex(codon);
  temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
  if (temp > max) {
    max = temp;
    max_amino_acid_index = guess_amino_acid_index;
    missing_base = 'C';
  }
  codon[missing_position] = 'G';
  guess_amino_acid_index = Codon2AminoAcidIndex(codon);
  temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
  if (temp > max) {
    max = temp;
    max_amino_acid_index = guess_amino_acid_index;
    missing_base = 'G';
  }
  codon[missing_position] = 'T';
  guess_amino_acid_index = Codon2AminoAcidIndex(codon);
  temp = match_emission[guess_amino_acid_index] - background[guess_amino_acid_index];
  if (temp > max) {
    max = temp;
    max_amino_acid_index = guess_amino_acid_index;
    missing_base = 'T';
  }
  return max;
}

// remove 'i' and change lowercase to uppercase.
string StripOutputSeq(const string& seq) {
  string out_seq;
  for (int i = 0; i < seq.length(); ++i) {
    char ch = seq[i];
    if (ch  == 'i') continue;
    else {
      out_seq.push_back(toupper(ch));
    }
  }
  return out_seq;
}

// output the alignment path.
void OutputPath(stack<string>& path) {
  string state;
  while (!path.empty()) {
    state = path.top();
    cout << state << " ";
    path.pop();
  }
  cout << endl;
}
