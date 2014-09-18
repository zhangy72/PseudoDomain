#include <string>
#include <vector>
#include <algorithm>
#include "DNA_seq.h"
#include "hmm_model.h"
#include "viterbi_model.h"

using namespace std;
int main(int argc, char* argv[]) {
  if (argc != 5) {
    printf("Usage: %s <hmm file> <fasta file> <output file> <error option>\n",
           argv[0]);
    return 1;
  }
  const char* hmm_file_name = argv[1];
  const char* fasta_file_name = argv[2];
  const char* out_file_name = argv[3];
  string error_para(argv[4]);
  // This block is for testing.
  ifstream hmm_file(hmm_file_name);
  if (!hmm_file) {
    printf("Cannot open HMM file!\n");
    return 1;
  }
  ifstream seq_file(fasta_file_name);
  if (!seq_file.is_open()) {
    printf("Cannot open fasta file!\n");
    hmm_file.close();
    return 1;
  }
  ofstream out_file(out_file_name);
  if (!out_file.is_open()) {
    printf("Cannot open output file!\n");
    hmm_file.close();
    seq_file.close();
    return 1;     
  }
  const float kNonHomoErrorRate = 0.0007f;
  const float kHomoErrorRate = 0.0044f;
  bool model_option;
  if (error_para == "0") 
    model_option = false;
  else 
    model_option = true;
  vector<HmmModel> hmms;
  while (1) {
    HmmModel hmm;
    if (ConstructHmmModel(hmm_file, hmm, kNonHomoErrorRate)) 
      hmms.push_back(hmm);
    else 
      break;
  }
  hmm_file.close();
  for (int i = 0; i < hmms.size(); ++i) {
    //Begin the Viterbi process. Each time store one sequence and process it. 
    string line;
    string seq;
    string seq_name;
    while (getline(seq_file, line)) {
      if (line.find_first_not_of(" \t") == string::npos) {
        continue;
      }
      if (line[0] == '>') {
        if (!seq_name.empty() && !seq.empty()) {
          DNASeq dna(seq_name, seq, kNonHomoErrorRate, kHomoErrorRate);
          dna.set_error_model(model_option);
          Viterbi(hmms[i], dna, out_file);
          seq.clear();
        }
        seq_name = line.substr(1); // seq_name will not include the prefix ">"
      } else {
        seq.append(line);
      }
    }
    // handle the last string
    DNASeq dna(seq_name, seq, kNonHomoErrorRate, kHomoErrorRate);
    dna.set_error_model(model_option);
    Viterbi(hmms[i], dna, out_file);
    
    // set the file stream back to the beginning.
    seq_file.clear();
    seq_file.seekg(0);
  }
  seq_file.close();
  out_file.close();
  return 0;
}

