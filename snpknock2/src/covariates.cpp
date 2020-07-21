#ifndef COVARIATES_CPP
#define COVARIATES_CPP

#define DEBUG 0

#include "covariates.h"

Covariates::Covariates() {
  K = 0;
  K_max = 10;
}

Covariates::Covariates(const string & pc_filename, const Metadata & metadata) : Covariates() {
  vector< vector<double> > Z_raw;
  Sample sample_raw;
  load_data(pc_filename, K_max, sample_raw, Z_raw);
  K = Z_raw[0].size();
  filter(metadata, sample_raw, Z_raw);
  mask_by_family(metadata);
}

Covariates::Covariates(const Covariates& obj) {
  K = obj.K;
  Z = obj.Z;
  sample = obj.sample;
}

Covariates::~Covariates() {
}

void Covariates::clear() {
  Z.clear();
}

void Covariates::mask_by_family(const Metadata & metadata) {
  if(DEBUG) cout << "[DEBUG] in mask_by_family()" << endl;

  vector< vector<double> > Z_masked = Z;

  for(int i=0; i<Z.size(); i++) {
    if(DEBUG) cout << "[DEBUG] i = " << i << ": ";
    const int fam_size = metadata.related_samples[i].size();
    vector<double> z_tmp = Z[i];
    for(int k=0; k<K; k++) {
      Z[i][k] /= (double) (fam_size+1);
    }
    for(int s=0; s<fam_size; s++) {
      int j = metadata.related_samples[i][s];
      if(DEBUG) cout << j << " ";
      for(int k=0; k<K; k++) {
        z_tmp[k] += Z[j][k] / (double) (fam_size+1);
      }
    }
    Z_masked[i] = z_tmp;
    if(DEBUG) cout << endl;
  }

  Z = Z_masked;
}

void Covariates::filter(const Metadata & metadata, Sample & sample_raw, vector< vector<double> > & Z_raw) {
  //cout << "[DEBUG] Filtering Z observations" << endl;

  vector<int> idx_match;
  right_join(metadata.sample_filter.ID, sample_raw.ID, idx_match);

  sample.clear();
  Z.clear();
  for(int i=0; i<idx_match.size(); i++) {
    int j = idx_match[i];    
    sample.push_back(sample_raw.row(j));
    Z.push_back(Z_raw[j]);
  }

  cout << "Matched principal components to " << idx_match.size() << " individuals." << endl;

  // Make sure list of filtered samples match metadata
  if(sample.size() != metadata.sample_filter.ID.size()) {
      cout << "Error in Z file: sample IDs do not match data IDs" << endl;
      exit(1);
  }
  for(int i=0; i<sample.size(); i++) {
    if(sample.ID[i] != metadata.sample_filter.ID[i]) {
      cout << "Error in Z file: sample IDs do not match data IDs" << endl;
      cout << "i " << i << " : " << sample.ID[i] << " " << metadata.sample_filter.ID[i] << endl;
      exit(1);
    }
  }
}

void Covariates::print() const {
  cout << "Printing covariates:" << endl;
  for(int i=0; i<Z.size(); i++) {
    for(int k=0; k<Z[i].size(); k++) {
      cout << " " << Z[i][k];
    }
    cout << endl;
  }
  cout << endl;
}

void load_data(const string & filename, int K_max, Sample & sample, vector< vector<double> > & Z) {
  // Load genetic principal components
  cout << "Loading genetic principal components from " << filename << endl;
  int K = 0;
  Z.clear();
  vector<double> pc_means;
  vector<double> pc_counts;
  vector<string> sample_row(4, "0");
  int line_num = 0;
  string buffer;
  std::vector<string> tokens;
  ifstream file(filename);
  while (std::getline(file, buffer)) {
    boost::split(tokens, buffer, boost::is_any_of(" ,\t"));
    if(line_num == 0) {
      K = tokens.size() - 2;
      if(K>K_max) K=K_max;
      if(K < 0) {
        cout << "Error in load_data(): incorrect header format." << endl;
        exit(1);
      }
      pc_means.resize(K,0);
      pc_counts.resize(K,0);
    } else {
      if(tokens.size() < (K+2)) {
        cout << "Error in load_data(): incorrect number of columns in row " << line_num + 1 << "." << endl;
        exit(1);
      }
      vector<double> pc_row(K);
      for(int k=0; k<K; k++) {
        if(tokens[2+k]=="NA") {
          // Impute missing value
          pc_row[k] = pc_means[k];          
          //cout << "Found NA on line "<< line_num + 1 << ", Z " << k << ", replaced with " << pc_means[k] << endl;
        } else {
          // Use real value
          pc_row[k] = std::stod(tokens[2+k]);
          // Update column mean
          pc_means[k] = ((pc_means[k]*pc_counts[k]) + pc_row[k])/(pc_counts[k]+1);
          pc_counts[k]++;
        }
      }
      Z.push_back(pc_row);
      sample_row[0] = tokens[0];
      sample_row[1] = tokens[1];
      sample.push_back(sample_row);
    }
    line_num++;
  }
  file.close();
  assert(Z.size()>0);
  cout << "Loaded " << K << " genetic principal components for " << Z.size() << " individuals." << endl;
  cout << endl;
}

#endif
