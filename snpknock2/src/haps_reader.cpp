#ifndef _HAPS_READER_CPP
#define _HAPS_READER_CPP

#include "haps_reader.h"

HapsReader::HapsReader(const string& _filename, const string& sample_filename, const string& legend_filename) {
  filename = _filename;

  string buffer;
  vector<string> tokens;
  int line = -2;

  // Read sample IDs from ".sample" file
  ifile fd_sample(sample_filename);
  while (getline(fd_sample, buffer, '\n')) {
    if(line>=0) {
      sutils::tokenize(buffer, tokens);
      sample.push_back(tokens);
    }
    line++;
  }
  fd_sample.close();

  // Read variant legend ".legend" file
  line = -1;
  ifile fd_legend(legend_filename);
  while (getline(fd_legend, buffer, '\n')) {
    if(line>=0) {
      sutils::tokenize(buffer, tokens);
      legend.push_back(tokens);
    }
    line++;
  }
  fd_legend.close();
}

void HapsReader::read(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H,
                      bool verbose) {
  
  // cout << "[DEBUG] in HapsReader::read" << endl;
  // cout << "Requested SNP ids: " << endl;
  // for(auto id : snp_ids) cout << " " << id;
  // cout << endl;
  // cout << "Variant legend: " << endl;
  // for(auto id : legend.ID) cout << " " << id;
  // cout << endl;

  // Initialize data container
  int num_samples = sample_ids.size();
  int num_snps = snp_ids.size();
  int num_haps = 2 * num_samples;

  // // Find HAPS row indices of requested samples
  // vector<int> rows_idx(num_snps,-1);
  // for(int j=0; j<num_snps; j++) {
  //   auto it = std::find(legend.ID.begin(), legend.ID.end(), snp_ids[j]);
  //   if(it != legend.ID.end()) {
  //     rows_idx[j] = std::distance(legend.ID.begin(), it);
  //   } else {
  //     cerr << "Error parsing HAPS file. Could not find variant " << snp_ids[j] << endl;
  //     exit(1);
  //   }
  // }
  // std::sort(rows_idx.begin(), rows_idx.end());

  // cout << "[DEBUG] will read the following rows: ";
  // for(int j=0; j<num_snps; j++) {
  //   cout << rows_idx[j] << " ";
  // }

  // Find HAPS col indices of requested samples
  vector<int> hap_idx(num_haps, -1);
  for(int i=0; i<num_samples; i++) {
    auto it = std::find(sample.ID.begin(), sample.ID.end(), sample_ids[i]);
    if(it != sample.ID.end()) {
      int i_haps = std::distance(sample.ID.begin(), it);
      hap_idx[2*i] = 2*i_haps;
      hap_idx[2*i+1] = 2*i_haps+1;
    } else {
      cerr << "Error parsing HAPS file. Could not find sample " << sample_ids[i] << endl;
      exit(1);
    }
  }

  //cout << "Found col indices" << endl;

  // TODO: show progress

  // Load requested rows and columns
  string buffer;
  vector <string> tokens;
  vector < vector<bool> > haps(num_haps, vector<bool>(num_snps, false));
  int row = 0;
  ifile fd_hap(filename);
  while (getline(fd_hap, buffer, '\n')) {
    // Skip this row if the variant was not requested
    string rsid = legend.ID[row];
    auto it = std::find(snp_ids.begin(), snp_ids.end(), rsid);
    if(it == snp_ids.end()) {
      row++;
      continue;
    }
    // Find output index of variant, if requested
    int j = std::distance(snp_ids.begin(), it);
    // Process row
    sutils::tokenize(buffer, tokens);
    for(int i=0; i<num_haps; i++) {
      int i_haps = hap_idx[i];      
      haps[i][j] = (tokens[i_haps] == "1");
    }
    row++;
  }
  fd_hap.close();
  // Store transposed haplotypes
  H.clear();
  for (int r = 0 ; r < haps.size() ; r ++){
    H.push_back(chaplotype(haps[r]));
  }

  // TODO: make sure all requested data have been read

  // cout << "[DEBUG] exiting HapsReader::read" << endl;
  // for(int i=0; i<20; i++) {
  //   for(int j=0; j<50; j++) {
  //     cout << H[i][j];
  //   }
  //   cout << endl;
  // }

}

#endif
