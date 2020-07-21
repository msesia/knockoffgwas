#ifndef _DATASET_CPP
#define _DATASET_CPP

#include "dataset.h"

#define DELIMITER "--------------------------------------------------------------------------------"

Dataset::Dataset() {  
}

Dataset::Dataset(const vector<Metadata>& _metadata, int _num_threads, int thinning) {
  num_chrs_ = _metadata.size();
  num_threads = _num_threads;

  // Load chromosomes one by one
  for(int chr=0; chr<num_chrs_; chr++) { 
    // Reduce number of SNPs (thinning)
    metadata.push_back(_metadata[chr].thin(thinning));

    cout << "Chromosome " << metadata[chr].get_chr_id() << " will be loaded from:" << endl;
    cout << "  haplotype file            : " << metadata[chr].df.data + "." + metadata[chr].df.format << endl;
    cout << "  sample file               : " << metadata[chr].df.sample << endl;
    cout << "  legend file               : " << metadata[chr].df.legend << endl;
    if(thinning!=1) {
      cout << "  thinning factor           : " << thinning << endl;
    }
    cout << "  sample filter file        : " << metadata[chr].df.keep << endl;
    cout << "  variant filter file       : " << metadata[chr].df.extract << endl;
    cout << "  number of SNPs            : " << metadata[chr].num_snps() << endl;
    cout << "  number of haplotypes      : " << metadata[chr].num_haps() << endl;
    cout << endl;

    // Insert new chromosome in list
    chromosomes.push_back(Haplotypes(metadata[chr]));

    // Load data
    chromosomes[chr].load_data(true, num_threads);
  }

}

void Dataset::load_worker(int chr_min, int chr_max, int wid, vector<int> &progress) {  
  int n_steps = num_chrs_;

  // Initialize progress bar (worker 0)
  if(wid==0) {
    cout << "|" << flush;
    for(int s=0; s<n_steps; s++) cout << ".";
    cout << "|" << endl;
    cout << "|" << flush;
  }

  //cout << "Worker " << wid << ": chromosomes " << chr_min << "-" << chr_max-1 << endl;
  // Load data for each chromosome
  for(int chr=chr_min; chr<chr_max; chr++) {
    chromosomes[chr].load_data(true, 1); // FIXME: this is now obsolete
    
    // Keep track of progress (all workers)
    progress[wid]++;
    cout << "=" << flush;
  }
  // Finalize progress bar (all workers)
  if( std::accumulate(progress.begin(),progress.end(),0) == n_steps) {
    cout << "|" << endl;
  }

}

int Dataset::count_lines(string filename) const {
  int num_lines = 0;
  string line;
  ifstream file(filename);
  while (std::getline(file, line)) num_lines++;
  file.close();
  return(num_lines);
}

int Dataset::num_chrs() const {
  return(chromosomes.size());
}

#endif
