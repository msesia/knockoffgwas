#ifndef _GENERATOR_MASTER_CPP
#define _GENERATOR_MASTER_CPP

#include "generator_master.h"

#define DELIMITER "--------------------------------------------------------------------------------"

KnockoffGenerator::KnockoffGenerator(const vector<Metadata>& _metadata, 
                                     const ivector3d& references_local, const ivector2d& references_global,
                                     int _num_threads, int _debug, int _seed, const string& _logfile) {

  metadata = _metadata;
  num_chrs = metadata.size();
  num_threads = _num_threads;
  K = references_local[0][0].size();
  debug = _debug;
  seed = _seed;

  // Number of haplotypes for samples that passed filters
  num_haps = metadata[0].num_samples() * 2;

  // Initialize each chromosome
  for(int chr=0; chr<num_chrs; chr++) {
    // Make sure the number of haplotypes is consistent
    assert(num_haps == metadata[chr].num_samples() * 2);

    // Make sure the number of windows is consistent
    const int num_windows = references_local[chr].size();
    assert(num_windows == metadata[chr].windows.num_windows);

    cout << "Chromosome " << metadata[chr].get_chr_id() << " will be loaded from:" << endl;
    cout << "  haplotype file            : " << metadata[chr].df.data << endl;
    cout << "  haplotype file format     : " << metadata[chr].df.format << endl;
    cout << "  sample file               : " << metadata[chr].df.sample << endl;
    cout << "  legend file               : " << metadata[chr].df.legend << endl;
    if(metadata[chr].df.map!="")
      cout << "  map file                  : " << metadata[chr].df.map << endl;
    else
      cout << "  map file                  : " << "NONE" << endl;
    cout << "  sample filter file        : " << metadata[chr].df.keep << endl;
    cout << "  variant filter file       : " << metadata[chr].df.extract << endl;
    cout << "  number of SNPs            : " << metadata[chr].num_snps() << endl;
    cout << "  number of windows         : " << metadata[chr].windows.num_windows << endl;
    cout << "  number of haplotypes      : " << metadata[chr].num_haps() << endl;
    cout << endl;

    // Add chromosome in list and load data
    cout << "Loading data for chromosome " << metadata[chr].get_chr_id() << endl;
    chromosomes.push_back(Knockoffs(metadata[chr], K, num_threads, debug, seed, _logfile));
  }

  // Make sure that the data looks reasonable
  check_sanity();

  // Assign references
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].set_references(references_local);
    chromosomes[chr].set_references_global(references_global);
  }

  cout << endl;
}

void KnockoffGenerator::set_partition(int r) {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].set_partition(r);
  }
}

void KnockoffGenerator::check_sanity() {
  // Make sure dataset is not empty
  assert(num_haps>0);
  // Make sure that at least one chromosome has been loaded
  assert(chromosomes.size() > 0);
  // Make sure that each chromosome has the same number of samples
  for(int chr=0; chr<num_chrs; chr++) {
    assert(chromosomes[chr].get_num_haps()==num_haps);
  }
  // Make sure that K is not too big
  if(K>num_haps) K = std::min(K,num_haps-1);
}

void KnockoffGenerator::tune_hmm() {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].CV(num_threads);
  }
}

void KnockoffGenerator::load_hmm(vector<string> hmm_files) {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].load_hmm(hmm_files[chr]);
  }
}

void KnockoffGenerator::init_hmm(double hmm_rho, double hmm_lambda) {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].init_hmm(K, hmm_rho, hmm_lambda);
  }
}

void KnockoffGenerator::fit_HMM(const double hmm_rho) {
  // Run EM chromosome by chromosome
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].init_hmm(K, hmm_rho);
    chromosomes[chr].EM(num_threads);
  }  
}

void KnockoffGenerator::generate() {
  // Generate knockoffs chromosome by chromosome
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].generate(num_threads);
  }
  cout << endl;
}

void KnockoffGenerator::writeKnockoffs(const vector<string> & out_file_name) const {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].write(out_file_name[chr], "BED", true);
    //chromosomes[chr].write(out_file_name[chr], "BED", false);
    //if(chromosomes[chr].get_num_snps() < 1000) {
    //chromosomes[chr].write(out_file_name[chr], "HAPS", false);
    //chromosomes[chr].write(out_file_name[chr], "HAPS", true);
    //}
  }
}

void KnockoffGenerator::writeAncestries(const vector<string> & out_file_name) const {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].writeAncestries(out_file_name[chr]);
  }
}

// void KnockoffGenerator::writeZ(const vector<string> & out_file_name) const {
//   for(int chr=0; chr<num_chrs; chr++) {
//     chromosomes[chr].writeZ(out_file_name[chr]);
//     chromosomes[chr].writeZk(out_file_name[chr]);
//   }
// }

void KnockoffGenerator::writeGroups(const vector<string> & out_file_name) const {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].writeGroups(out_file_name[chr]);
  }
}

void KnockoffGenerator::writeWindows(const vector<string> & out_file_name) const {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].writeWindows(out_file_name[chr]);
  }
}

void KnockoffGenerator::writeHMM(const vector<string> & out_file_name) const {
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes[chr].writeHMM(out_file_name[chr]);
  }
  cout << endl;
}

int KnockoffGenerator::num_partitions() const {
  return chromosomes.back().num_partitions();
}

#endif
