#ifndef _GENERATOR_MASTER_H
#define _GENERATOR_MASTER_H

//#include "filter_reader.h"
#include "knockoffs.h"
#include "utils.h"
#include "metadata.h"

typedef vector< vector<unsigned int> > ivector2d;
typedef vector< ivector2d > ivector3d;

class KnockoffGenerator {
public:
  // Constructors/desctructors
  KnockoffGenerator(const vector<Metadata>& _metadata, 
                    const ivector3d& references_local, const ivector2d& references_global,
                    unsigned int _num_threads, unsigned int _debug, unsigned int _seed, const string& _logfile);


  // Methods
  void generate();
  void fit_HMM(const double hmm_rho);
  void tune_hmm();
  void load_hmm(vector<string> hmm_files);
  void init_hmm(double hmm_rho, double hmm_lambda);
  void set_partition(unsigned int r);
  unsigned int num_partitions() const;

  // Output
  void writeKnockoffs(const vector<string> & out_file_name) const;
  void writeGroups(const vector<string> & out_file_name) const;
  void writeWindows(const vector<string> & out_file_name) const;
  void writeAncestries(const vector<string> & out_file_name) const;
  void writeHMM(const vector<string> & out_file_name) const;
  // void writeZ(const vector<string> & out_file_name) const;

private:
  // Data
  vector<Knockoffs> chromosomes;
  vector<Metadata> metadata;
  unsigned int K;
  unsigned int num_chrs;
  unsigned int num_haps;
  unsigned int num_threads;
  unsigned int debug;
  unsigned int seed;

  // Methods
  void generate_worker(unsigned int chr, const vector<unsigned int> & i_list, boost::random::taus88 & rng, unsigned int wid, 
                       vector<unsigned int> &progress, unsigned int progress_total);
  void check_sanity();
};

#endif
