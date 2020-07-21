#ifndef _GENERATOR_MASTER_H
#define _GENERATOR_MASTER_H

//#include "filter_reader.h"
#include "knockoffs.h"
#include "utils.h"
#include "metadata.h"

typedef vector< vector<int> > ivector2d;
typedef vector< ivector2d > ivector3d;

class KnockoffGenerator {
public:
  // Constructors/desctructors
  KnockoffGenerator(const vector<Metadata>& _metadata, 
                    const ivector3d & references_local, const ivector2d & references_global, 
                    int _num_threads, int _debug, int _seed, const string & _logfile);

  // Methods
  void generate();
  void fit_HMM(const double hmm_rho);
  void tune_hmm();
  void load_hmm(vector<string> hmm_files);
  void init_hmm(double hmm_rho, double hmm_lambda);
  void set_partition(int r);
  int num_partitions() const;

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
  int K;
  int num_chrs;
  int num_haps;
  int num_threads;
  int debug;
  int seed;

  // Methods
  void generate_worker(int chr, const vector<int> & i_list, boost::random::taus88 & rng, int wid, 
                       vector<int> &progress, int progress_total);
  void check_sanity();
};

#endif
