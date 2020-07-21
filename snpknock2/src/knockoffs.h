#ifndef _KNOCKOFFS_H
#define _KNOCKOFFS_H

#include <fstream>
#include <functional>
#include <boost/thread.hpp>
#include <boost/date_time.hpp>
#include "haplotypes.h"
#include "utils.h"
#include "metadata.h"
#include "plink_interface.h"
#include "cross_validation.h"
#include "family_bp.h"
#include "family_knockoffs.h"
//#include </home/groups/candes/Software/usr/include/gperftools/profiler.h>

typedef vector< vector<int> > ivector2d;
typedef vector<ivector2d> ivector3d;

class Knockoffs : public Haplotypes {

public:
  // Constructors
  Knockoffs(const Metadata& _metadata, int _K, int nthreads, int _debug, int _seed, const string & _logfile);
  Knockoffs(const Metadata& _metadata, int _K, bool _load_data, int nthreads, int _debug, int _seed,
            const string & _logfile);
  void load(const int nthreads);

  // Get/set
  void set_reference(const int i, const ivector2d & input);
  void set_references(const ivector3d& input);
  void get_reference(const int i, ivector2d& output) const;
  void set_reference_global(const int i, const vector<int>& input);
  void set_references_global(const ivector2d& input);
  void get_reference_global(const int i, vector<int>& output) const;
  void set_ancestry(const int i, const int w, const vector<double> & input);
  void get_ancestry(const int i, const int w, vector<double> & output) const;
  void set_ancestry_global(const int i, const vector<double> & input);
  void get_ancestry_global(const int i, vector<double> & output) const;
  void set_Hk(const int i, vector<int> & hk);
  void set_partition(const int r);
  int num_partitions() const;
  int num_groups() const;
  void load_hmm(const string hmm_file);
  void init_hmm(const int K, const double _rho, const double _lambda);
  void init_hmm(const int K, const double _rho);
  void update_hmm(const double _rho, const vector<double>& _mut_rates);

  // Methods
  void write(const string& basename, const string& format, bool include_genotypes) const;
  void reset_knockoffs();
  void CV(const int num_threads);
  void EM(const int num_threads);

  // Methods for knockoff generation
  void generate(const int num_threads);
  void generate_unrelated(const int num_threads);
  void generate_related(const int num_threads);
  void generate_related_cluster(const vector<int> & related_cluster, boost::random::taus88 & rng);

  // Debug
  void print_hmm() const;
  void writeGroups(const string _out_file_name) const;
  void writeWindows(const string _out_file_name) const;
  void writeAncestries(const string _out_file_name) const;
  void writeHMM(const string _out_file_name) const;
  // void writeZ(const string _out_file_name) const;
  // void writeZk(const string _out_file_name) const;

private:
  const int seed, debug;

  // HMM
  double rho;
  vector <double> genetic_dist, phys_dist;
  vector <double> rec_rates, b;
  vector <double> mut_rates;
  void set_groups(const vector<int> & _groups);
  inline double mutate(const int h, const double mut_rate) const;

  // Output
  vector <chaplotype> Hk;

  // // DEBUG: store Z
  // vector< vector<int> > Z, Zk;  

  // EM methods
  double EM_step(const int num_threads, const bool show_progress);
  void EM_worker(const vector<int> & i_list, boost::random::taus88 & rng, const int wid, 
                 vector<int> &progress, const int n_steps_each, const bool show_progress,
                 vector<double>& gamma, vector<double>& Xi) const;
  double EM_Psi_rho(const double _rho, const vector<double>& Xi) const;
  double EM_Phi_rho(const double _rho) const;
  double EM_solve_rho(const double _rho, const vector<double>& Xi) const;
  void EM_solve_alpha(const vector<double>& eta, const vector<double>& xi, 
                      vector<double>& omega, const int j_start, const int j_end) const;
  
  // Methods for knockoff generation (related)
  void generate_related_worker(const vector< vector<int> > & related_families, boost::random::taus88 & rng,
                               const int wid, vector<int> &progress, const int n_steps_each);

  // Methods for knockoff generation (unrelated)
  void generate_unrelated_worker(const vector<int> & i_list, boost::random::taus88 & rng,
                                 const int wid, vector<int> &progress, const int n_steps_each);

  // Internal algorithms for knockoff generation
  void sample_emission(const vector<int> & z, const vector<int>& ref, vector<int> & output,
                       boost::random::taus88 & rng) const;
  void sample_emission(const vector<int> & z, const vector< vector<int> >& ref, 
                       vector<int> & output, boost::random::taus88 & rng) const;
  void sample_emission(const vector<int> & z, const vector<int>& ref, vector<int> & output,
                       const int j_start, const int j_end, boost::random::taus88 & rng) const;
  void sample_emission_copula(const int i, const vector<int> & z, const vector<int> & zk, vector<int> & output,
                              double lambda, boost::random::taus88 & rng) const;
  void sample_emission_copula(const int i, const vector<int> & z, const vector<int> & zk, const vector<int>& ref,
                              vector<int> & output, double lambda, const int j_start, const int j_end,
                              boost::random::taus88 & rng) const;
  void compute_hmm_forward(const int i, boost::random::taus88 & rng, const vector<int>& ref, 
                           vector< vector<double> >& forward) const;
  void compute_hmm_forward(const int i, boost::random::taus88 & rng, const vector<int>& ref, 
                           vector< vector<double> >& forward, 
                           const int w, const int j_start, const int j_end) const;
  void compute_hmm_backward(const int i, boost::random::taus88 & rng, const vector<int>& ref,
                            vector< vector<double> >& backward) const;
  void compute_hmm_backward(const int i, boost::random::taus88 & rng, const vector<int>& ref,
                            vector< vector<double> >& backward, 
                            const int w, const int j_start, const int j_end) const;
  void sample_posterior_MC(const int i, vector<int> & z, vector< vector<double> >& backward, const int w,
                           boost::random::taus88 & rng) const;
  void sample_knockoff_MC(const int i, const vector<int> & z, vector<int> & output, 
                          const int w, boost::random::taus88 & rng) const;
  void sample_knockoff_MC(const int i, const vector<int> & z, vector<int> & output,
                          boost::random::taus88 & rng) const;
  void estimate_motif_prevalence(const int i, vector<int> & z, 
                                 const vector< vector<double> >& backward, vector< vector<double> >& forward,
                                 vector< vector<double> >& new_alpha, boost::random::taus88 & rng) const;
  
  // Reference set
  const int K;
  ivector3d references_local;
  vector< vector< vector<double> > > alpha;
  ivector2d references_global;
  vector< vector<double> > alpha_global;

  // Genomic windows
  int num_windows;

  // Variant groups
  int num_groups_;
  vector< vector<int> > elements;
  vector<int> groups;
  vector< vector<int> > partitions;

  // DEBUG
  const string logfile;
  void print_debug_log(const int wid, const vector<int> & related_cluster) const;
  void print_debug_log(const int wid, const vector<int> & related_cluster, const bool finished) const;
};

#endif
