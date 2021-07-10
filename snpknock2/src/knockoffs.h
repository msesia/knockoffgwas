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

typedef vector< vector<unsigned int> > ivector2d;
typedef vector<ivector2d> ivector3d;

class Knockoffs : public Haplotypes {

public:
  // Constructors
  Knockoffs(const Metadata& _metadata, unsigned int _K, unsigned int nthreads, unsigned int _debug, unsigned int _seed, const string & _logfile);
  Knockoffs(const Metadata& _metadata, unsigned int _K, bool _load_data, unsigned int nthreads, unsigned int _debug, unsigned int _seed,
            const string & _logfile);
  void load(const unsigned int nthreads);

  // Get/set
  void set_reference(const unsigned int i, const ivector2d & input);
  void set_references(const ivector3d& input);
  void get_reference(const unsigned int i, ivector2d& output) const;
  void set_reference_global(const unsigned int i, const vector<unsigned int>& input);
  void set_references_global(const ivector2d& input);
  void get_reference_global(const unsigned int i, vector<unsigned int>& output) const;
  void set_ancestry(const unsigned int i, const unsigned int w, const vector<double> & input);
  void get_ancestry(const unsigned int i, const unsigned int w, vector<double> & output) const;
  void set_ancestry_global(const unsigned int i, const vector<double> & input);
  void get_ancestry_global(const unsigned int i, vector<double> & output) const;
  void set_Hk(const unsigned int i, vector<unsigned int> & hk);
  void set_partition(const unsigned int r);
  unsigned int num_partitions() const;
  unsigned int num_groups() const;
  void load_hmm(const string hmm_file);
  void init_hmm(const unsigned int K, const double _rho, const double _lambda);
  void init_hmm(const unsigned int K, const double _rho);
  void update_hmm(const double _rho, const vector<double>& _mut_rates);

  // Methods
  void write(const string& basename, const string& format, bool include_genotypes) const;
  void reset_knockoffs();
  void CV(const unsigned int num_threads);
  void EM(const unsigned int num_threads);

  // Methods for knockoff generation
  void generate(const unsigned int num_threads);
  void generate_unrelated(const unsigned int num_threads);
  void generate_related(const unsigned int num_threads);
  void generate_related_cluster(const vector<unsigned int> & related_cluster, boost::random::taus88 & rng);

  // Debug
  void print_hmm() const;
  void writeGroups(const string _out_file_name) const;
  void writeWindows(const string _out_file_name) const;
  void writeAncestries(const string _out_file_name) const;
  void writeHMM(const string _out_file_name) const;
  // void writeZ(const string _out_file_name) const;
  // void writeZk(const string _out_file_name) const;

private:
  const unsigned int seed, debug;

  // HMM
  double rho;
  vector <double> genetic_dist, phys_dist;
  vector <double> rec_rates, b;
  vector <double> mut_rates;
  void set_groups(const vector<unsigned int> & _groups);
  inline double mutate(const unsigned int h, const double mut_rate) const;

  // Output
  vector <chaplotype> Hk;

  // // DEBUG: store Z
  // vector< vector<unsigned int> > Z, Zk;  

  // EM methods
  double EM_step(const unsigned int num_threads, const bool show_progress);
  void EM_worker(const vector<unsigned int> & i_list, boost::random::taus88 & rng, const unsigned int wid, 
                 vector<unsigned int> &progress, const unsigned int n_steps_each, const bool show_progress,
                 vector<double>& gamma, vector<double>& Xi) const;
  double EM_Psi_rho(const double _rho, const vector<double>& Xi) const;
  double EM_Phi_rho(const double _rho) const;
  double EM_solve_rho(const double _rho, const vector<double>& Xi) const;
  void EM_solve_alpha(const vector<double>& eta, const vector<double>& xi, 
                      vector<double>& omega, const unsigned int j_start, const unsigned int j_end) const;
  
  // Methods for knockoff generation (related)
  void generate_related_worker(const vector< vector<unsigned int> > & related_families, boost::random::taus88 & rng,
                               const unsigned int wid, vector<unsigned int> &progress, const unsigned int n_steps_each);

  // Methods for knockoff generation (unrelated)
  void generate_unrelated_worker(const vector<unsigned int> & i_list, boost::random::taus88 & rng,
                                 const unsigned int wid, vector<unsigned int> &progress, const unsigned int n_steps_each);

  // Internal algorithms for knockoff generation
  void sample_emission(const vector<unsigned int> & z, const vector<unsigned int>& ref, vector<unsigned int> & output,
                       boost::random::taus88 & rng) const;
  void sample_emission(const vector<unsigned int> & z, const vector< vector<unsigned int> >& ref, 
                       vector<unsigned int> & output, boost::random::taus88 & rng) const;
  void sample_emission(const vector<unsigned int> & z, const vector<unsigned int>& ref, vector<unsigned int> & output,
                       const unsigned int j_start, const unsigned int j_end, boost::random::taus88 & rng) const;
  void sample_emission_copula(const unsigned int i, const vector<unsigned int> & z, const vector<unsigned int> & zk, vector<unsigned int> & output,
                              double lambda, boost::random::taus88 & rng) const;
  void sample_emission_copula(const unsigned int i, const vector<unsigned int> & z, const vector<unsigned int> & zk, const vector<unsigned int>& ref,
                              vector<unsigned int> & output, double lambda, const unsigned int j_start, const unsigned int j_end,
                              boost::random::taus88 & rng) const;
  void compute_hmm_forward(const unsigned int i, boost::random::taus88 & rng, const vector<unsigned int>& ref, 
                           vector< vector<double> >& forward) const;
  void compute_hmm_forward(const unsigned int i, boost::random::taus88 & rng, const vector<unsigned int>& ref, 
                           vector< vector<double> >& forward, 
                           const unsigned int w, const unsigned int j_start, const unsigned int j_end) const;
  void compute_hmm_backward(const unsigned int i, boost::random::taus88 & rng, const vector<unsigned int>& ref,
                            vector< vector<double> >& backward) const;
  void compute_hmm_backward(const unsigned int i, boost::random::taus88 & rng, const vector<unsigned int>& ref,
                            vector< vector<double> >& backward, 
                            const unsigned int w, const unsigned int j_start, const unsigned int j_end) const;
  void sample_posterior_MC(const unsigned int i, vector<unsigned int> & z, vector< vector<double> >& backward, const unsigned int w,
                           boost::random::taus88 & rng) const;
  void sample_knockoff_MC(const unsigned int i, const vector<unsigned int> & z, vector<unsigned int> & output, 
                          const unsigned int w, boost::random::taus88 & rng) const;
  void sample_knockoff_MC(const unsigned int i, const vector<unsigned int> & z, vector<unsigned int> & output,
                          boost::random::taus88 & rng) const;
  void estimate_motif_prevalence(const unsigned int i, vector<unsigned int> & z, 
                                 const vector< vector<double> >& backward, vector< vector<double> >& forward,
                                 vector< vector<double> >& new_alpha, boost::random::taus88 & rng) const;
  
  // Reference set
  const unsigned int K;
  ivector3d references_local;
  vector< vector< vector<double> > > alpha;
  ivector2d references_global;
  vector< vector<double> > alpha_global;

  // Genomic windows
  unsigned int num_windows;

  // Variant groups
  unsigned int num_groups_;
  vector< vector<unsigned int> > elements;
  vector<unsigned int> groups;
  vector< vector<unsigned int> > partitions;

  // DEBUG
  const string logfile;
  void print_debug_log(const unsigned int wid, const vector<unsigned int> & related_cluster) const;
  void print_debug_log(const unsigned int wid, const vector<unsigned int> & related_cluster, const bool finished) const;
};

#endif
