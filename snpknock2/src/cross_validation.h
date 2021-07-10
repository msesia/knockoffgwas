#ifndef _CROSS_VALIDATION_H
#define _CROSS_VALIDATION_H

#include <fstream>
#include <functional>
#include <boost/algorithm/clamp.hpp>
#include "utils.h"
#include "haplotypes.h"

class CrossVal {

 public:
  // Constructors
  CrossVal(const vector <chaplotype> & _H, const vector< vector<unsigned int> > & _ref,
           const vector< vector<double> > & _alpha, const vector <double> & _b,
           const vector <double> & _phys_dist, unsigned int _seed);

  // Methods
  pair<double,double> fit(unsigned int num_threads);

 private:

  // Data
  unsigned int num_haps, num_snps, K, seed;
  const vector <chaplotype> & H;
  const vector< vector<unsigned int> > & ref;
  const vector< vector<double> > & alpha;
  const vector <double> & b;
  const vector <double> & phys_dist;

  // Methods  
  void CV_errors(vector<double> & errors,
                 const vector<bool> & mask, const vector< vector<unsigned int> > & masked_blocks,
                 const vector <double> & b_masked, const vector <double> & mut_rates,
                 unsigned int num_threads, boost::random::taus88 & rng) const;

  void CV_errors_worker(const vector<unsigned int> & i_list, vector<double> & errors,
                        const vector<bool> & mask, const vector< vector<unsigned int> > & masked_blocks,
                        const vector <double> & b_masked, const vector <double> & mut_rates,
                        boost::random::taus88 & rng
                        ) const;

  void sample_posterior_MC(unsigned int i, vector<unsigned int> & z, const vector<bool> & mask, 
                           const vector <double> & b_masked, const vector <double> & mut_rates,
                           boost::random::taus88 & rng) const;

  void impute_MC(unsigned int i, vector<unsigned int> & z,
                 const vector<bool> & mask, const vector< vector<unsigned int> > & masked_blocks,
                 boost::random::taus88 & rng) const;

  void impute_HMM(unsigned int i, vector<unsigned int> & x, const vector<unsigned int> &z, const vector<bool> & mask,
                  const vector <double> & mut_rates, boost::random::taus88 & rng) const;

  void init_mask(vector<bool> & mask, vector< vector<unsigned int> > & masked_blocks,
                 double q_start, double q_stop, boost::random::taus88 & rng) const;

  void apply_mask(const vector<bool> & mask, double rho, double lambda,
                  vector <double> & b_masked, vector <double> & mut_rates);

  double mutate(unsigned int h, double mut_rate) const;
};

void log_range(vector<double> & values, double min_value, double max_value, unsigned int length);
double compute_sd(const vector<double> & data);
double compute_mean(const vector<double> & data);
unsigned int min_1sd(const vector<double> & mean, const vector<double> & sd, bool prefer_larger);
pair<unsigned int,unsigned int> range_1sd(const vector<double> & mean, const vector<double> & sd);
double init_mut_rate(unsigned int K);

#endif
