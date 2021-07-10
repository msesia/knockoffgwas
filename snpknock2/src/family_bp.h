#ifndef _FAMILY_BP_H
#define _FAMILY_BP_H

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include "haplotypes.h"
#include "utils.h"
#include "metadata.h"
#include "ibd.h"

typedef vector< vector<unsigned int> > vector2i;
typedef vector< vector2i > vector3i;
typedef vector< vector<double> > vector2d;
typedef vector< vector2d > vector3d;
typedef vector< vector3d > vector4d;

using namespace std;

class FamilyBP {
public:
  // Constructors/destructors
  FamilyBP(const vector<unsigned int> & _indices, const IbdCluster& ibd_cluster,
           const vector <chaplotype> & _H, const vector< vector<unsigned int> > & _ref,
           const vector< vector<double> > & _alpha, const vector <double> & _mut_rates,
           const vector<double> & _b, const Metadata& _metadata, boost::random::taus88 & _rng);

  int run(vector2i & Z);

private:
  // Data
  const unsigned int K, num_haps, num_snps;
  const vector<unsigned int> indices;
  const vector<double> b;
  vector3i neighbors; // (id1,j,id2)
  set<pair<unsigned int,unsigned int>> inactive_nodes;
  vector <double> genetic_dist;

  // Shadow data
  const vector <chaplotype> & H;
  const vector2i & ref;
  const vector2d & alpha;
  const vector <double> & mut_rates;
  const Metadata& metadata; // DEBUG

  // BP messages
  vector3d messages_right; // (id, j, k)
  vector3d messages_left;  // (id, j, k)

  // Random number generator
  boost::random::taus88 & rng;

  // Methods
  void init_messages();
  int run_bp(bool verbose);
  double send_messages_right(const unsigned int i1, const unsigned int j, vector<double> & message, vector<double>& f_j,
                             vector<unsigned int>& neigh_diff_left, vector<unsigned int>& neigh_diff_right) const;
  double send_messages_left(const unsigned int i1, const unsigned int j, vector<double> & message, vector<double>& f_j,
                            vector<unsigned int>& neigh_diff_left, vector<unsigned int>& neigh_diff_right) const;
  void condition_on(const unsigned int i, const unsigned int j, const unsigned int k, const bool forward);
  void condition_on(const unsigned int i, const unsigned int j, const unsigned int k);
  void sample_posterior(vector2i & Z);
  unsigned int sample_posterior_joint(const unsigned int i, const unsigned int j, const vector<unsigned int>& neigh_diff_left,
                             const vector<unsigned int>& neigh_diff_right) const;
  void get_marginals(const unsigned int i1, const unsigned int j, vector<double>& f_j, vector<double>& marginals) const;
  bool is_inactive(const unsigned int i, const unsigned int j) const;

  // Debug
  void print_Z(const vector2i & Z) const;
};

inline double emission_prob(const unsigned int h1, const unsigned int h2, const double mut_rate);

#endif
