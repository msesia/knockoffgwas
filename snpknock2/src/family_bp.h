#ifndef _FAMILY_BP_H
#define _FAMILY_BP_H

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include "haplotypes.h"
#include "utils.h"
#include "metadata.h"
#include "ibd.h"

typedef vector< vector<int> > vector2i;
typedef vector< vector2i > vector3i;
typedef vector< vector<double> > vector2d;
typedef vector< vector2d > vector3d;
typedef vector< vector3d > vector4d;

using namespace std;

class FamilyBP {
public:
  // Constructors/destructors
  FamilyBP(const vector<int> & _indices, const IbdCluster& ibd_cluster,
           const vector <chaplotype> & _H, const vector< vector<int> > & _ref,
           const vector< vector<double> > & _alpha, const vector <double> & _mut_rates,
           const vector<double> & _b, const Metadata& _metadata, boost::random::taus88 & _rng);

  int run(vector2i & Z);

private:
  // Data
  const int K, num_haps, num_snps;
  const vector<int> indices;
  const vector<double> b;
  vector3i neighbors; // (id1,j,id2)
  set<pair<int,int>> inactive_nodes;
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
  double send_messages_right(const int i1, const int j, vector<double> & message, vector<double>& f_j,
                             vector<int>& neigh_diff_left, vector<int>& neigh_diff_right) const;
  double send_messages_left(const int i1, const int j, vector<double> & message, vector<double>& f_j,
                            vector<int>& neigh_diff_left, vector<int>& neigh_diff_right) const;
  void condition_on(const int i, const int j, const int k, const bool forward);
  void condition_on(const int i, const int j, const int k);
  void sample_posterior(vector2i & Z);
  int sample_posterior_joint(const int i, const int j, const vector<int>& neigh_diff_left,
                             const vector<int>& neigh_diff_right) const;
  void get_marginals(const int i1, const int j, vector<double>& f_j, vector<double>& marginals) const;
  bool is_inactive(const int i, const int j) const;

  // Debug
  void print_Z(const vector2i & Z) const;
};

inline double emission_prob(const int h1, const int h2, const double mut_rate);

#endif
