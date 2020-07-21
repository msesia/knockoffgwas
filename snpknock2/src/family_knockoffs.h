#ifndef _FAMILY_KNOCKOFFS_H
#define _FAMILY_KNOCKOFFS_H

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include "haplotypes.h"
#include "utils.h"
#include "ibd.h"

class FamilyKnockoffs {
public:
  // Constructors/destructors
  FamilyKnockoffs(const vector<int> & _indices, const IbdCluster& ibd_cluster,
                  const vector< vector<int> > & _Z, const vector< vector<double> > & _alpha, 
                  const vector<double> & _b, boost::random::taus88 & _rng);

  // Methods
  void run(const vector< vector<int> > & elements, const vector<int> & groups,
           vector< vector<int> > & Zk);

private:

  // Data
  vector<int> indices;
  vector<double> b;
  int num_haps, num_snps, K;
  vector< vector< vector<int > > > neighbors; // (id1,j,id2)
  vector<int> segments_j_min, segments_j_max;
  vector< vector<int> > segments_i;

  // Shadow data 
  const vector< vector<int> > & Z;
  const vector< vector<double> > & alpha;

  // Random number generator
  boost::random::taus88 & rng;
  
  // HMM methods
  void sample_knockoff_MC(int i, const vector<int> & z, vector<int> & zk, int g_begin, int g_end,
                          const vector< vector<int> > & elements, const vector<int> & groups) const;

  // Debug
  void print_Z(vector< vector<int> > & Zk) const;
};

#endif
