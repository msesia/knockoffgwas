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
  FamilyKnockoffs(const vector<unsigned int> & _indices, const IbdCluster& ibd_cluster,
                  const vector< vector<unsigned int> > & _Z, const vector< vector<double> > & _alpha, 
                  const vector<double> & _b, boost::random::taus88 & _rng);

  // Methods
  void run(const vector< vector<unsigned int> > & elements, const vector<unsigned int> & groups,
           vector< vector<unsigned int> > & Zk);

private:

  // Data
  vector<unsigned int> indices;
  vector<double> b;
  unsigned int num_haps, num_snps, K;
  vector< vector< vector<unsigned int > > > neighbors; // (id1,j,id2)
  vector<unsigned int> segments_j_min, segments_j_max;
  vector< vector<unsigned int> > segments_i;

  // Shadow data 
  const vector< vector<unsigned int> > & Z;
  const vector< vector<double> > & alpha;

  // Random number generator
  boost::random::taus88 & rng;
  
  // HMM methods
  void sample_knockoff_MC(unsigned int i, const vector<unsigned int> & z, vector<unsigned int> & zk, unsigned int g_begin, unsigned int g_end,
                          const vector< vector<unsigned int> > & elements, const vector<unsigned int> & groups) const;

  // Debug
  void print_Z(vector< vector<unsigned int> > & Zk) const;
};

#endif
