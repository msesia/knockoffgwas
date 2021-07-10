#ifndef _HAPLOTYPE_H
#define _HAPLOTYPE_H

#include "utils.h"
#include "chaplotype.h"
#include "filter_reader.h"
#include "metadata.h"
#include "haps_reader.h"
#include "bgen_reader.h"

class Haplotypes {
public:
  // Constructors
  Haplotypes(const Metadata& _metadata);
  ~Haplotypes();

  // Methods
  void load_data(bool verbose, int nthreads);
  void writeHaplotypes(string filename) const;
  unsigned int rpartition(double w, vector<unsigned int> & output);
  unsigned int distance(const unsigned int i1, const unsigned int i2) const;
  unsigned int distance(const unsigned int i1, const unsigned int i2, const unsigned int j_start, const unsigned int j_end) const;
  void compute_maf(vector<double>& maf) const;

  // Get/set
  unsigned int get_num_snps() const;
  unsigned int get_num_haps() const;
  void push_back(const vector<unsigned int> & input);
  void push_back(const vector<bool> & input);

  // Data
  vector <chaplotype> H;

  friend class KinshipComputer;
	
protected:	
  unsigned int num_snps;
  unsigned int num_haps;
  Metadata metadata;

};

#endif
