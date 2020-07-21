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
  int rpartition(double w, vector<int> & output);
  int distance(const int i1, const int i2) const;
  int distance(const int i1, const int i2, const int j_start, const int j_end) const;
  void compute_maf(vector<double>& maf) const;

  // Get/set
  int get_num_snps() const;
  int get_num_haps() const;
  void push_back(const vector<int> & input);
  void push_back(const vector<bool> & input);

  // Data
  vector <chaplotype> H;

  friend class KinshipComputer;
	
protected:	
  int num_snps;
  int num_haps;
  Metadata metadata;

};

#endif
