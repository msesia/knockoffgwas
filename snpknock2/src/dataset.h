#ifndef _DATASET_H
#define _DATASET_H

#include <numeric>
#include <boost/thread.hpp>
#include "filter_reader.h"
#include "haplotypes.h"
#include "metadata.h"

class Dataset {
public:
  // Constructors/desctructors
  Dataset();
  Dataset(const vector<Metadata>& _metadata, unsigned int _num_threads, unsigned int _thinning);

  // Methods
  unsigned int num_chrs() const;

  friend class KinshipComputer;

private:
  // Data
  vector<Haplotypes> chromosomes;
  vector<Metadata> metadata;
  unsigned int num_haps;
  unsigned int num_chrs_;
  unsigned int num_threads;

  // Methods
  void load_worker(unsigned int chr_min, unsigned int chr_max, unsigned int wid, vector<unsigned int> &progress);
  void check_sanity() const;
  unsigned int count_lines(string filename) const;

};

#endif
