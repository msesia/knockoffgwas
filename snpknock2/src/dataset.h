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
  Dataset(const vector<Metadata>& _metadata, int _num_threads, int _thinning);

  // Methods
  int num_chrs() const;

  friend class KinshipComputer;

private:
  // Data
  vector<Haplotypes> chromosomes;
  vector<Metadata> metadata;
  int num_haps;
  int num_chrs_;
  int num_threads;

  // Methods
  void load_worker(int chr_min, int chr_max, int wid, vector<int> &progress);
  void check_sanity() const;
  int count_lines(string filename) const;

};

#endif
