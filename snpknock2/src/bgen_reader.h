#ifndef _BGEN_READER_H
#define _BGEN_READER_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <boost/thread.hpp>
#include <numeric>
#include "bgen_utils.h"
#include "chaplotype.h"
#include "metadata.h"
#include "utils.h"

using namespace std;

class BgenReader {
public:
  BgenReader(const string& filename, const string& sample_filename);
  void summarise();
  void read(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H, 
            bool verbose, int nthreads);
  void read(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H, 
            bool verbose);
  void read_worker(const vector<string>& sample_ids, const vector<string>& abs_snp_ids, 
                   const vector<string>& snp_ids, vector<chaplotype>& H, 
                   bool verbose, int wid, vector<int> &progress);
  void read_mt(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H, 
               bool verbose, int nthreads);
  void print(const vector<chaplotype>& H) const;

private:
  std::pair<bool, bool> prob_to_hap(const vector<double>& probs) const;
  string filename;
  Sample sample;
};

#endif
