#ifndef COVARIATES_H
#define COVARIATES_H

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "stdafx.h"
#include "utils.h"
#include "metadata.h"

using namespace std;

class Covariates {
public:
  // Constructors/destructors
  Covariates();
  Covariates(const string & pc_filename, const Metadata & metadata);
  Covariates(const Covariates& obj);
  ~Covariates();
  // Data
  vector< vector<double> > Z;
  Sample sample;
  int K, K_max;

  void print() const;

private:
  // Methods (TODO)
  void clear();
  void filter(const Metadata & metadata, Sample & sample_raw, vector< vector<double> > & Z_raw);
  void mask_by_family(const Metadata & metadata);

};

void load_data(const string & filename, int K_max, Sample & sample_raw, vector< vector<double> > & Z_raw);

#endif
