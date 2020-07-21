#ifndef _HAPS_READER_H
#define _HAPS_READER_H

#include <iostream>
#include <fstream>
#include "chaplotype.h"
#include "metadata.h"

using namespace std;

class HapsReader {
public:
  HapsReader(const string& filename, const string& sample_filename, const string& legend_filename);
  void read(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H, bool verbose);

private:
  Sample sample;
  Legend legend;
  string filename;
};

#endif
