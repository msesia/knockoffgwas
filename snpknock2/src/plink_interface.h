#ifndef PLINK_INTERFACE_H
#define PLINK_INTERFACE_H

#include <vector>
#include <string>
#include <fstream>
#include <random>
#include "metadata.h"
#include "chaplotype.h"

using namespace std;

namespace plink {
  void write_binary(const string& basename, const Metadata& metadata, const vector<chaplotype>& Hk);
  void write_binary(const string& basename, const Metadata& metadata, 
                    const vector<chaplotype>& H, const vector<chaplotype>& Hk, bool random);
  void write_haps(const string& basename, const Metadata& metadata, const vector<chaplotype>& Hk);  
  void writeSAMPLE(const string& basename, const Metadata& metadata);
  }

#endif
