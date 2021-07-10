#ifndef WINDOWS_H
#define WINDOWS_H

#include <vector>
#include <string>
#include "utils.h"

class Windows {
  public:
  Windows();
  Windows(unsigned int _num_snps);
  Windows(const Windows& obj);

  void load(const vector<unsigned int>& start_points, const unsigned int _num_snps);
  void print() const;
    
  unsigned int num_snps, num_windows;
  vector<unsigned int> start, end; 
  vector<unsigned int> window;
};

#endif
