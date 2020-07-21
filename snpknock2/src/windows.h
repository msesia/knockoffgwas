#ifndef WINDOWS_H
#define WINDOWS_H

#include <vector>
#include <string>
#include "utils.h"

class Windows {
  public:
  Windows();
  Windows(int _num_snps);
  Windows(const Windows& obj);

  void load(const vector<int>& start_points, const int _num_snps);
  void print() const;
    
  int num_snps, num_windows;
  vector<int> start, end; 
  vector<int> window;
};

#endif
