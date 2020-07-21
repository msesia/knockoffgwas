#ifndef WINDOWS_CPP
#define WINDOWS_CPP

#include "windows.h"

Windows::Windows() {
}

Windows::Windows(int _num_snps) {
  num_snps = _num_snps;
  num_windows = 1;
  start = {0};
  end = {num_snps};
  window.resize(num_snps, 0);
}

Windows::Windows(const Windows& obj) {
  num_snps = obj.num_snps;
  num_windows = obj.num_windows;
  start = obj.start;
  end = obj.end;
  window = obj.window;
}


void Windows::load(const vector<int>& start_points, const int _num_snps) {
  num_snps = _num_snps;
  num_windows = start_points.size();
  window.resize(num_snps, 0);

  for(int w=0; w<num_windows; w++) {
    start.push_back(start_points[w]);
    if(w<(num_windows-1)) end.push_back(start_points[w+1]);
    else end.push_back(num_snps);
    for(int j=start.back(); j<end.back(); j++) {
      window[j] = w;
    }        
  } 
}


void Windows::print() const {
  cout << "Printing summary of " << num_windows << " windows:" << endl;
  for(int w=0; w<num_windows; w++) {
    cout << std::setw(6) << w << ": " << start[w] << "--" << end[w] << endl;
  }
}

#endif
