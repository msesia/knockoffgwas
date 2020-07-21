//$Id: chaplotype.h 597 2012-07-18 14:00:06Z koskos $

#ifndef _CHAPLOTYPE_H
#define _CHAPLOTYPE_H

#include "utils.h"

class chaplotype {
public:
  //DATA
  vector < unsigned char > data;
  unsigned int n;

  //CONSTRUCTOR/DESTRUCTOR
  chaplotype();
  chaplotype(unsigned int n);
  chaplotype(const chaplotype & b);
  chaplotype(const vector <bool> & V);
  chaplotype(const vector <int> & V);
  ~chaplotype();
  void clear();

  //METHODS
  void operator = (chaplotype b);
  inline bool operator [] (const int i) const {return(data[i >> 3] & (1 << (i & 7)));}
  void set(int i, bool b);
  void push_back(vector <bool> & V);
  unsigned int size() const;

  inline bool get(int i) const { return data[i >> 3] & (1 << (i & 7)); }

  int hamming(const chaplotype & b) const;
  int hamming(const chaplotype & b, int from, int to) const;
  int hamming(const chaplotype & b, int from, int to, int max) const;
	
};


#endif
