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
  chaplotype(const vector <unsigned int> & V);
  ~chaplotype();
  void clear();

  //METHODS
  void operator = (chaplotype b);
  inline bool operator [] (const unsigned int i) const {return(data[i >> 3] & (1 << (i & 7)));}
  void set(unsigned int i, bool b);
  void push_back(vector <bool> & V);
  unsigned int size() const;

  inline bool get(unsigned int i) const { return data[i >> 3] & (1 << (i & 7)); }

  unsigned int hamming(const chaplotype & b) const;
  unsigned int hamming(const chaplotype & b, unsigned int from, unsigned int to) const;
  unsigned int hamming(const chaplotype & b, unsigned int from, unsigned int to, unsigned int max) const;
	
};


#endif
