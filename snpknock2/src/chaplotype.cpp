#include "chaplotype.h"

unsigned int distance_lookup[256] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };

chaplotype::chaplotype() {
  this-> n = 0;
}

chaplotype::chaplotype(unsigned int n) {
  this-> n = n;
  data = vector < unsigned char > ((unsigned int)ceil(n * 1.0 / 8), 0);
}

chaplotype::chaplotype(const chaplotype & b) {
  n = b.n;
  data = b.data;
}

chaplotype::chaplotype(const vector <bool> & V) {
  n = V.size();
  data = vector < unsigned char > ((unsigned int)ceil(n * 1.0 / 8), 0);
  for (unsigned int i = 0 ; i < n ; i ++) set(i, V[i]);
}

chaplotype::chaplotype(const vector <unsigned int> & V) {
  n = V.size();
  data = vector < unsigned char > ((unsigned int)ceil(n * 1.0 / 8), 0);
  for (unsigned int i = 0 ; i < n ; i ++) set(i, (bool)(V[i]));
}

chaplotype::~chaplotype() {
  data.clear();
  n = 0;
}

void chaplotype::clear() {
  vector < unsigned char > ().swap(data);
  n = 0;
}

void chaplotype::operator = (chaplotype b) {
  n = b.n;
  data = b.data;
}

void chaplotype::set(unsigned int i, bool b) {
  if (b) data[i >> 3] |= (1 << (i & 7));
  else data[i >> 3] &= ~(1 << (i & 7));
}

void chaplotype::push_back(vector < bool > & V) {
  unsigned int offset = n;
  n += V.size();
  data.resize((int)ceil(n * 1.0 / 8), 0);
  for (unsigned int i = 0 ; i < V.size() ; i ++) set(i+offset, V[i]);
}

unsigned int chaplotype::size() const {
  return n;
}

unsigned int chaplotype::hamming(const chaplotype & b) const {
  unsigned int c = 0;
  for (unsigned int i = 0 ; i < data.size() ; i ++) c+= distance_lookup[data[i] ^ b.data[i]];
  return c;
}


unsigned int chaplotype::hamming(const chaplotype & b, unsigned int from, unsigned int to) const{
  unsigned int c = 0;
  unsigned int bfrom = from >> 3;
  unsigned int bto =  to >> 3;
  for (unsigned int i = bfrom ; i <= bto ; i ++) c+= distance_lookup[data[i] ^ b.data[i]];
  return c;
}

unsigned int chaplotype::hamming(const chaplotype & b, unsigned int from, unsigned int to, unsigned int max) const {
  unsigned int c = 0;
  unsigned int bfrom = from >> 3;
  unsigned int bto =  to >> 3;
  for (unsigned int i = bfrom ; i <= bto ; i ++) {
    c+= distance_lookup[data[i] ^ b.data[i]];
    if (c > max) return BIG_POS_INT;
  }
  return c;
}
