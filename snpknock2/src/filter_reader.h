#ifndef _FILTER_READER_H
#define _FILTER_READER_H

#include "utils.h"

class filter_reader {
public:
  //DATA
  map < string , bool > ind_excluded;
  map < string , bool > ind_included;
  map < int , bool > snp_excluded;
  map < int , bool > snp_included;
  int min_pos;
  int max_pos;

  //CONSTRUCTOR/DESTRUCTOR
  filter_reader();
  ~filter_reader();

  //METHODS
  void initialise(int, int);
  bool checkInd(string);
  bool checkSnp(int);

  //IO
  void readIndExcludeList(string);
  void readIndIncludeList(string);
  void readSnpExcludeList(string);
  void readSnpIncludeList(string);
};

#endif
