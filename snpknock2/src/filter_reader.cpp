#include "filter_reader.h"

filter_reader::filter_reader() {
}

filter_reader::~filter_reader() {
  ind_excluded.clear();
  ind_included.clear();
  snp_excluded.clear();
  snp_included.clear();
}

void filter_reader::initialise(int min_pos, int max_pos) {
  this->min_pos = min_pos;
  this->max_pos = max_pos;
}

bool filter_reader::checkInd(string id) {
  map< string, bool>::iterator it = ind_excluded.find(id);
  if (it != ind_excluded.end()) return false;

  if (ind_included.size() > 0) {
    it = ind_included.find(id);
    if (it == ind_included.end()) return false;
  }
  return true;
}

bool filter_reader::checkSnp(int pos) {
  if (snp_excluded.size() > 0 && snp_excluded.find(pos) != snp_excluded.end())
    return false;
  if (snp_included.size() > 0 && snp_included.find(pos) == snp_included.end())
    return false;
  if (pos >= min_pos && pos < max_pos) return true;
  return false;
}


void filter_reader::readIndExcludeList(string file_ind) {
  string buffer;
  //LOG.println("\nReading individuals to exclude from input file in ["  + file_ind + "]");
  ifile fd_ind(file_ind);
  int cpt_name = 0;
  while (getline(fd_ind, buffer, '\n')) {
    ind_excluded[buffer] = false;
    cpt_name++;
  }
  fd_ind.close();
  //if (cpt_name > 0) LOG.println("  * " + sutils::int2str(cpt_name) + " individuals found in the exclude list");
}

void filter_reader::readIndIncludeList(string file_ind) {
  string buffer;
  //LOG.println("\nReading individuals to include from input file in ["  + file_ind + "]");
  ifile fd_ind(file_ind);
  int cpt_name = 0;
  while (getline(fd_ind, buffer, '\n')) {
    ind_included[buffer] = false;
    cpt_name++;
  }
  fd_ind.close();
  //if (cpt_name > 0) LOG.println("  * " + sutils::int2str(cpt_name) + " individuals found in the include list");
}

void filter_reader::readSnpExcludeList(string file_snp) {
  string buffer;
  //LOG.println("\nReading SNPs to exclude from input file in ["  + file_snp + "]");
  ifile fd_snp(file_snp);
  int cpt_pos = 0;
  while (getline(fd_snp, buffer, '\n')) {
    snp_excluded[atoi(buffer.c_str())] = true;
    cpt_pos ++;
  }
  fd_snp.close();
  //if (cpt_pos > 0) LOG.println("  * " + sutils::int2str(cpt_pos) + " snps found in the exclude list");
}

void filter_reader::readSnpIncludeList(string file_snp) {
  string buffer;
  //LOG.println("\nReading SNPs to include from input file in ["  + file_snp + "]");
  ifile fd_snp(file_snp);
  int cpt_pos = 0;
  while (getline(fd_snp, buffer, '\n')) {
    snp_included[atoi(buffer.c_str())] = true;
    cpt_pos ++;
  }
  fd_snp.close();
  //if (cpt_pos > 0) LOG.println("  * " + sutils::int2str(cpt_pos) + " snps found in the include list");
}
