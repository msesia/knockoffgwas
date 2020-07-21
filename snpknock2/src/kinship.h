#ifndef _KINSHIP_H
#define _KINSHIP_H

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include "kinship_computer.h"
#include "utils.h"
#include "metadata.h"
#include "covariates.h"

typedef vector< vector<int> > ivector2d;
typedef vector<ivector2d> ivector3d;
typedef vector<ivector3d> ivector4d;

class Kinship {
public:
  // Constructors and destructors
  Kinship();
  Kinship(const vector<Metadata>& metadata, vector<string> output_files,
          int compression, int _cluster_size_min, int _cluster_size_max, int _K, int _num_threads);
  Kinship(const vector<Metadata>& metadata, const Covariates & covariates, vector<string> output_files,
          int compression, int _cluster_size_min, int _cluster_size_max, int _K, int _num_threads);
  Kinship(const vector<string>& ref_files_local, const vector<string>& ref_files_global);
  ~Kinship();

  // Methods
  const ivector3d & get_references(int chr) const;
  const ivector2d & get_references_global(int chr) const;

  int get_K() const;
  void writeReferences(const vector<Metadata>& metadata, const vector<string> & out_file_names) const;

private:
  int num_haps;
  int num_chrs;
  int num_windows;
  int K;
  ivector4d references_local;  // CHR, id, window, K
  ivector3d references_global; // CHR, id, K

  void sanitycheck_references() const;
  void load_references_local(const vector<string>& ref_files);
  void load_references_global(const vector<string>& ref_files);
  void writeReferencesLocal(const vector<Metadata>& metadata, const vector<string> & out_file_names) const;
  void writeReferencesGlobal(const vector<Metadata>& metadata, const vector<string> & out_file_names) const;

};

#endif
