#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include "utils.h"
#include "metadata.h"

using namespace std;

class Arguments {
public:
  //Constructors
  Arguments(); //Default constructor
  explicit Arguments(unsigned int arg_num, const char* argslist[]);
  void check() const;

  //Accessors
  vector<DataFiles> get_filenames() const;
  DataFiles get_filenames(unsigned int chr) const;
  string get_filename(string filetype) const;

  string get_pc_file() const;
  string get_output_file(unsigned int chr) const;
  string get_log_file(unsigned int chr) const;
  vector<string> get_output_files() const;
  vector<string> get_log_files() const;
  vector<string> get_hmm_files() const;
  string get_hmm_file(unsigned int chr) const;
  vector<string> get_ibd_files() const;
  string get_ibd_file(unsigned int chr) const;
  vector<string> get_ref_files() const;
  string get_ref_file(unsigned int chr) const;
  vector<string> get_lref_files() const;
  string get_lref_file(unsigned int chr) const;
  unsigned int num_threads() const;
  unsigned int get_debug() const;
  unsigned int get_cluster_size_min() const;
  unsigned int get_cluster_size_max() const;
  unsigned int get_K() const;
  unsigned int get_seed() const;
  double get_hmm_rho() const;
  double get_hmm_lambda() const;
  bool get_compute_kinship() const;
  bool get_estimate_hmm() const;
  bool get_generate_knockoffs() const;
  unsigned int num_chrs() const;
  unsigned int get_resolution() const;
  double get_window_size() const;
private:

  string keep_file, pc_file;
  string data_format;
  vector<string> data_files;
  vector<string> map_files, out_files, log_files, part_files, hmm_files, ref_files, lref_files, ibd_files;
  vector<string> extract_files;
  unsigned int cluster_size_min, cluster_size_max;
  const unsigned int cluster_size_min_default = 1000;
  const unsigned int cluster_size_max_default = 2500;
  double hmm_rho, hmm_lambda;
  unsigned int n_threads, debug;
  unsigned int K;
  unsigned int seed;
  const unsigned int K_default = 100;
  bool compute_kinship, estimate_hmm, generate_knockoffs;
  set<string> operations;
  unsigned int num_chrs_;
  unsigned int resolution;
  double window_size;

  void parse_args(unsigned int arg_num,const char* argslist[]);
  void parse_file_str(string file_name, string state);
  void parse_wildcard_str(string file_name, string state);
  void parse_n_threads_option(string str);
  void parse_resolution_option(string str);
  void parse_debug_option(string str);
  void parse_cluster_size_option(string str, string state);
  void parse_K_option(string K_str);
  void parse_seed_option(string seed_str);
  void parse_hmm_rho_option(string input);
  void parse_hmm_lambda_option(string input);
  void parse_window_size_option(string input);
  void check_file(string input_file) const;
  void print_help() const;
  
};

#endif
