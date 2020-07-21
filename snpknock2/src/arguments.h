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
  explicit Arguments(int arg_num, const char* argslist[]);
  void check() const;

  //Accessors
  vector<DataFiles> get_filenames() const;
  DataFiles get_filenames(int chr) const;
  string get_filename(string filetype) const;

  string get_pc_file() const;
  string get_output_file(int chr) const;
  string get_log_file(int chr) const;
  vector<string> get_output_files() const;
  vector<string> get_log_files() const;
  vector<string> get_hmm_files() const;
  string get_hmm_file(int chr) const;
  vector<string> get_ibd_files() const;
  string get_ibd_file(int chr) const;
  vector<string> get_ref_files() const;
  string get_ref_file(int chr) const;
  vector<string> get_lref_files() const;
  string get_lref_file(int chr) const;
  int num_threads() const;
  int get_debug() const;
  int get_cluster_size_min() const;
  int get_cluster_size_max() const;
  int get_K() const;
  int get_seed() const;
  double get_hmm_rho() const;
  double get_hmm_lambda() const;
  bool get_compute_kinship() const;
  bool get_estimate_hmm() const;
  bool get_generate_knockoffs() const;
  int num_chrs() const;
  int get_resolution() const;
  double get_window_size() const;
private:

  string keep_file, pc_file;
  string data_format;
  vector<string> data_files;
  vector<string> map_files, out_files, log_files, part_files, hmm_files, ref_files, lref_files, ibd_files;
  vector<string> extract_files;
  int cluster_size_min, cluster_size_max;
  const int cluster_size_min_default = 1000;
  const int cluster_size_max_default = 2500;
  double hmm_rho, hmm_lambda;
  int n_threads, debug;
  int K;
  int seed;
  const int K_default = 100;
  bool compute_kinship, estimate_hmm, generate_knockoffs;
  set<string> operations;
  int num_chrs_;
  int resolution;
  double window_size;

  void parse_args(int arg_num,const char* argslist[]);
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
