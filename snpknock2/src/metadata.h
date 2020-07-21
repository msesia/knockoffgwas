#ifndef METADATA_H
#define METADATA_H

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "stdafx.h"
#include "interpolation.h"
#include "ibd.h"
#include "utils.h"
#include "utils_graph.h"
#include "windows.h"

using namespace std;

struct DataFiles {
  string data = "";         // File prefix containing haplotype data
  string format = "haps";   // Format of the haplotype file ("haps" or "bgen")
  string sample = "";       // File containing sample and family IDs
  string legend = "";       // File containing legend of variants
  string map = "";          // File containing map of variant positions
  string keep = "";         // File containing list of sample IDs for filter
  string extract = "";      // File containing list of variant IDs for filter
  string partitions = "";   // File containing list of variant partitions at different resolutions
  string ibd = "";          // File containing list of IBD segments
};

class Sample {
public:
  // Constructors/destructors
  Sample();
  Sample(const Sample& obj);
  ~Sample();

  // Methods
  void clear();
  void push_back(vector<string> row);
  vector<string> row(int idx) const;
  int size() const;
  void filter(const vector<string>& id_list);
  Sample extract_filtered() const;

  // Data
  vector<string> ID, famID;
  vector<bool> missing;
  vector<int> sex;
  vector<bool> keep;
};


class Legend {
public:
  // Constructors/destructors
  Legend();
  Legend(const Legend& obj);
  ~Legend();

  // Methods
  void clear();
  void push_back(vector<string> row);
  vector<string> row(int idx) const;
  int size() const;
  void filter(const vector<string>& id_list);
  Legend extract_filtered() const;
  void print() const;

  // Data
  vector<string> chr, ID;
  vector<int> bp;
  vector<double> cm;
  vector<string> A0, A1;
  vector<bool> extract;
};

class Metadata {
public:
  // Constructors/destructors
  Metadata(const DataFiles& _df, string _chr_id, int _window_size);
  Metadata(const Metadata& obj);
  ~Metadata();

  // Methods
  int num_snps() const;
  int num_samples() const;
  int num_haps() const;
  string get_chr_id() const;
  Metadata thin(int compression) const;
  void print_segments() const;

  // Data
  DataFiles df;                           // Data files
  Legend legend, legend_filter;           // Variant legend
  Sample sample, sample_filter;           // Sample information
  vector< vector<int> > partitions;       // Variant partitions
  vector< vector<IbdSeg> > ibd_segments;  // IBD segments
  vector< vector<int> > related_samples;  // Groups of individuals sharing at least one IBD segment
  vector< vector<int> > related_families; // Indices of groups of haplotypes sharing at least one IBD segment
  set<int> unrelated;                     // List of unrelated individuals
  set<int> related;                       // List of unrelated individuals
  Windows windows;

protected:
  string chr_id;
  int num_snps_;
  int num_samples_;
  int window_size;

  void load_sample();
  void load_legend();
  void apply_filters();
  void read_matrix(string filename, int num_columns, int skip, vector< vector<string> >& output) const;
  void load_partitions();
  void load_genetic_map();
  int load_ibd_segments();
  int find_variant(string _id) const;
  int find_variant(int _bp) const;
  int find_sample(string _id) const;
  void interpolate(const vector<double>& x, const vector<double>& y,
                   const vector<double>& x_new, vector<double>& y_new) const;
  void define_windows(int _window_size);
};


#endif
