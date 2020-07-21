#ifndef _KINSHIP_COMPUTER_H
#define _KINSHIP_COMPUTER_H

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include "dataset.h"
#include "plink_interface.h"
#include "metadata.h"
#include "covariates.h"
#include "windows.h"

typedef std::multimap<int,int> distmap;
typedef vector< vector<double> > vector2d;
typedef vector< vector<int> > ivector2d;
typedef vector<ivector2d> ivector3d;
typedef vector<ivector3d> ivector4d;

class KinshipComputer {
public:
  // Constructors and destructors
  KinshipComputer(const vector<Metadata>& _metadata);
  KinshipComputer(const vector<Metadata>& _metadata, int compression,
                  int _cluster_size_min, int _cluster_size_max, int _num_threads);
  KinshipComputer(const vector<Metadata>& _metadata, int compression, const Covariates & covariates_,
                  int _cluster_size_min, int _cluster_size_max, int _num_threads);
  ~KinshipComputer();

  // Methods
  void findNeighbors(const int i, const int K, const int chr, const Windows& windows, ivector2d & output) const;
  void writeClusters(const vector<string> & out_file_names) const;
  ivector3d assign_references(const int K, const int chr, const Windows& windows) const;


private:
  // Data
  int num_haps;
  int num_chrs;
  int num_threads;
  Dataset genomes;
  vector<const Haplotypes*> chromosomes;
  bool covariates_available;
  Covariates covariates;
  const vector<Metadata>& metadata;

  // Clustering
  vector< vector<int> > assignments;
  vector< vector< vector<int> > > clusters;
  int cluster_size_min, cluster_size_max;

  // Methods for reference assignment
  void assign_references_worker(int i_min, int i_max, int K, int chr, ivector3d & output, 
                                const Windows& windows, int wid, vector<int> & progress) const;
  ivector3d combine_references_families(const int chr, const ivector3d& _ref) const;

  // Methods for distances
  double distance(int i, int j) const;
  double distance_leaveout(int i, int j, int chr_out) const;
  double distance_hap_window(const int i1, const int i2, const int chr, const Windows& windows, const int w) const;
  double distance(int i, int j, int chr) const;
  void distance(int i, int j, vector<int> & output) const;
  double distance(int i, const vector2d & mu) const;
  double distance(int i, int chr_out, const vector2d & mu) const;
  double distance(const vector2d & mu0, const vector2d & mu1) const;
  double distance(const vector<double> & mu0, const vector<double> & mu1) const;

  // Methods for K-means
  void cluster(bool separate_chrs);
  void cluster_worker(int chr_min, int chr_max, int wid, vector<int> &progress);
  void bifurcating_kmeans(int chr_out, vector<int> & assign_chr, vector< vector<int> > & clust_chr,
                          bool verbose) const;
  tuple<bool,int,int> kmeans(int cluster, int chr_out,
                             vector< vector<int> > & clust_chr, vector<int> & assign_chr,
                             vector<int> & unfinished_clusters) const;

  void kmeans_core(int c, int chr_out, vector< vector<int> > & clust_chr, vector<int> & new_assign_chr) const;
  void init_centroids(int c, int chr_out, const vector< vector<int> > & clust_chr,
                        vector2d & mu0, vector2d & mu1) const;
  void update_assignments(int c, int chr_out, const vector< vector<int> > & clust_chr,
                            const vector2d & mu0, const vector2d & mu1,
                            vector<int> & new_assign_chr) const;
  void update_centroids(int c, int chr_out,
                          const vector< vector<int> > & clust_chr, const vector<int> & assign_chr,
                          vector2d & mu0, vector2d & mu1) const;

};

#endif
