#ifndef _KINSHIP_COMPUTER_H
#define _KINSHIP_COMPUTER_H

#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include "dataset.h"
#include "plink_interface.h"
#include "metadata.h"
#include "covariates.h"
#include "windows.h"

typedef std::multimap<unsigned int, unsigned int> distmap;
typedef vector< vector<double> > vector2d;
typedef vector< vector<unsigned int> > ivector2d;
typedef vector<ivector2d> ivector3d;
typedef vector<ivector3d> ivector4d;

class KinshipComputer {
public:
  // Constructors and destructors
  KinshipComputer(const vector<Metadata>& _metadata);
  KinshipComputer(const vector<Metadata>& _metadata, unsigned int compression,
                  unsigned int _cluster_size_min, unsigned int _cluster_size_max, unsigned int _num_threads);
  KinshipComputer(const vector<Metadata>& _metadata, unsigned int compression, const Covariates & covariates_,
                  unsigned int _cluster_size_min, unsigned int _cluster_size_max, unsigned int _num_threads);

  ~KinshipComputer();

  // Methods
  void findNeighbors(const unsigned int i, const unsigned int K, const unsigned int chr, const Windows& windows, ivector2d & output) const;
  void writeClusters(const vector<string> & out_file_names) const;
  ivector3d assign_references(const unsigned int K, const unsigned int chr, const Windows& windows) const;


private:
  // Data
  unsigned int num_haps;
  unsigned int num_chrs;
  unsigned int num_threads;
  Dataset genomes;
  vector<const Haplotypes*> chromosomes;
  bool covariates_available;
  Covariates covariates;
  const vector<Metadata>& metadata;

  // Clustering
  vector< vector<unsigned int> > assignments;
  vector< vector< vector<unsigned int> > > clusters;
  unsigned int cluster_size_min, cluster_size_max;

  // Methods for reference assignment
  void assign_references_worker(unsigned int i_min, unsigned int i_max, unsigned int K, unsigned int chr, ivector3d & output, 
                                const Windows& windows, unsigned int wid, vector<unsigned int> & progress) const;
  ivector3d combine_references_families(const unsigned int chr, const ivector3d& _ref) const;

  // Methods for distances
  double distance(unsigned int i, unsigned int j) const;
  double distance_leaveout(unsigned int i, unsigned int j, int chr_out) const;
  double distance_hap_window(const unsigned int i1, const unsigned int i2, const unsigned int chr, const Windows& windows, const unsigned int w) const;
  double distance(unsigned int i, unsigned int j, unsigned int chr) const;
  void distance(unsigned int i, unsigned int j, vector<unsigned int> & output) const;
  double distance(unsigned int i, const vector2d & mu) const;
  double distance(unsigned int i, int chr_out, const vector2d & mu) const;
  double distance(const vector2d & mu0, const vector2d & mu1) const;
  double distance(const vector<double> & mu0, const vector<double> & mu1) const;

  // Methods for K-means
  void cluster(bool separate_chrs);
  void cluster_worker(unsigned int chr_min, unsigned int chr_max, unsigned int wid, vector<unsigned int> &progress);
  void bifurcating_kmeans(int chr_out, vector<unsigned int> & assign_chr, vector< vector<unsigned int> > & clust_chr,
                          bool verbose) const;
  tuple<bool,int,int> kmeans(unsigned int cluster, int chr_out,
                             vector< vector<unsigned int> > & clust_chr, vector<unsigned int> & assign_chr,
                             vector<unsigned int> & unfinished_clusters) const;

  void kmeans_core(unsigned int c, int chr_out, vector< vector<unsigned int> > & clust_chr, vector<unsigned int> & new_assign_chr) const;
  void init_centroids(unsigned int c, int chr_out, const vector< vector<unsigned int> > & clust_chr,
                        vector2d & mu0, vector2d & mu1) const;
  void update_assignments(unsigned int c, int chr_out, const vector< vector<unsigned int> > & clust_chr,
                            const vector2d & mu0, const vector2d & mu1,
                            vector<unsigned int> & new_assign_chr) const;
  void update_centroids(unsigned int c, int chr_out,
                          const vector< vector<unsigned int> > & clust_chr, const vector<unsigned int> & assign_chr,
                          vector2d & mu0, vector2d & mu1) const;

};

#endif
