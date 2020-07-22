#ifndef _KINSHIP_COMPUTER_CPP
#define _KINSHIP_COMPUTER_CPP

#include "kinship_computer.h"

#define DEBUG false
#define DELIMITER "--------------------------------------------------------------------------------"

KinshipComputer::KinshipComputer(const vector<Metadata>& _metadata) : metadata(_metadata){
}

KinshipComputer::KinshipComputer(const vector<Metadata>& _metadata, int compression,
                                 int _cluster_size_min, int _cluster_size_max, int _num_threads):
  metadata(_metadata) {

  covariates_available = false;

  // Load data
  genomes = Dataset(metadata, _num_threads, compression);
  num_chrs = genomes.num_chrs();
  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes.push_back(&(genomes.chromosomes[chr]));
  }

  num_haps = chromosomes[0]->get_num_haps();
  cluster_size_min = _cluster_size_min;
  cluster_size_max = _cluster_size_max;

  num_threads = _num_threads;
  //num_threads = 1;

  // K-means hierarchical clustering
  cluster(true);

}

KinshipComputer::KinshipComputer(const vector<Metadata>& _metadata, int compression,
                                 const Covariates & covariates_,
                                 int _cluster_size_min, int _cluster_size_max, int _num_threads) :
  metadata(_metadata) {

  //cout << "[DEBUG] in KinshipComputer::KinshipComputer()" << endl;

  // Load data
  genomes = Dataset(metadata, _num_threads, compression);
  num_chrs = genomes.num_chrs();

  for(int chr=0; chr<num_chrs; chr++) {
    chromosomes.push_back(&(genomes.chromosomes[chr]));
  }

  // Load covariates
  covariates_available = true;
  covariates = covariates_;

  num_chrs = metadata.size();
  num_haps = chromosomes[0]->get_num_haps();
  cluster_size_min = _cluster_size_min;
  cluster_size_max = _cluster_size_max;

  num_threads = _num_threads;
  //num_threads = 1;

  // Sanity checks
  assert(num_haps == 2 * covariates_.sample.size());

  // K-means hierarchical clustering
  cluster(false);
}

KinshipComputer::~KinshipComputer() {
}

double KinshipComputer::distance(int i, int j) const {
  return(distance_leaveout(i, j, -1));
}

double KinshipComputer::distance_hap_window(const int i1, const int i2, const int chr,
                                            const Windows& windows, const int w) const {
  // cout << "[DEBUG] in distance_hap_window(" << i1   << ", " << i2  << ", "  << chr << ", w=" << w << flush;
  //cout << ", j_start = " << windows.start[w] << ", j_end = " << windows.end[w] << ")" << endl;
  // Compute distance using genetic data around a window
  int result = 0;

  const int num_windows = windows.num_windows;

  if(num_windows==1) {
    result = chromosomes[chr]->distance(i1, i2);
  } else {
    if(w>0) {
      // Compute distance using windows on the left
      const int j_start = windows.start[0];
      const int j_end = windows.end[w-1];
      result += chromosomes[chr]->distance(i1, i2, j_start, j_end);
    }
    if(w<(num_windows-1)) {
      // Compute distance using windows on the right
      const int j_start = windows.start[w+1];
      const int j_end = windows.end.back();
      result += chromosomes[chr]->distance(i1, i2, j_start, j_end);
    }
  }
  return(result);
}

double KinshipComputer::distance_leaveout(int i, int j, int chr_out) const {
  int result = 0;
  if(covariates_available) {
    // Compute distance using covariates
    result = distance(covariates.Z[i/2], covariates.Z[j/2]);
  } else {
    // Compute distance using genetic data
    const int i1a = 2*(i/2);
    const int i1b = 2*(i/2)+1;
    const int i2a = 2*(j/2);
    const int i2b = 2*(j/2)+1;
    for(int chr=0; chr<num_chrs; chr++) {
      if(chr!=chr_out) {
        result += chromosomes[chr]->distance(i1a, i2a) + chromosomes[chr]->distance(i1b, i2b);
      }
    }
  }
  return(result);
}

double KinshipComputer::distance(int i, int j, int chr) const {
  int result = 0;
  if(covariates_available) {
    // Compute distance using covariates
    result = distance(covariates.Z[i/2], covariates.Z[j/2]);
  } else {
    result = chromosomes[chr]->distance(i, j);
  }
  return(result);
}

void KinshipComputer::distance(int i, int j, vector<int> & output) const {
  int i0 = 2*(i/2);
  int i1 = 2*(i/2)+1;
  int j0 = 2*(j/2);
  int j1 = 2*(j/2)+1;
  // Compute distances on each chromosome
  for(int chr=0; chr<num_chrs; chr++) {
    output[chr] = chromosomes[chr]->distance(i0, j0) + chromosomes[chr]->distance(i1, j1);
  }
  // Compute total distance
  int dist_total = std::accumulate(output.begin(), output.end(), 0);
  // Compute hold-one-out distances
  // If only one chromosome is present, do not hold out
  for(int chr=0; chr<num_chrs; chr++) {
    if(num_chrs>1) output[chr] = dist_total-output[chr];
    else output[chr] = dist_total;
  }
}

double KinshipComputer::distance(const vector2d & mu0, const vector2d & mu1) const {
  assert(mu0.size()==mu1.size());
  double dist = 0.0;
  for(int chr=0; chr<mu0.size(); chr++) {
    dist += distance(mu0[chr], mu1[chr]);
  }
  return(dist);
}

double KinshipComputer::distance(const vector<double> & mu0, const vector<double> & mu1) const {
  assert(mu0.size()==mu1.size());
  double dist = 0.0;
  for(int j=0; j<mu0.size(); j++) {
    dist += std::pow(mu0[j]-mu1[j], 2);
  }
  return(dist);
}

void KinshipComputer::findNeighbors(const int i, const int K, const int chr, const Windows& windows,
                                    ivector2d& output) const {
  //cout << "[DEBUG] in findNeighbors(" << i <<", " << K << ")" << endl;

  // Find IBD-sharing family for this sample
  const vector<int>& ibd_family = metadata[chr].related_families[i];

  // Compute distance leaving out one window
  const int num_windows = windows.num_windows;

  for(int w=0; w<num_windows; w++) {
    distmap D;

    int ci = assignments[chr][i];

    for(int l=0; l<clusters[chr][ci].size(); l++) {
      const int j = clusters[chr][ci][l];
      if(i==j) continue;

      // Check whether the individuals are in the same IBD-sharing family
      bool is_related = false;
      if(ibd_family.size()>1) {
        is_related = std::count(ibd_family.begin(), ibd_family.end(), j) > 0;
      }
      double dij;

      if(is_related) {
        // Prevent IBD-related samples to be assigned as references
        dij = 1e6;
      } else {
        dij = distance_hap_window(i, j, chr, windows, w);
      }

      //cout << "\t computed d(" << i << ", "<< j << ", " << w << "): " << dij << endl;

      // Only add new measurement if within top K
      if(D.size() >= K) {
        // Find largest distance
        int top_key = D.rbegin()->first;
        // Check whether the new distance should be in the top K
        if(dij < top_key) {
          // Insert new distance
          D.emplace(dij,j);
          // Remove largest distance
          D.erase(std::next(D.begin(),K));
        }
      } else {
        D.emplace(dij,j);
      }
    }

    // Copy list of neighbors to output
    int k = 0;
    for (auto it=D.begin(); it!=D.end(); ++it) {
      output[w][k] = it->second;
      k++;
    }

  }
}

double KinshipComputer::distance(int i, const vector2d & mu) const {
  return(distance(i, -1, mu));
}

double KinshipComputer::distance(int i, int chr_out, const vector2d & mu) const {
  double tot_dist = 0.0;
  for(int chr=0; chr<num_chrs; chr++) {
    if(chr==chr_out) continue;
    int num_snps = chromosomes[chr]->get_num_snps();
    for(int j=0; j<num_snps; j++) {
      double hj = (double)(chromosomes[chr]->H[i][j]);
      tot_dist += std::pow(hj-mu[chr][j],2);
    }
  }
  return(tot_dist);
}

void KinshipComputer::init_centroids(int c, int chr_out, const vector< vector<int> > & clust_chr,
                                     vector2d & mu0, vector2d & mu1) const {
  // Choose haplotypes for new centroids
  // Draw 100 pairs, take the pair with the greatest distance from each other.
  int clust_size = clust_chr[c].size();
  int dil_max = 0;
  int i_best = 0;
  int l_best = 0;
  int num_pairs = 100;
  for(int p=0; p<num_pairs; p++) {
    int i = clust_chr[c][putils::getRandom(clust_size)];
    int i0 = 2*(i/2);
    int i1 = 2*(i/2)+1;
    int l = i;
    while( (l==i0) || (l==i1) ) {
      l = clust_chr[c][putils::getRandom(clust_size)];
    }
    int dil = distance_leaveout(i, l, chr_out);
    if(dil>dil_max) {
      i_best = i;
      l_best = l;
      dil_max = dil;
    }
  }
  if(DEBUG) cout << endl << "[DEBUG] " << " Initial centroids: " << i_best << ", " << l_best << endl;

  // Compute centroids
  mu0.resize(num_chrs);
  mu1.resize(num_chrs);
  for(int chr=0; chr<num_chrs; chr++) {
    if(covariates_available) {
      int dim = covariates.K;
      mu0[chr].resize(dim, 0);
      mu1[chr].resize(dim, 0);
      for(int j=0; j<dim; j++) {
        mu0[chr][j] = (double) covariates.Z[i_best/2][j];
        mu1[chr][j] = (double) covariates.Z[l_best/2][j];
      }
    } else {
      int num_snps = chromosomes[chr]->get_num_snps();
      if(chr==chr_out) {
        mu0[chr].resize(num_snps, 0);
        mu1[chr].resize(num_snps, 0);
      } else {
        mu0[chr].resize(num_snps);
        mu1[chr].resize(num_snps);
        for(int j=0; j<num_snps; j++) {
          mu0[chr][j] = (double)(chromosomes[chr]->H[i_best][j]);
          mu1[chr][j] = (double)(chromosomes[chr]->H[l_best][j]);
        }
      }
    }
  }
}

void print_multimap(const std::multimap<int,int> & pairs) {
  cout << "Key\tValue\n";
  // use const_iterator to walk through elements of pairs
  for ( auto iter = pairs.begin(); iter != pairs.end(); ++iter )
    cout << iter->first << '\t' << iter->second << '\n';
}

void KinshipComputer::update_assignments(int c, int chr_out, const vector< vector<int> > & clust_chr,
                                         const vector2d & mu0, const vector2d & mu1,
                                         vector<int> & new_assignments) const {

  std::multimap<double,int> distances0, distances1;
  int csize0 = 0, csize1 = 0;
  // Update assignments based on centroids (old code)
  for(int idx=0; idx<clust_chr[c].size(); idx++) {
    int i = clust_chr[c][idx];
    double dist_0, dist_1;
    if(covariates_available) {
      dist_0 = distance(covariates.Z[i/2], mu0[0]);
      dist_1 = distance(covariates.Z[i/2], mu1[0]);
    } else {
      int i0 = 2*(i/2);
      int i1 = 2*(i/2)+1;
      dist_0 = distance(i0, chr_out, mu0) + distance(i1, chr_out, mu0);
      dist_1 = distance(i0, chr_out, mu1) + distance(i1, chr_out, mu1);
    }
    if(dist_0<dist_1) {
      new_assignments[idx] = 0;
      csize0++;
      distances0.emplace(dist_0,idx);
    }
    else {
      new_assignments[idx] = 1;
      csize1++;
      distances1.emplace(dist_1,idx);
    }
  }

  // Check cluster sizes
  if(DEBUG) cout<<"[DEBUG] assignments for cluster "<<c<<" ("<<csize0<<", "<<csize1<<")"<<endl;

  // Make sure that cluster size >= cluster_size_min (this is optional, but recommended)
  if(csize0 < cluster_size_min) {
    int n_missing = cluster_size_min - csize0;
    // cout << endl << "Sub-cluster 0 is too small, missing " << n_missing << endl;
    for (auto it=distances1.rbegin(); it!=std::next(distances1.rbegin(),n_missing); ++it) {
      // std::cout << " (" << (*it).first << ", " << (*it).second << ") ";
      new_assignments[(*it).second] = 0;
    }
    // cout << endl;
  }
  if(csize1 < cluster_size_min) {
    int n_missing = cluster_size_min - csize1;
    // cout << endl << "Sub-cluster 1 is too small, missing" << n_missing << endl;
    for (auto it=distances0.rbegin(); it!=std::next(distances0.rbegin(),n_missing); ++it) {
      // std::cout << " (" << (*it).first << ", " << (*it).second << ")";
      new_assignments[(*it).second] = 1;
    }
    // cout << endl;
  }

  //exit(0);

}

void KinshipComputer::update_centroids(int c, int chr_out, const vector< vector<int> > & clust_chr,
                                       const vector<int> & assign_chr, vector2d & mu0, vector2d & mu1) const {
  // Compute centroids based on assignments
  for(int chr=0; chr<num_chrs; chr++) {
    // Reset centroids
    std::fill(mu0[chr].begin(),mu0[chr].end(),0.0);
    std::fill(mu1[chr].begin(),mu1[chr].end(),0.0);
    // Skip left-out chromosome
    if(chr==chr_out) continue;
    // Update means
    int num_haps_0 = 0, num_haps_1 = 0;
    for(int idx=0; idx<clust_chr[c].size(); idx++) {
      int i = clust_chr[c][idx];
      int ci = assign_chr[idx];
      if(ci==0) num_haps_0++;
      if(ci==1) num_haps_1++;
      if(covariates_available) {
        for(int j=0; j<mu0[chr].size(); j++) {
          if(ci==0) mu0[chr][j] += (double) covariates.Z[i/2][j];
          if(ci==1) mu1[chr][j] += (double) covariates.Z[i/2][j];
        }

      } else {
        for(int j=0; j<mu0[chr].size(); j++) {
          if(ci==0) mu0[chr][j] += (double)(chromosomes[chr]->H[i][j]);
          if(ci==1) mu1[chr][j] += (double)(chromosomes[chr]->H[i][j]);
        }
      }
    }
    // Normalize centroids
    for(int j=0; j<mu0[chr].size(); j++) {
      mu0[chr][j] /= num_haps_0;
      mu1[chr][j] /= num_haps_1;
    }
  }
}

void KinshipComputer::cluster(bool separate_chrs) {
  //cout << "[DEBUG] in KinshipComputer::cluster" << endl;
  if(separate_chrs) {

    // Do not use more threads than chromosomes to load the data
    int num_threads_chr = std::min(num_threads, num_chrs);

    // Assign chromosomes to workers
    vector< vector<int> > worker_assignments(num_threads_chr);
    int nchrs_remaining = num_chrs;
    int nthreads_remaining = num_threads_chr;
    int share = 1 + (nchrs_remaining-1) / nthreads_remaining;
    int wid = 0;
    for(int chr=0; chr<num_chrs; chr++) {
      worker_assignments[wid].push_back(chr);
      nchrs_remaining--;
      if((worker_assignments[wid].size() >= share) && (nthreads_remaining>1)) {
        wid++;
        nthreads_remaining--;
        share = 1 + (nchrs_remaining-1) / nthreads_remaining;
      }
    }

    // // DEBUG: check whether this is correct
    // for(int w=0; w<num_threads_chr; w++) {
    //   cerr << "[DEBUG] w " << w << " : ";
    //   for(int i=0; i<worker_assignments[w].size(); i++) {
    //     cerr << worker_assignments[w][i] + 1 << " ";
    //   }
    //   cerr << endl;
    // }

    cout << "Solving " << num_chrs << " bifurcating K-means problems using " << num_threads_chr;
    cout << " threads." << endl;

    // Bifurcating K-means clustering (leaving out one chromosome at a time)
    assignments.resize(num_chrs);
    clusters.resize(num_chrs);
    vector<int> progress(num_threads_chr,0);
    if(num_chrs>1) {
      // If multiple chrs
      if(num_threads_chr>1) {
        vector<boost::thread> workers;
        for(int w=0; w<num_threads_chr; w++ ) {
          // Divide task
          int chr_start = worker_assignments[w].front();
          int chr_end = worker_assignments[w].back();
          // Create worker
          workers.push_back(boost::thread(&KinshipComputer::cluster_worker, this, chr_start, chr_end,
                                          w, boost::ref(progress)));
        }
        // Launch workers
        for(int w=0; w<num_threads_chr; w++ ) {
          workers[w].join();
        }
      } else {
        cluster_worker(0, num_chrs-1, 0, progress);
      }
    } else {
      // If only one chromosome
      bifurcating_kmeans(-1, assignments[0], clusters[0], true);
    }
    cout << endl;
  } else {
    // Run on one chromosome only, copy results to other chromosomes
    assignments.resize(num_chrs);
    clusters.resize(num_chrs);
    bifurcating_kmeans(-1, assignments[0], clusters[0], true);
    for(int chr=1; chr<num_chrs; chr++) {
      assignments[chr] = assignments[0];
      clusters[chr] = clusters[0];
    }
  }
}

void KinshipComputer::cluster_worker(int chr_min, int chr_max, int wid, vector<int> &progress) {
  int n_steps = num_chrs;

  //cout << "[DEBUG] in cluster_worker(" << chr_min << "," << chr_max << "," << wid << ")" << endl;

  // Initialize progress bar (worker 0)
  if(wid==0) {
    cout << "|" << flush;
    for(int s=0; s<n_steps; s++) cout << ".";
    cout << "|" << endl;
    cout << "|" << flush;
  }

  bool verbose = wid==0;

  for(int chr=chr_min; chr<=chr_max; chr++) {
    bifurcating_kmeans(chr, assignments[chr], clusters[chr], verbose);
    // Keep track of progress (all workers)
    progress[wid]++;
    cout << "=" << flush;
  }
  // Finalize progress bar (last worker to finish)
  if( std::accumulate(progress.begin(),progress.end(),0) == n_steps) {
    cout << "|" << endl;
  }
}

void KinshipComputer::bifurcating_kmeans(int chr_out, vector<int> & assign_chr,
                                         vector< vector<int> > & clust_chr, bool verbose) const {
  assert(num_haps%2==0);

  // Initialize clustering
  assign_chr.resize(num_haps);
  std::fill(assign_chr.begin(),assign_chr.end(),0);
  vector<int> clust_chr_tmp(num_haps);
  std::iota(clust_chr_tmp.begin(), clust_chr_tmp.end(), 0);
  clust_chr.clear();
  clust_chr.push_back(clust_chr_tmp);
  vector<int> unfinished_clusters;
  unfinished_clusters.push_back(0);

  if(verbose) {
    cerr << endl << "Bifurcating K-means";
    if(chr_out!=-1) cerr << " (leaving out chromosome " << chr_out+1 << ")";
    cerr << endl << "Smallest allowed cluster size: ";
    cerr << cluster_size_min << endl;
    cerr << setfill(' ') << setw(5) << "step" << "\t";
    cerr << setfill(' ') << setw(8) << "cluster" << "\t";
    cerr << setfill(' ') << setw(8) << "size" << "\t";
    cerr << setfill(' ') << setw(8) << "left" << "\t";
    cerr << setfill(' ') << setw(8) << "right" << "\t";
    cerr << setfill(' ') << setw(8) << "accepted";
    cerr << endl;
  }

  // Bifurcating k-means
  int step = 0;
  int max_steps = 10000;
  while((unfinished_clusters.size()>0) && (step++<max_steps)) {
    int c = unfinished_clusters[0];

    // If the cluster is small enough, no further splitting
    if(clust_chr[c].size() <= cluster_size_max) {
      auto it = find(unfinished_clusters.begin(), unfinished_clusters.end(), c);
      unfinished_clusters.erase(it);
      continue;
    }

    if(verbose) {
      cerr << setfill(' ') << setw(5) << step << "\t" << flush;
      cerr << setfill(' ') << setw(8) << c << "\t" << flush;
      cerr << setfill(' ') << setw(8) << clust_chr[c].size() << "\t" << flush;

      tuple<bool,int,int> split = kmeans(c, chr_out, clust_chr, assign_chr, unfinished_clusters);

      cerr << setfill(' ') << setw(8) << std::get<1>(split) << "\t" << flush;
      cerr << setfill(' ') << setw(8) << std::get<2>(split) << "\t" << flush;
      if(std::get<0>(split)) {
        cerr << setfill(' ') << setw(8) << "yes" << "\t" << flush;
      } else {
        cerr << setfill(' ') << setw(8) << "no" << "\t" << flush;
      }
      cerr << endl;

    } else {
      kmeans(c, chr_out, clust_chr, assign_chr, unfinished_clusters);
    }
  }

  if(verbose) {
    if(step>=max_steps) {
      cerr << "Warning: not completed after " << step-1 << " steps." << endl;
    } else {
      cerr << "Bifurcating K-means completed after " << step-1 << " steps." << endl;
      cerr << "Number of clusters: " << clust_chr.size() << "." << endl << endl;
    }
  }

  // // DEBUG
  // cout << "Printing assignments" << endl;
  // for(auto x : assign_chr) {
  //   cout << " " << x;
  // }
  // cout << endl;
 
  // // Return list of cluster assignments in terms of single haplotypes
  // assign_chr_.resize(num_haps);
  // for(int i=0; i<num_haps/2; i++) {
  //   assign_chr_[2*i] = assign_chr[i];
  //   assign_chr_[2*i+1] = assign_chr[i];
  // }

  // // Return list of clusters in terms of single haplotypes
  // clust_chr_.clear();
  // for(int c=0; c<clust_chr.size(); c++) {
  //   vector<int> tmp(2*clust_chr[c].size());
  //   for(int i=0; i<clust_chr[c].size(); i++) {
  //     tmp[2*i] = 2*clust_chr[c][i];
  //     tmp[2*i+1] = 2*clust_chr[c][i]+1;
  //   }
  //   clust_chr_.push_back(tmp);
  // }
}

void KinshipComputer::kmeans_core(int c, int chr_out, vector< vector<int> > & clust_chr,
                                  vector<int> & new_assign_chr) const {
  // Initialize centroids
  vector2d mu0, mu1;
  init_centroids(c, chr_out, clust_chr, mu0, mu1);
  vector2d mu0_old = mu0;
  vector2d mu1_old = mu1;

  // Update assignments and centroids until convergence
  double delta_mu0=1, delta_mu1=1;
  double mu_eps = 1e-3;
  int it_max = 50;
  for(int it=0; it<it_max; it++) {
    // Update assignments and centroids
    update_assignments(c, chr_out, clust_chr, mu0, mu1, new_assign_chr);
    update_centroids(c, chr_out, clust_chr, new_assign_chr, mu0, mu1);
    // Check convergence
    delta_mu0 = distance(mu0, mu0_old);
    delta_mu1 = distance(mu1, mu1_old);
    mu0_old = mu0;
    mu1_old = mu1;
    if((delta_mu0<mu_eps)&&(delta_mu1<mu_eps)) {
      break;
    }
  }
}


tuple<bool,int,int> KinshipComputer::kmeans(int c, int chr_out, vector< vector<int> > & clust_chr,
                                            vector<int> & assign_chr, vector<int> & unfinished_clusters) const {

  int clust_size = clust_chr[c].size();

  vector<int> new_assign_chr(clust_size, 0);
  kmeans_core(c, chr_out, clust_chr, new_assign_chr);

  // Compute sizes of candidate clusters to decide whether to accept them
  int clust_size_0 = 0, clust_size_1 = 0;
  for(int idx=0; idx<clust_size; idx++) {
    if(new_assign_chr[idx]==0) {
      clust_size_0++;
    } else {
      clust_size_1++;
    }
  }

  // Define new clusters and assignments, if accepted
  if(std::min(clust_size_0,clust_size_1) >= cluster_size_min) {
    int num_clust = clust_chr.size();

    // Define new assignments and clusters
    vector<int> new_cluster_0, new_cluster_1;
    for(int idx=0; idx<clust_size; idx++) {
      int i = clust_chr[c][idx];
      if(new_assign_chr[idx]==0) {
        new_assign_chr[idx] = c;
        new_cluster_0.push_back(i);
      } else {
        new_assign_chr[idx] = num_clust;
        new_cluster_1.push_back(i);
      }
      assign_chr[i] = new_assign_chr[idx];
    }

    // Remove old cluster
    clust_chr.erase(clust_chr.begin()+c);

    // Insert new clusters
    clust_chr.insert(clust_chr.begin()+c, new_cluster_0);
    clust_chr.push_back(new_cluster_1);

    // Update list of unfinished clusters
    if(clust_size_0<=cluster_size_min) {
      auto it = find (unfinished_clusters.begin(), unfinished_clusters.end(), c);
      unfinished_clusters.erase(it);
    }
    if(clust_size_1>=cluster_size_min) {
      unfinished_clusters.push_back(num_clust);
    }
    return(std::make_tuple(true,clust_size_0,clust_size_1));
  } else {
    return(std::make_tuple(false,clust_size_0,clust_size_1));
  }
}

void KinshipComputer::writeClusters(const vector<string> & out_file_names) const {
  for(int chr=0; chr<num_chrs; chr++) {
    string out_file_name = out_file_names[chr] + "_clust.txt";
    ofstream outfile(out_file_name.c_str());
    if (!outfile.is_open()){
      cout << "Problem creating the output file: " << out_file_name;
      cout <<"Either the directory does not exist or you do not have write permissions." << endl;
    }
    for(int c=0; c<clusters[chr].size(); c++) {
      for(int i=0; i<clusters[chr][c].size(); i++) {
        outfile << clusters[chr][c][i];
        if(i+1<clusters[chr][c].size()) outfile <<" ";
      }
      outfile << endl;
    }
    outfile.close();

    // Write list of samples
    plink::writeSAMPLE(out_file_names[chr], metadata[chr]);

    cout << "Kinship clusters ";
    if(num_chrs>1) cout << "(leaving out chromosome " << chr+1 << ") ";
    cout << "written to:" << endl;
    cout << "  " << out_file_name << endl;
    cout << "  " << out_file_names[chr] + ".sample" << endl;
  }
  cout << endl;
}

void KinshipComputer::assign_references_worker(int i_min, int i_max, int K, int chr, ivector3d & output,
                                               const Windows& windows, int wid, vector<int> & progress) const {

  int step = 0;
  const int n_steps = 100;
  const int progress_period = std::max(1, (i_max-i_min+1)/n_steps);
  const int progress_total = std::accumulate(progress.begin(), progress.end(), 0);

  // Initialize progress bar (worker 0)
  if(wid==0) {
    cout << "|" << flush;
    for(int s=0; s<n_steps; s++) cout << ".";
    cout << "|" << endl;
    cout << "|" << flush;
  }

  for(int i=i_min; i<i_max; i++) {
    findNeighbors(i, K, chr, windows, output[i]);
    // Keep track of progress (all workers)
    progress[wid]--;
    // Update progress bar (worker 0)
    if(wid==0) {
      if(i % progress_period == 0) {
        const int progress_remaining = std::accumulate(progress.begin(),progress.end(),0);
        int new_step = (n_steps * (progress_total-progress_remaining)) / progress_total;
        for(int s=step; s<new_step; s++) {
          cout << "=" << flush;
          step = new_step;
        }
      }
    }
  }

  // Finalize progress bar (worker 0)
  if(wid==0) {
    for(int s=step; s<n_steps; s++) {
      cout << "=" << flush;
    }
    cout << "|" << endl;
  }
}


ivector3d KinshipComputer::combine_references_families(const int chr, const ivector3d& _ref) const {
  ivector3d ref = _ref;
  const int num_windows = ref[0].size();
  for(int w=0; w<num_windows; w++) {
    for(int i=0; i<num_haps; i++) {
      std::sort(ref[i][w].begin(), ref[i][w].end());
    }
  }

  for(int w=0; w<num_windows; w++) {
    const int K = ref[0][w].size();
    vector<bool> already_processed(num_haps, false);
    for(int i1=0; i1<num_haps; i1++) {
      // Check whether this haplotype has already been processed
      if(already_processed[i1]) continue;

      // If not, add it and its family to the list
      const vector<int>& ibd_family = metadata[chr].related_families[i1];
      already_processed[i1] = true;
      for(int i2 : ibd_family) already_processed[i2] = true;

      // Skip family if it only contains one haplotype
      const int fam_size = ibd_family.size();
      if(fam_size<=1) continue;

      // Find intersection of all references within family
      vector<int> v_intersection = ref[i1][w];
      for(int i2 : ibd_family) {
        if(i2==i1) continue;
        vector<int> v_tmp;
        std::set_intersection(v_intersection.begin(), v_intersection.end(),
                              ref[i2][w].begin(), ref[i2][w].end(),
                              back_inserter(v_tmp));
        v_intersection = v_tmp;
      }

      // cout << "Family of haplotype " << i1 << ":";
      // for(int i : ibd_family) cout << " " << i;
      // cout << endl;
      // cout << "Intersection has size " << v_intersection.size() << endl;

      // Add top references for each family element until full
      int i = 0;
      int k = 0;
      while(v_intersection.size()<K) {
        int id = _ref[ibd_family[i]][w][k];
        if(!std::binary_search(v_intersection.begin(), v_intersection.end(), id)) {
          v_intersection.push_back(id);
          sort(v_intersection.begin(), v_intersection.end());
        }
        if(i+1<fam_size) {
          i++;
        } else {
          i=0;
          k++;
        }
      }

      // cout << "New reference set:" << endl;
      // for(int i : v_intersection) cout << " " << i;
      // cout << endl;

      // Store new set of references
      for(int i : ibd_family) {
        ref[i][w] = v_intersection;
      }
    }
  }

  return(ref);
}

ivector3d KinshipComputer::assign_references(const int K, const int chr, const Windows& windows) const {

  // Find number of windows
  int num_windows = windows.num_windows;

  // Initialize references
  ivector3d ref(num_haps, ivector2d(num_windows, vector<int>(K,-1)));

  if(num_windows>1) {
    cout<<"Assigning references within "<<num_windows<<" windows using "<<num_threads<<" threads: "<<endl;
  } else {
    cout<<"Assigning references for the whole chromosome using "<<num_threads<<" threads: "<<endl;
  }

  vector<int> progress(num_threads);
  for(int wid=0; wid<num_threads; wid++ ) {
    int i_start = wid*num_haps/num_threads;
    int i_end = (wid+1)*num_haps/num_threads;
    progress[wid] = i_end - i_start + 1;
  }

  if(num_threads>1) {
    // Compute distances for blocks of individuals in parallel
    vector<boost::thread> workers;
    for(int wid=0; wid<num_threads; wid++ ) {
      // Divide task
      int i_start = wid*num_haps/num_threads;
      int i_end = (wid+1)*num_haps/num_threads;
      if(wid==(num_threads-1)) i_end = num_haps;
      // Create worker
      workers.push_back(boost::thread(&KinshipComputer::assign_references_worker, this, i_start, i_end, K, chr,
                                      boost::ref(ref), boost::ref(windows), wid, boost::ref(progress)));
    }
    // Launch workers
    for(int wid=0; wid<num_threads; wid++ ) {
      workers[wid].join();
    }
  } else {
    progress[0] = num_haps;
    assign_references_worker(0, num_haps, K, chr, ref, windows, 0, progress);
  }
  cout << endl;

  // Combine references within IBD-sharing families
  ref = combine_references_families(chr, ref);

  // Sanity check: make sure references are constant within families
  for(int w=0; w<num_windows; w++) {
    for(int i1=0; i1<num_haps; i1++) {
      const vector<int>& ibd_family = metadata[chr].related_families[i1];
      if(ibd_family.size()>1) {
        for(int i2 : ibd_family) {
          if(ref[i1][w] != ref[i2][w]) {
            cerr << "Error: references for haplotypes " << i1 << " and " << i2 << " do not match!" << endl;
            cout << "Family:";
            for(int i : ibd_family) cout << " " << i;
            cout << endl;
            for(int k : ref[i1][w]) cout << " " << k;
            cout << endl;
            for(int k : ref[i2][w]) cout << " " << k;
            cout << endl;
            exit(-1);
          }
        }
      }
    }
  }

  // Sanity check: make sure all references have been correctly assigned
  for(int w=0; w<num_windows; w++) {
    for(int i=0; i<num_haps; i++) {
      const vector<int>& ibd_family = metadata[chr].related_families[i];
      for(int k=0; k<K; k++) {
        int j = ref[i][w][k];
        if(j < 0) {
          cerr << "Error: could not assign all references for lack of sufficient samples. ";
          cerr << "Try decreasing K." << endl;
          exit(-1);
        }
        bool is_related = false;
        if(ibd_family.size()>1) {
          is_related = std::count(ibd_family.begin(), ibd_family.end(), j) > 0;
        }
        if(is_related)  {
          cerr << "Error: could not assign all references for lack of sufficient unrelated samples. ";
          cerr << "Try decreasing K." << endl;
          exit(-1);
        }
      }
    }
  }

  return(ref);
}


#endif
