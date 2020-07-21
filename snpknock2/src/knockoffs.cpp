#ifndef _KNOCKOFFS_CPP
#define _KNOCKOFFS_CPP

#define DEBUG 0

#include "knockoffs.h"

//////////////////////////
// Basic implementation //
//////////////////////////

Knockoffs::Knockoffs(const Metadata& _metadata, int _K, bool _load_data, int nthreads,
                     int _debug, int _seed, const string& _logfile) :
  Haplotypes(_metadata), K(_K), seed(_seed), debug(_debug), logfile(_logfile) {

  // Genomic windows
  num_windows = metadata.windows.num_windows;

  // Initialize local reference set
  references_local.resize(num_haps, ivector2d(num_windows, vector<int>(K,0)));
  alpha.resize(num_haps, vector< vector<double> > (num_windows, vector<double>(K,1.0/K)));

  // Initialize global reference set
  references_global.resize(num_haps, vector<int>(K,0));
  alpha_global.resize(num_haps, vector<double>(K,1.0/K));

  // Reset knockoff generator
  reset_knockoffs();

  // Load data from files
  if(_load_data) load(nthreads);
}

Knockoffs::Knockoffs(const Metadata& _metadata, int _K, int nthreads, int _debug, int _seed,
                     const string& _logfile) :
  Knockoffs(_metadata, _K, true, nthreads, _debug, _seed, _logfile) {}

void Knockoffs::load(const int nthreads) {
  bool verbose = true;

  // Initialize HMM
  genetic_dist.resize(metadata.num_snps(),0);
  phys_dist.resize(metadata.num_snps(),0);
  for(int j=1; j<metadata.num_snps(); j++) {
    double cm0 = metadata.legend_filter.cm[j-1];
    double cm1 = metadata.legend_filter.cm[j];
    if((cm0>0)&&(cm1>0)) {
      genetic_dist[j] = std::max(1e-9, cm1 - cm0);
    } else {
      genetic_dist[j] = 1e-3;
    }
    double bp0 = metadata.legend_filter.bp[j-1];
    double bp1 = metadata.legend_filter.bp[j];
    if((bp0>0)&&(bp1>0)) {
      phys_dist[j] = bp1 - bp0;
    } else {
      phys_dist[j] = 1e3;
    }
  }
  init_hmm(K, 1.0, 1.0e-4);

  // Load SNP partitions
  partitions = metadata.partitions;

  // Initialize group assignments
  set_partition(0);

  // Load haplotypes
  load_data(verbose, nthreads);
}

void Knockoffs::reset_knockoffs() {
  // Initialize HMM ancestry weights to default values
  alpha.clear();
  alpha.resize(num_haps, vector< vector<double> > (num_windows, vector<double>(K,1.0/K)));
  alpha_global.clear();
  alpha_global.resize(num_haps, vector<double>(K,1.0/K));

  // Initialize output
  vector<int> hk_placeholder(num_snps,0);
  chaplotype chk_placeholder(hk_placeholder);
  Hk = vector <chaplotype> (num_haps, chk_placeholder);

  // // DEBUG: store Z
  // Z.clear();
  // Z.resize(num_haps, vector<int>(num_snps, -1));
  // Zk.clear();
  // Zk.resize(num_haps, vector<int>(num_snps, -1));
}

void Knockoffs::set_partition(const int r) {
  assert(r<partitions.size());
  set_groups(partitions[r]);
}

void Knockoffs::set_groups(const vector<int> & _groups) {
  //cout << "Setting partition:" << endl;
  groups.clear();
  groups.resize(num_snps);
  elements.clear();
  num_groups_ = _groups.back()+1;
  elements.resize(num_groups_);
  for(int j=0; j<num_snps; j++) {
    groups[j] = _groups[j];
    elements[groups[j]].push_back(j);
  }
}

void Knockoffs::set_reference(const int i, const ivector2d& input) {
  assert(input.size() == num_windows);
  for(int w=0; w<num_windows; w++) {
    assert(input[w].size() == K);
    for(int k=0; k<K; k++) {
      references_local[i][w][k] = input[w][k];
    }
    std::sort(references_local[i][w].begin(), references_local[i][w].end());
  }
}

void Knockoffs::set_references(const ivector3d& input) {
  assert(input.size()==num_haps);
  for(int i=0; i<num_haps; i++) {
    set_reference(i, input[i]);
  }
}

void Knockoffs::get_reference(const int i, ivector2d & output) const {
  for(int w=0; w<num_windows; w++) {
    for(int k=0; k<K; k++) {
      output[w][k] = references_local[i][w][k];
    }
  }
}

void Knockoffs::set_reference_global(const int i, const vector<int>& input) {
  assert(input.size() == K);
  for(int k=0; k<K; k++) {
    references_global[i][k] = input[k];
  }
  std::sort(references_global[i].begin(), references_global[i].end());
}

void Knockoffs::set_references_global(const ivector2d& input) {
  assert(input.size()==num_haps);
  for(int i=0; i<num_haps; i++) {
    set_reference_global(i, input[i]);
  }
}

void Knockoffs::get_reference_global(const int i, vector<int>& output) const {
  for(int k=0; k<K; k++) {
    output[k] = references_global[i][k];
  }
}

void Knockoffs::set_ancestry(const int i, const int w, const vector<double> & input) {
  for(int k=0; k<K; k++) {
    alpha[i][w][k] = input[k];
  }
}

void Knockoffs::get_ancestry(const int i, const int w, vector<double> & output) const {
  for(int k=0; k<K; k++) {
    output[k] = alpha[i][w][k];
  }
}

void Knockoffs::set_ancestry_global(const int i, const vector<double> & input) {
  for(int k=0; k<K; k++) {
    alpha_global[i][k] = input[k];
  }
}

void Knockoffs::get_ancestry_global(const int i, vector<double> & output) const {
  for(int k=0; k<K; k++) {
    output[k] = alpha_global[i][k];
  }
}

void Knockoffs::write(const string& basename, const string& format, bool include_genotypes) const {
  if (format=="BED") {
    if(include_genotypes) {
      plink::write_binary(basename, metadata, H, Hk, true);
    } else {
      plink::write_binary(basename+"_tilde", metadata, Hk);
    }
  } else if (format=="HAPS") {
    if(include_genotypes) {
      plink::write_haps(basename+"_original", metadata, H);
    } else {
      plink::write_haps(basename, metadata, Hk);
    }
  } else {
    cerr << "Error: unrecognized format " << format << endl;
  }
}

void Knockoffs::writeAncestries(const string _out_file_name) const {
  string out_file_name = _out_file_name;
  out_file_name.append("_anc.txt");
  ofstream outfile(out_file_name.c_str());
  if (!outfile.is_open()) {
    cout << "Problem creating the output file: " << out_file_name;
    cout <<"Either the directory does not exist or you do not have write permissions." << endl;
  }
  const int num_windows = metadata.windows.num_windows;
  outfile << "# " << num_haps << endl;
  outfile << "# " << num_windows << endl;
  for(int i=0; i<num_haps; i++) {
    for(int w=0; w<num_windows; w++) {
      outfile << "# ID " << metadata.sample_filter.ID[i/2] << " " << (i % 2);
      outfile << " window " << w << " ";
      outfile << metadata.windows.start[w] << " " << metadata.windows.end[w] << endl;
      for(int k=0; k<K; k++) {
        outfile << alpha_global[i][k];
        if(k+1<K) outfile <<" ";
      }
      outfile << endl;
    }
  }
  outfile.close();
  cout << "Reconstructed ancestries written to:" << endl << "  " << out_file_name << endl;
}

// void Knockoffs::writeZ(const string _out_file_name) const {
//   string out_file_name = _out_file_name;
//   out_file_name.append("_Z.txt");
//   ofstream outfile(out_file_name.c_str());
//   if (!outfile.is_open()){
//     cout << "Problem creating the output file: " << out_file_name;
//     cout <<"Either the directory does not exist or you do not have write permissions." << endl;
//   }
//   for(int i=0; i<std::min(20000,num_haps); i++) {
//     for(int j=0; j<num_snps; j++) {
//       outfile << Z[i][j];
//       if(j+1<num_snps) outfile <<" ";
//     }
//     outfile << endl;
//   }
//   outfile.close();
//   cout << "Latent Markov chains written to:" << endl << "  " << out_file_name << endl;
// }

// void Knockoffs::writeZk(const string _out_file_name) const {
//   string out_file_name = _out_file_name;
//   out_file_name.append("_Zk.txt");
//   ofstream outfile(out_file_name.c_str());
//   if (!outfile.is_open()){
//     cout << "Problem creating the output file: " << out_file_name;
//     cout <<"Either the directory does not exist or you do not have write permissions." << endl;
//   }
//   for(int i=0; i<std::min(20000,num_haps); i++) {
//     for(int j=0; j<num_snps; j++) {
//       outfile << Zk[i][j];
//       if(j+1<num_snps) outfile <<" ";
//     }
//     outfile << endl;
//   }
//   outfile.close();
//   cout << "Latent Markov chains written to:" << endl << "  " << out_file_name << endl;
// }

void Knockoffs::writeGroups(const string _out_file_name) const {
  string out_file_name = _out_file_name;
  out_file_name.append("_grp.txt");
  ofstream outfile(out_file_name.c_str());
  if (!outfile.is_open()){
    cout << "Problem creating the output file: " << out_file_name;
    cout <<"Either the directory does not exist or you do not have write permissions." << endl;
  }
  outfile << "SNP Group" << endl;
  for(int j=0; j<num_snps; j++) {
    outfile << metadata.legend_filter.ID[j] << " " << groups[j] << endl;
  }
  outfile.close();
  cout << "SNP groups written to:" << endl << "  " << out_file_name << endl;
}

void Knockoffs::writeWindows(const string _out_file_name) const {
  string out_file_name = _out_file_name;
  out_file_name.append("_windows.txt");
  ofstream outfile(out_file_name.c_str());
  if (!outfile.is_open()){
    cout << "Problem creating the output file: " << out_file_name;
    cout <<"Either the directory does not exist or you do not have write permissions." << endl;
  }
  outfile << "SNP Window" << endl;
  for(int j=0; j<num_snps; j++) {
    outfile << metadata.legend_filter.ID[j] << " " << metadata.windows.window[j] << endl;
  }
  outfile.close();
  cout << "SNP windows written to:" << endl << "  " << out_file_name << endl;
}

int Knockoffs::num_partitions() const {
  return partitions.size();
}

void Knockoffs::set_Hk(const int i, vector<int> & hk) {
  Hk[i] = chaplotype(hk);
}

int Knockoffs::num_groups() const {
  return(groups.back()+1);
}

//////////////////////////////////
// Master methods for knockoffs //
//////////////////////////////////

void Knockoffs::generate(const int num_threads) {
  // Reset knockoff generator
  reset_knockoffs();

  // Process related individuals
  if(metadata.related.size()>0) {
    //ProfilerStart("/home/users/msesia/Workspace/popstruct/analysis_related/knockoffs/profile.log");
    generate_related(num_threads);
    //ProfilerStop();
  }

  // Process unrelated individuals
  if(metadata.unrelated.size()>0) {
    generate_unrelated(num_threads);
  }

}

void Knockoffs::generate_related(const int num_threads) {
  // List non-trivial IBD-sharing haplotype clusters
  set< vector<int> > related_families;
  for(int f=0; f<metadata.related_families.size(); f++) {
    // Skip trivial clusters
    if(metadata.related_families[f].size() <= 1) continue;
    // List of haplotypes
    vector<int> cluster_tmp;
    for(int i=0; i<metadata.related_families[f].size(); i++) {
      cluster_tmp.push_back(metadata.related_families[f][i]);
    }
    related_families.insert(cluster_tmp);
  }
  vector< vector<int> > related_families_vec(related_families.begin(), related_families.end());
  int num_clusters = related_families_vec.size();

  // Print header
  cout << endl << "Generating related knockoffs for chromosome " << metadata.get_chr_id();
  cout <<" (" << metadata.related.size() << " haplotypes in "<< num_clusters << " families, ";
  cout << num_snps << " variants in " << num_groups() << " groups)" << ":" << endl;

  if(DEBUG) {
    cout << "[DEBUG] IBD-sharing haplotype clusters:" << endl;
    int clust_id = 0;
    for(auto cluster : related_families) {
      cout << "Cluster " << clust_id << ": ";
      for(int i=0; i<cluster.size(); i++) {
        cout << cluster[i] << " ";
      }
      cout << endl;
      clust_id++;
    }
  }

  // Divide the clusters among the workers evenly by size
  std::sort(related_families_vec.begin(), related_families_vec.end(),
            [](const vector<int> & a, const vector<int> & b){ return a.size() > b.size(); });
  vector< vector< vector<int> > > related_cluster_assignments(num_threads);
  int w = 0;
  for(int i=0; i<related_families.size(); i++) {
    if(w>=num_threads) w=0;
    related_cluster_assignments[w].push_back(related_families_vec[i]);
    w++;
  }

  // // Debug: print assignments
  // cout << "[DEBUG] Worker assignments:" << endl;
  // for(int w=0; w<num_threads; w++ ) {
  //   cout << "- Worker " << w << ":" << endl;
  //   for(auto cluster : related_cluster_assignments[w]) {
  //     cout << "Cluster: ";
  //     for(int i=0; i<cluster.size(); i++) {
  //       cout << cluster[i] << " ";
  //     }
  //     cout << endl;
  //   }
  // }

  // Initialize progress bar
  vector<int> progress(num_threads,0);
  int n_steps = std::max(std::min(100, num_clusters), num_threads);
  int n_steps_each = std::max(1, n_steps / num_threads);
  cout << "|" << flush;
  for(int s=0; s<n_steps; s++) cout << ".";
  cout << "|" << endl;
  cout << "|" << flush;

  if(num_threads>1) {

    // Initialize workers
    vector<boost::thread> workers;
    // Random number generators for each worker
    vector<boost::random::taus88> rng(num_threads);

    // Generate knockoffs for blocks of individuals in parallel
    for(int w=0; w<num_threads; w++ ) {
      // Seed the random number generator
      rng[w].seed(seed + w);

      // Create worker
      workers.push_back(boost::thread(&Knockoffs::generate_related_worker, this,
                                      boost::ref(related_cluster_assignments[w]), boost::ref(rng[w]),
                                      w, boost::ref(progress), n_steps_each));
    }
    // Launch workers
    for(int w=0; w<num_threads; w++ ) {
      workers[w].join();
    }
  } else {
    // Initialize random number generator
    boost::random::taus88 rng;
    rng.seed(seed);
    // Generate related knockoffs, cluster by cluster
    generate_related_worker(related_cluster_assignments[0], rng, 0, boost::ref(progress), n_steps_each);
  }

  // Finalize progress bar
  for(int s=std::accumulate(progress.begin(),progress.end(),0); s<n_steps; s++) {
    cout << "=";
  }
  cout << "|" << endl;

}

void Knockoffs::generate_related_worker(const vector< vector<int> > & related_families,
                                        boost::random::taus88 & rng,
                                        const int wid, vector<int> &progress, const int n_steps_each) {

  // Initialize progress bar
  int progress_period = std::max(1, (int)std::ceil((double)(related_families.size()) / (double)(n_steps_each)));

  for(int i=0; i<related_families.size(); i++) {
    // Debug
    print_debug_log(wid, related_families[i], false);

    generate_related_cluster(related_families[i], rng);

    // Update progress bar
    if((i+1) % progress_period == 0) {
      if(progress[wid]<n_steps_each) {
        cout << "=" << flush;
        progress[wid]++;
      }
    }
  }
  print_debug_log(wid, related_families.back(), true);

  // Finalize progress bar
  for(int s=progress[wid]; s<n_steps_each; s++) {
    cout << "=" << flush;
    progress[wid]++;
  }

  //cout << "[DEBUG] worker " << wid << " finished." << endl;
}

void Knockoffs::generate_related_cluster(const vector<int> & related_cluster, boost::random::taus88 & rng) {
  if(DEBUG) {
    cout << endl << "[DEBUG] in generate_related_cluster( ";
    for(int i : related_cluster) cout << i <<" " << metadata.sample_filter.ID[i/2] << ";  ";
    cout << ")" << endl;
    //metadata.print_segments();
  }

  // Initialize IBD cluster
  set<IbdSeg> ibd_segments;
  for(auto i : related_cluster) {
    for(auto ibd : metadata.ibd_segments[i]) {
      ibd_segments.insert(ibd);
    }
  }
  IbdCluster ibd_cluster(ibd_segments);
  
  // Sample from posterior using BP
  FamilyBP familyBP(related_cluster, ibd_cluster, H, references_global, alpha_global, mut_rates, b, 
                    metadata, rng);
  vector< vector<int> > _Z;
  familyBP.run(_Z);

  // Generate related MC knockoffs
  FamilyKnockoffs familyKnockoffs(related_cluster, ibd_cluster, _Z, alpha_global, b, rng);
  vector< vector<int> > _Zk;
  familyKnockoffs.run(elements, groups, _Zk);

  // // DEBUG: trivial knockoffs
  // vector< vector<int> > _Zk = _Z;

  // // DEBUG: store Z
  // for(int i=0; i<related_cluster.size(); i++) {
  //   int i_abs = related_cluster[i];
  //   for(int j=0; j<num_snps; j++) {
  //     Z[i_abs][j] = _Z[i][j];
  //     Zk[i_abs][j] = _Zk[i][j];
  //   }
  // }
    
  // Generate emissions
  vector< vector<int> > _Hk(related_cluster.size(), vector<int>(num_snps,-1));
  for(int i=0; i<related_cluster.size(); i++) {
    const vector<int> ref = references_global[related_cluster[i]];
    sample_emission(_Zk[i], ref, _Hk[i], rng);
  }

  // // DEBUG: Make sure Hk in consistent on IBD segments

  // const vector<int> indices = related_cluster;
  // vector3i neighbors; // (id1,j,id2)
  // // Initialize list of neighbors
  // for(int i=0; i<num_haps; i++) {
  //   vector< vector<int> > tmp_v;
  //   vector<int> tmp_empty;
  //   for(int j=0; j<num_snps; j++) {
  //     tmp_v.emplace_back(tmp_empty);
  //   }
  //   neighbors.push_back(tmp_v);
  // }

  // // Fill list of neighbors
  // for(int s=0; s<ibd_cluster.size(); s++) {
  //   IbdSeg segment;
  //   ibd_cluster.get(s,segment);
  //   vector<int> segment_indices;
  //   for(int id : segment.indices) {
  //     // Find id in list of ids for this family
  //     auto it = std::find(indices.begin(), indices.end(), id);
  //     int i = std::distance(indices.begin(), it);
  //     segment_indices.push_back(i);
  //   }
  //   for(int i1 : segment_indices) {
  //     for(int j=segment.j_min; j<=segment.j_max; j++) {
  //       for(int i2 : segment_indices) {
  //         if(i1 != i2) {
  //           neighbors[i1][j].push_back(i2);
  //         }
  //       }
  //     }
  //   }
  // }

  // for(int i=0; i<related_cluster.size(); i++) {
  //   double agreement_num = 0;
  //   double agreement_den = 0;
  //   for(int j=0; j<num_snps; j++) {
  //     for(int i2 : neighbors[i][j]) {
  //       int i2_abs = indices[i2];
  //       if(_Hk[i][j]==_Hk[i2][j]) {
  //         agreement_num++;
  //       }
  //       agreement_den++;
  //     }
  //   }
  //   if(agreement_den==0) agreement_den = -1;
  //   int i_abs = indices[i];
  //   const string i_name = metadata.sample_filter.ID[i_abs/2];
  //   cout << "[DEBUG] IBD agreement (Hk) for haplotype " << i << "(" << i_name << "): ";
  //   cout  << agreement_num/agreement_den << endl;
  // }

  // Store results
  for(int i=0; i<related_cluster.size(); i++) {
    set_Hk(related_cluster[i], _Hk[i]);
  }
}

void Knockoffs::generate_unrelated(const int num_threads) {
  // Print header
  cout << endl << "Generating unrelated knockoffs for chromosome " << metadata.get_chr_id();
  cout <<" (" << metadata.unrelated.size() << " haplotypes; " << num_snps << " variants in ";
  cout << num_groups() << " groups)" << ":" << endl;

  // Assign individuals to workers
  vector<int> unrelated;
  std::copy(metadata.unrelated.begin(), metadata.unrelated.end(), std::back_inserter(unrelated));
  int num_unrelated = unrelated.size();
  vector< vector<int> > i_list;
  for(int w=0; w<num_threads; w++ ) {
    int i_start = w*num_unrelated/num_threads;
    int i_end = (w+1)*num_unrelated/num_threads;
    if(w==(num_threads-1)) i_end = num_unrelated;
    vector<int> i_list_tmp;
    for(int i=i_start; i<i_end; i++) i_list_tmp.push_back(unrelated[i]);
    i_list.push_back(i_list_tmp);
  }

  // Initialize workers
  vector<boost::thread> workers;
  // Random number generators for each worker
  vector<boost::random::taus88> rng(num_threads);

  // Initialize progress bar
  vector<int> progress(num_threads,0);
  int n_steps = std::min(100, num_unrelated);
  int n_steps_each = std::max(1, n_steps / num_threads);
  cout << "|" << flush;
  for(int s=0; s<n_steps; s++) cout << ".";
  cout << "|" << endl;
  cout << "|" << flush;

  // Progress monitor
  if(num_threads>1) {
    // Generate knockoffs for blocks of individuals in parallel
    for(int w=0; w<num_threads; w++ ) {
      // Seed the random number generator
      rng[w].seed(seed + w);

      // Create worker
      workers.push_back(boost::thread(&Knockoffs::generate_unrelated_worker, this, i_list[w],
                                      boost::ref(rng[w]), w, boost::ref(progress), n_steps_each));
    }
    // Launch workers
    for(int w=0; w<num_threads; w++ ) {
      workers[w].join();
    }
  } else {
    // Generate knockoffs on current thread
    generate_unrelated_worker(i_list[0], rng[0], 0, boost::ref(progress), n_steps_each);
  }

  // Finalize progress bar
  for(int s=std::accumulate(progress.begin(),progress.end(),0); s<n_steps; s++) {
    cout << "=";
  }
  cout << "|" << endl;

  cout << endl;
}

void Knockoffs::generate_unrelated_worker(const vector<int> & i_list, boost::random::taus88 & rng,
                                          const int wid, vector<int> &progress, const int n_steps_each) {

  // Initialize progress bar
  int progress_period = std::max(1, (int)std::ceil((double)(i_list.size()) / (double)(n_steps_each)));
  // cout << " steps_each = " << n_steps_each << endl;
  // cout << " i_list.size() = " << i_list.size() << endl;
  // cout << " progress_period = " << progress_period << endl;

  vector<int> z(num_snps);
  vector<int> zk(num_snps);
  vector<int> hk_tmp(num_snps);
  vector<int> hk(num_snps);

  vector< vector<double> > backward(num_snps, vector<double> (K,1.0/K));
  vector< vector<double> > forward(num_snps, vector<double> (K,1.0/K));
  vector< vector<double> > new_alpha(num_windows, vector<double>(K,0));

  for(int i : i_list) {
    for(int w=0; w<num_windows; w++) {     
      // Local references within this window
      const vector<int>& ref = references_local[i][w];
      const int j_start = metadata.windows.start[w];
      const int j_end = metadata.windows.end[w];

      // Sample Z|H
      sample_posterior_MC(i, z, backward, w, rng);

      // Sample Zk|Z
      sample_knockoff_MC(i, z, zk, w, rng);

      // // DEBUG: trivial knockoffs
      // zk = z;

      // // DEBUG: store Z
      // for(int j=j_start; j<j_end; j++) {
      //   Z[i][j] = z[j];
      //   Zk[i][j] = zk[j];
      // }

      // Sample Hk|Zk
      sample_emission_copula(i, z, zk, ref, hk_tmp, 0.5, j_start, j_end, rng);
      for(int j=j_start; j<j_end; j++) {
        hk[j] = hk_tmp[j];
      }
    }

    set_Hk(i, hk);

    // Update progress bar
    if((i+1) % progress_period == 0) {
      if(progress[wid]<n_steps_each) {
        cout << "=" << flush;
        progress[wid]++;
      }
    }
  }

  // Finalize progress bar
  for(int s=progress[wid]; s<n_steps_each; s++) {
    cout << "=" << flush;
    progress[wid]++;
  }

}

///////////////////
// Model fitting //
///////////////////

void Knockoffs::EM(const int num_threads) {
  cout << endl << "Estimating HMM parameters with EM for chromosome " << metadata.get_chr_id();
  cout <<" (" << metadata.sample_filter.size() << " haplotypes; " << num_snps << " variants)" << ":" << endl;

  cout << "--------------------------------------------------" << endl;
  cout << std::setw(10) << "Iteration" << "\t" << std::setw(10) << "Rho" << "\t";
  cout << std::setw(10) << "Delta" << endl;
  cout << "--------------------------------------------------" << endl;

  const double delta_min = 5e-2;
  const int it_max = 20;
  bool converged = false;
  int it = 0;
  for(it=0; it<it_max; it++) {
    const bool print_period = (it % 1 == 0);
    if(print_period) cout << std::setw(10) << std::fixed << it+1 << "\t" << flush;
    const double delta_parameters = EM_step(num_threads, false);
    if(print_period) {
      cout << std::setw(10) << setprecision(3) << std::scientific << rho << "\t";
      cout << std::setw(10) << setprecision(3) << delta_parameters << endl;
    }
    if(delta_parameters < delta_min) {
      converged = true;
      if(!print_period) {
        cout << std::setw(10) << setprecision(0) << std::fixed << it+1 << "\t" << flush;
        cout << std::setw(10) << setprecision(3) << std::scientific << rho << "\t";
        cout << std::setw(10) << setprecision(3) << delta_parameters << endl;
      }
      break;
    }
  }

  cout << "--------------------------------------------------" << endl;

  if(converged) {
    cout << "EM converged after " << it+1 << " iterations." << endl;
  } else {
    cout << "EM did not converge after " << it << " iterations." << endl;
  }
  cout << endl;
}

void Knockoffs::EM_worker(const vector<int> & i_list, boost::random::taus88 & rng,
                          const int wid, vector<int> &progress, const int n_steps_each,
                          const bool show_progress, vector<double>& gamma, vector<double>& Xi) const {

  // Initialize progress bar
  int progress_period = std::max(1, (int)std::ceil((double)(i_list.size()) / (double)(n_steps_each)));

  vector< vector<double> > backward(num_snps, vector<double> (K,1.0/K));
  vector< vector<double> > forward(num_snps, vector<double> (K,1.0/K));
  vector< vector<double> > marginal(num_snps, vector<double> (K,1.0/K));

  gamma.resize(num_snps, 0);
  Xi.resize(num_snps, 0);

  for(int i : i_list) {

    for(int w=0; w<num_windows; w++) {
      // Local references within this window
      const vector<int>& ref = references_local[i][w];
      const int j_start = metadata.windows.start[w];
      const int j_end = metadata.windows.end[w];

      // Window padding
      int j_start_pad = j_start;
      int j_end_pad = j_end;
      if(w>0) j_start_pad = metadata.windows.start[w-1];
      if(w<(num_windows-1)) j_end_pad = metadata.windows.end[w+1];

      // Compute forward and backward weights
      compute_hmm_backward(i, rng, ref, backward, w, j_start_pad, j_end_pad);
      compute_hmm_forward(i, rng, ref, forward, w, j_start_pad, j_end_pad);

      // Compute posterior marginal
      for(int j=j_start; j<j_end; j++) {
        double marginal_sum = 0;
        for(int k=0; k<K; k++) {
          marginal[j][k] = backward[j][k] * forward[j][k];
          marginal_sum += marginal[j][k];
        }
        for(int k=0; k<K; k++) {
          marginal[j][k] /= marginal_sum;
        }
      }

      // Compute gamma (expected mutations under this model)
      for(int j=j_start; j<j_end; j++) {
        const int h_real = H[i][j];
        for(int k=0; k<K; k++) {
          if(h_real != H[ref[k]][j]) {
            gamma[j] += marginal[j][k];
          }
        }
      }

      // Compute the normalization constant for xi
      vector<double> xi_norm(num_snps, 0);
      for(int j=j_start; j<j_end; j++) {
        const double bj = b[j];
        const double aj = (1.0-bj) / (double)(K);
        const int h_real = H[i][j];
        const double mut_j = mut_rates[j];
        for(int k=0; k<K; k++) {
          const int h_ref = H[ref[k]][j];
          const double fjk = (h_real==h_ref) ? (1.0-mut_j) : mut_j;
          xi_norm[j] += aj * fjk * backward[j][k];
          if(j>0) xi_norm[j] += bj * forward[j-1][k] * fjk * backward[j][k];
        }
      }

      // Compute the unnormalized diagonal of Xi (expected transitions into same state)
      for(int j=j_start; j<j_end; j++) {
        const double bj = b[j];
        const double aj = (1.0-bj) / (double)(K);
        const int h_real = H[i][j];
        const double mut_j = mut_rates[j];
        for(int k=0; k<K; k++) {
          const int h_ref = H[ref[k]][j];
          const double fjk = (h_real==h_ref) ? (1.0-mut_j) : mut_j;
          if(j>0) Xi[j] += (aj+bj) * forward[j-1][k] * backward[j][k] * fjk / xi_norm[j];
          else Xi[j] += aj * backward[j][k] * fjk / xi_norm[j];
        }
      }

    }

    // Update progress bar
    if(show_progress) {
      if((i+1) % progress_period == 0) {
        if(progress[wid]<n_steps_each) {
          cout << "=" << flush;
          progress[wid]++;
        }
      }
    }

    }

  // Finalize progress bar
  if(show_progress) {
    for(int s=progress[wid]; s<n_steps_each; s++) {
      cout << "=" << flush;
      progress[wid]++;
    }
  }

}

double Knockoffs::EM_step(const int num_threads, const bool show_progress) {
  //cout << "[DEBUG] in Knockoffs::EM()" << endl;

  vector< vector<double> > gamma_list(num_threads, vector<double> (num_snps, 0));
  vector< vector<double> > Xi_list(num_threads, vector<double> (num_snps, 0));

  // Determine which samples to use for EM
  const int n = num_haps;
  vector<double> i_list_full(num_haps);
  std::iota(i_list_full.begin(), i_list_full.end(), 0);
  const double n_inv = 1.0 / (double)(n);

  // Assign samples to workers
  vector< vector<int> > i_list;
  for(int w=0; w<num_threads; w++ ) {
    int i_start = w*n/num_threads;
    int i_end = (w+1)*n/num_threads;
    if(w==(num_threads-1)) i_end = n;
    vector<int> i_list_tmp;
    for(int i=i_start; i<i_end; i++) i_list_tmp.push_back(i_list_full[i]);
    i_list.push_back(i_list_tmp);
  }

  // Initialize workers
  vector<boost::thread> workers;
  // Random number generators for each worker
  vector<boost::random::taus88> rng(num_threads);

  // Initialize progress bar
  vector<int> progress(num_threads,0);
  int n_steps = std::min(100, n);
  int n_steps_each = std::max(1, n_steps / num_threads);
  if(show_progress) {
    cout << "|" << flush;
    for(int s=0; s<n_steps; s++) cout << ".";
    cout << "|" << endl;
    cout << "|" << flush;
  }

  // Progress monitor
  if(num_threads>1) {
    // Generate knockoffs for blocks of individuals in parallel
    for(int w=0; w<num_threads; w++ ) {
      // Seed the random number generator
      rng[w].seed(seed + w);

      // Create worker
      workers.push_back(boost::thread(&Knockoffs::EM_worker, this, i_list[w], boost::ref(rng[w]),
                                      w, boost::ref(progress), n_steps_each, show_progress,
                                      boost::ref(gamma_list[w]), boost::ref(Xi_list[w])));
    }
    // Launch workers
    for(int w=0; w<num_threads; w++ ) {
      workers[w].join();
    }
  } else {
    // Generate knockoffs on current thread
    EM_worker(i_list[0], rng[0], 0, boost::ref(progress), n_steps_each, show_progress,
              boost::ref(gamma_list[0]), boost::ref(Xi_list[0]));
  }

  // Finalize progress bar
  if(show_progress) {
    for(int s=std::accumulate(progress.begin(),progress.end(),0); s<n_steps; s++) {
      cout << "=";
    }
    cout << "|" << endl;
    cout << endl;
  }

  // Collect results
  vector<double> gamma(num_snps, 0);
  vector<double> Xi(num_snps, 0);
  for(int j=0; j<num_snps; j++) {
    for(int w=0; w<gamma_list.size(); w++) {
      gamma[j] += gamma_list[w][j] * n_inv;
      Xi[j] += Xi_list[w][j] * n_inv;
    }
  }

  // Make sure the mutation rates are not too small or too large
  for(int j=0; j<num_snps; j++) {
    gamma[j] = std::max(1e-6, gamma[j]);
    gamma[j] = std::min(1e-3, gamma[j]);
  }

  // Solve for the recombination rate multiplier
  //const double new_rho = EM_solve_rho(rho, Xi);
  const double new_rho = rho;

  // Compare new parameters with old parameters
  const double delta_recomb = std::abs(new_rho - rho);
  double delta_mutation = 0;
  for(int j=0; j<num_snps; j++) {
    const double delta_j = std::abs(gamma[j]-mut_rates[j]) / std::max(1e-5, mut_rates[j]);
    delta_mutation += delta_j*delta_j;
  }
  delta_mutation = std::sqrt(delta_mutation / (double)(num_snps));
  const double delta_parameters = std::max(delta_recomb, delta_mutation);

  // Update HMM parameters
  update_hmm(new_rho, gamma);

  // Return
  return(delta_parameters);
}

double Knockoffs::EM_Psi_rho(const double _rho, const vector<double>& Xi) const {
  //cout << endl << "[DEBUG]: in Knockoffs::EM_Psi_rho(" << _rho << ")" << endl;
  // TODO: make this robust to divisions by zero
  double psi = 0;
  for(int j=1; j<num_snps; j++) {
    const double dj = genetic_dist[j];
    const double bj = std::exp(-_rho * dj);
    //cout << "j = " << j << ", dj = " << dj << ", bj = " << bj << endl;
    const double new_psi = dj * bj * ( (K-1.0)/(1.0+(K-1.0)*bj) + 1.0/(1.0-bj) ) * Xi[j];
    //const double new_psi = dj * bj * ( (K-1.0)/(1.0+(K-1.0)*bj) + 1.0/(1.0-bj) );
    //cout << "new_psi = " << new_psi << endl;
    psi += new_psi;
  }
  return(psi);
}

double Knockoffs::EM_Phi_rho(const double _rho) const {
  double d_bar = 0;
  for(int j=1; j<num_snps; j++) {
    d_bar += genetic_dist[j];
  }
  d_bar /= (double)(num_snps);

  // TODO: make this robust to divisions by zero
  double phi = 0;
  for(int j=1; j<num_snps; j++) {
    const double dj = genetic_dist[j];
    const double bj_bar = std::exp(-_rho * (dj-d_bar));
    const double bj = std::exp(-_rho * dj);
    phi += dj * bj_bar / (1.0-bj);
  }
  return(phi);
}

double Knockoffs::EM_solve_rho(const double _rho, const vector<double>& Xi) const {
  //cout << "[DEBUG] in EM_solve_rho(" << _rho << ")" << endl;
  double d_bar = 0;
  for(int j=1; j<num_snps; j++) {
    d_bar += genetic_dist[j];
  }
  d_bar /= (double)(num_snps);
  double new_rho = _rho;
  for(int rep=0; rep<50; rep++) {
    const double psi = EM_Psi_rho(new_rho, Xi);
    const double phi = EM_Phi_rho(new_rho);
    new_rho = - std::log(psi/phi) / d_bar;
    //cout << "  new_rho " << new_rho << endl;
  }
  return(new_rho);
}

void Knockoffs::EM_solve_alpha(const vector<double>& eta, const vector<double>& xi,
                               vector<double>& new_alpha, const int j_start, const int j_end) const {

  // cout << "[DEBUG] in EM_solve_new_alpha()" << endl;
  const double eta_bar = std::accumulate(xi.begin()+j_start, xi.begin()+j_end, 0.0);
  const double K_inv = 1.0/(double)(K);

  for(int rep=0; rep<10; rep++) {
    double new_alpha_sum = 0;
    for(int k=0; k<K; k++) {
      double sum_tmp = 0;
      for(int j=j_start; j<j_end; j++) {
        const double bj = b[j];
        sum_tmp += (1.0-bj) / ((1.0-bj)*K_inv + bj) * xi[k];
      }
      new_alpha[k] = eta[k] - eta_bar + K_inv * sum_tmp;
      new_alpha_sum += new_alpha[k];
    }
    for(int k=0; k<K; k++) {
      new_alpha[k] /= new_alpha_sum;
    }
  }
}

////////////////
// Algorithms //
////////////////

void Knockoffs::compute_hmm_backward(const int i, boost::random::taus88 & rng, const vector<int>& ref,
                                     vector< vector<double> >& backward, const int w,
                                     const int j_start, const int j_end) const {

  // Compute backward weights
  std::fill(backward[j_end-1].begin(), backward[j_end-1].end(), 1.0/K);
  for(int j=j_end-2; j>=j_start; j--) {
    double backward_const = 0.0;
    const int h = H[i][j+1];
    const double bj = b[j+1];
    const double mut_j = mut_rates[j+1];
    for(int k=0; k<K; k++) {
      const double a_j_k = (1.0-bj) * alpha[i][w][k];                  // a[j+1][k]
      const int h_ref = H[ref[k]][j+1];
      if(h==h_ref) {
        backward_const += backward[j+1][k] * a_j_k * (1.0-mut_j);
      } else {
        backward_const += backward[j+1][k] * a_j_k * mut_j;
      }
    }
    double backwardSum = 0.0;
    for(int k=0; k<K; k++) {
      const int h_ref = H[ref[k]][j+1];
      if(h==h_ref) {
        backward[j][k] = backward_const + bj * backward[j+1][k] * (1.0-mut_j);
      } else {
        backward[j][k] = backward_const + bj * backward[j+1][k] * mut_j;
      }
      backwardSum += backward[j][k];
    }
    for(int k=0; k<K; k++) {
      backward[j][k] /= backwardSum;
    }
  }
}

void Knockoffs::compute_hmm_backward(const int i, boost::random::taus88 & rng, const vector<int>& ref,
                                     vector< vector<double> >& backward) const {
  compute_hmm_backward(i, rng, ref, backward, 0, 0, num_snps);
}

void Knockoffs::compute_hmm_forward(const int i, boost::random::taus88 & rng, const vector<int>& ref,
                                    vector< vector<double> >& forward) const {
  compute_hmm_forward(i, rng, ref, forward, 0, 0, num_snps);
}

void Knockoffs::compute_hmm_forward(const int i, boost::random::taus88 & rng, const vector<int>& ref,
                                    vector< vector<double> >& forward,
                                    const int w, const int j_start, const int j_end) const {


  // Compute forward weights
  double forwardSum_0 = 0.0;
  const int h0 = H[i][j_start];
  for(int k=0; k<K; k++) {
    const int h0_ref = H[ref[k]][j_start];
    if(h0==h0_ref) forward[j_start][k] = alpha[i][w][k] * (1.0-mut_rates[j_start]);
    else forward[j_start][k] = alpha[i][w][k] * mut_rates[j_start];
    forwardSum_0 += forward[j_start][k];
  }
  for(int k=0; k<K; k++) {
    forward[j_start][k] /= forwardSum_0;
  }

  for(int j=j_start+1; j<j_end; j++) {
    const int h = H[i][j];
    const double bj = b[j];
    const double mut_j = mut_rates[j];
    double forwardSum = 0.0;
    for(int k=0; k<K; k++) {
      const double a_j_k = (1.0-bj) * alpha[i][w][k];                  // a[j][k]
      forward[j][k] = a_j_k + bj * forward[j-1][k];
      const int h_ref = H[ref[k]][j];
      if(h==h_ref) {
        forward[j][k] *= (1.0-mut_j);
      } else {
        forward[j][k] *= mut_j;
      }
      forwardSum += forward[j][k];
    }
    for(int k=0; k<K; k++) {
      forward[j][k] /= forwardSum;
    }
  }
}

void Knockoffs::estimate_motif_prevalence(const int i, vector<int> & z,
                                          const vector< vector<double> >& backward,
                                          vector< vector<double> >& forward,
                                          vector< vector<double> >& new_alpha,
                                          boost::random::taus88 & rng) const {

  //cout << "[DEBUG] in estimate_motif_prevalence()" << endl;

  vector<double> xi_norm(num_snps, 0);
  vector<double> xi(num_snps,0);
  vector<double> eta(K,0);
  vector< vector<double> > marginal(num_snps, vector<double> (K,1.0/K));

  for(int w=0; w<num_windows; w++) {
    // Local references within this window
    const vector<int>& ref = references_local[i][w];
    const int j_start = metadata.windows.start[w];
    const int j_end = metadata.windows.end[w];

    // Window padding
    int j_start_pad = j_start;
    int j_end_pad = j_end;
    if(w>0) j_start_pad = metadata.windows.start[w-1];
    if(w<(num_windows-1)) j_end_pad = metadata.windows.end[w+1];

    // Compute forward weights
    compute_hmm_forward(i, rng, ref, forward, w, j_start_pad, j_end_pad);

    // Compute posterior marginal
    for(int j=j_start; j<j_end; j++) {
      double marginal_sum = 0;
      for(int k=0; k<K; k++) {
        marginal[j][k] = backward[j][k] * forward[j][k];
        marginal_sum += marginal[j][k];
      }
      for(int k=0; k<K; k++) {
        marginal[j][k] /= marginal_sum;
      }
    }

    // Compute the normalization constant for xi
    for(int j=j_start; j<j_end; j++) {
      const double bj = b[j];
      const double aj = (1.0-bj) / (double)(K);
      const int h_real = H[i][j];
      const double mut_j = mut_rates[j];
      for(int k=0; k<K; k++) {
        const int h_ref = H[ref[k]][j];
        const double fjk = (h_real==h_ref) ? (1.0-mut_j) : mut_j;
        xi_norm[j] += aj * fjk * backward[j][k];
        if(j>0) xi_norm[j] += bj * forward[j-1][k] * fjk * backward[j][k];
      }
    }

    // Compute the unnormalized diagonal of Xi (expected transitions into same state)
    for(int j=j_start; j<j_end; j++) {
      const double bj = b[j];
      const double aj = (1.0-bj) / (double)(K);
      const int h_real = H[i][j];
      const double mut_j = mut_rates[j];
      for(int k=0; k<K; k++) {
        const int h_ref = H[ref[k]][j];
        const double fjk = (h_real==h_ref) ? (1.0-mut_j) : mut_j;
        if(j>0) xi[j] += (aj+bj) * forward[j-1][k] * backward[j][k] * fjk / xi_norm[j];
        else xi[j] += aj * backward[j][k] * fjk / xi_norm[j];
      }
    }

    // Compute eta (expected frequency of transitions into state k)
    for(int j=j_start; j<j_end; j++) {
      const double bj = b[j];
      const double aj = (1.0-bj) / (double)(K);
      const int h_real = H[i][j];
      const double mut_j = mut_rates[j];
      for(int k=0; k<K; k++) {
        const int h_ref = H[ref[k]][j];
        const double fjk = (h_real==h_ref) ? (1.0-mut_j) : mut_j;
        if(j>0) eta[k] += (aj+bj * forward[j-1][k]) * backward[j][k] * fjk / xi_norm[j];
        else eta[k] += aj * backward[j][k] * fjk / xi_norm[j];
      }
    }

    // Estimate ancestry
    EM_solve_alpha(eta, xi, new_alpha[w], j_start, j_end);
  }

}

void Knockoffs::sample_posterior_MC(const int i, vector<int> & z, vector< vector<double> >& backward,
                                    const int w, boost::random::taus88 & rng) const {

  // Local references within this window
  const vector<int>& ref = references_local[i][w];
  const int j_start = metadata.windows.start[w];
  const int j_end = metadata.windows.end[w];
  // Window padding
  int j_start_pad = j_start;
  int j_end_pad = j_end;
  if(w>0) j_start_pad = metadata.windows.start[w-1];
  if(w<(num_windows-1)) j_end_pad = metadata.windows.end[w+1];

  // Compute backward weights
  compute_hmm_backward(i, rng, ref, backward, w, j_start_pad, j_end_pad);

  // Forward sampling within this window
  vector<double> weights_mc(K);
  for(int k=0; k<K; k++) {
    const double a_0_k = (1.0-b[j_start_pad]) * alpha[i][w][k];
    const double theta_0_k = mutate(H[ref[k]][j_start_pad], mut_rates[j_start_pad]);
    if ( H[i][j_start_pad] == 1) {
      weights_mc[k] = a_0_k * theta_0_k * backward[j_start_pad][k];
    }
    else {
      weights_mc[k] = a_0_k * (1.0-theta_0_k) * backward[j_start_pad][k];
    }
  }
  z[j_start_pad] = weighted_choice(weights_mc, rng);

  for(int j=j_start_pad+1; j<j_end_pad; j++) {
    const double bj = b[j];
    const int h = H[i][j];
    const double mut_j = mut_rates[j];
    for(int k=0; k<K; k++) {
      const double a_j_k = (1.0-bj) * alpha[i][w][k];                    // a[j][k]
      const double theta_j_k = mutate(H[ref[k]][j], mut_j);           // theta[j][k]
      if (h == 1) {
        weights_mc[k] = (a_j_k + bj*(double)(k==z[j-1])) * theta_j_k * backward[j][k];
      }
      else {
        weights_mc[k] = (a_j_k + bj*(double)(k==z[j-1])) * (1.0-theta_j_k) * backward[j][k];
      }
    }
    z[j] = weighted_choice(weights_mc, rng);
  }
}

void Knockoffs::sample_knockoff_MC(const int i, const vector<int> & z, vector<int> & zk,
                                   boost::random::taus88 & rng) const {
  for(int w=0; w<num_windows; w++) {
    sample_knockoff_MC(i, z, zk, w, rng);
  }
}

void Knockoffs::sample_knockoff_MC(const int i, const vector<int> & z, vector<int> & zk,
                                   const int w, boost::random::taus88 & rng) const {

  // Local window
  const int j_start = metadata.windows.start[w];
  const int j_end = metadata.windows.end[w];
  // Window padding
  int j_start_pad = j_start;
  int j_end_pad = j_end;
  if(w>0) j_start_pad = metadata.windows.start[w-1];
  if(w<(num_windows-1)) j_end_pad = metadata.windows.end[w+1];
  // Boundary groups
  int g_start = groups[j_start_pad];
  int g_end = groups[j_end_pad-1]+1;

  // Intermediate variables for knockoff generation
  vector<double> weights_mc(K);
  vector<double> N(K);
  vector<double> N_old(K);
  std::vector<double> vstar;
  vector< vector<double> > vbar;

  // Initalize normalization functions
  std::fill(N.begin(), N.end(), 1.0);
  std::fill(N_old.begin(), N_old.end(), 1.0);
  const double N_min = 1.0e-10;

  for(int g=g_start; g<g_end; g++) {
    // Compute vstar
    int groupSize = elements[g].size();
    vstar.resize(groupSize,0.0);
    for(int j=groupSize-1; j>=0; j--) {
      if(j < groupSize-1) {
        vstar[j] = vstar[j+1] * b[elements[g][j+1]];
      }
      else {
        if(g < g_end-1) {
          vstar[j] = b[elements[g][0]];
        }
        else {
          vstar[j] = 0.0;
        }
      }
    }

    // Compute vbar matrix
    vbar = vector< vector<double> >(groupSize, std::vector<double> (K,0));
    for(int j=groupSize-1; j>=0; j--) {
      const int w = metadata.windows.window[j];
      double sum_a = 0;
      if(j < groupSize-1) {
        for(int l=0; l<K; l++) {
          sum_a += (1.0-b[elements[g][j+1]]) * alpha[i][w][l];
        }
      }
      for(int z=0; z<K; z++) {
        if(j < groupSize-1) {
          vbar[j][z] = vbar[j+1][z] * (sum_a + b[elements[g][j+1]]);
          vbar[j][z] += vstar[j+1] * (1.0-b[elements[g][j+1]]) * alpha[i][w][z];
        }
        else {
          if(g < g_end-1) {
            vbar[j][z] = (1.0-b[elements[g+1][0]]) * alpha[i][w][z];
          }
          else {
            vbar[j][z] = 1.0;
          }
        }
      }
    }

    // First variant in the group
    const int jg0 = elements[g][0];
    const int w0 = metadata.windows.window[0];

    // Precompute sum for partition function
    double N_sum = 0;
    for(int k=0; k<K; k++) {
      if( g==g_start ) {
        N_sum += (1.0-b[jg0]) * alpha[i][w0][k];
      }
      else {
        const int z0 = z[elements[g-1].back()];
        const int z0k = zk[elements[g-1].back()];
        double tmp = ((1.0-b[jg0])* alpha[i][w0][k] + b[jg0]*(double)(k==z0));
        tmp  = tmp * ((1.0-b[jg0])* alpha[i][w0][k] + b[jg0]*(double)(k==z0k));
        N_sum += tmp / N_old[k];
      }
    }

    // Compute partition function
    for(int k=0; k<K; k++) {
      if(g < g_end-1) {
        N[k] += vbar[0][k] * N_sum;
        if(g==0) {
          N[k] += vstar[0] * (1.0-b[jg0]) * alpha[i][w0][k];
        }
        else {
          const int z0 = z[elements[g-1].back()];
          const int z0k = zk[elements[g-1].back()];
          double tmp = ((1.0-b[jg0]) * alpha[i][w0][k] + b[jg0]*(double)(k==z0));
          tmp  = tmp * ((1.0-b[jg0]) * alpha[i][w0][k] + b[jg0]*(double)(k==z0k));
          N[k] += vstar[0] * tmp / N_old[k];
        }
      }
    }

    // Normalize partition function and make sure we avoid division by zero
    double N_norm = 0;
    for(int k=0; k<K; k++) {
      if(N[k] < N_min) {
        N[k] = N_min;
      }
      N_norm += N[k];
    }
    for(int k=0; k<K; k++) {
      N[k] /= N_norm;
      if(N[k] < N_min) {
        N[k] = N_min;
      }
    }

    // Compute sampling weights
    for(int j=0; j<groupSize; j++) {
      const int w = metadata.windows.window[j];
      std::fill(weights_mc.begin(), weights_mc.end(), 1.0);
      for(int k=0; k<K; k++) {
        if(j>0) {
          const int z0k = zk[elements[g][j-1]];
          weights_mc[k]*=((1.0-b[elements[g][j]])*alpha[i][w][k]+b[elements[g][j]]*(double)(k==z0k));
        }
        else {
          if( g==0 ) {
            weights_mc[k] *= (1.0-b[jg0]) * alpha[i][w][k];
          }
          else {
            const int z0 = z[elements[g-1].back()];
            const int z0k = zk[elements[g-1].back()];
            weights_mc[k] *= ((1.0-b[jg0]) * alpha[i][w][k] + b[jg0]*(double)(k==z0));
            weights_mc[k] *= ((1.0-b[jg0]) * alpha[i][w][k] + b[jg0]*(double)(k==z0k));
          }
          weights_mc[k] /= N_old[k];
        }
        if(g == g_end-1) {
          weights_mc[k] *= (vbar[j][0]);
        }
        else {
          const int z1 = z[elements[g+1][0]];
          weights_mc[k] *= (vbar[j][z1] + vstar[j] * (double)(k==z1));
        }
      }
      zk[elements[g][j]] = weighted_choice(weights_mc, rng);
    }

    for(int k=0; k<K; k++) {
      N_old[k] = N[k];
    }
  }
}

void Knockoffs::sample_emission(const vector<int> & z, const vector<int>& ref, vector<int> & output,
                                const int j_start, const int j_end, boost::random::taus88 & rng) const {
  std::vector<double> weights_binary(2,1.0);
  // Sample all snps using the provided global reference
  for(int j=j_start; j<j_end; j++) {
    const int motif_id = ref[z[j]];
    const int h_predicted = H[motif_id][j];
    weights_binary[1] = mutate(h_predicted, mut_rates[j]);
    weights_binary[0] = 1.0 - weights_binary[1];
    output[j] = weighted_choice(weights_binary, rng);
  }
}

void Knockoffs::sample_emission(const vector<int> & z, const vector<int>& ref_global, vector<int> & output,
                                boost::random::taus88 & rng) const {
  sample_emission(z, ref_global, output, 0, num_snps, rng);
}

void Knockoffs::sample_emission(const vector<int> & z, const vector< vector<int> >& ref_local, 
                                vector<int> & output, boost::random::taus88 & rng) const {
  for(int w=0; w<num_windows; w++) {
    // Local references within this window
    const int j_start = metadata.windows.start[w];
    const int j_end = metadata.windows.end[w];
    sample_emission(z, ref_local[w], output, 0, num_snps, rng);
  }
}

void Knockoffs::sample_emission_copula(const int i, const vector<int> & z, const vector<int> & zk,
                                       vector<int> & output, double lambda,
                                       boost::random::taus88 & rng) const {
  for(int w=0; w<num_windows; w++) {
    const vector<int>& ref = references_local[i][w];
    const int j_start = metadata.windows.start[w];
    const int j_end = metadata.windows.end[w];
    sample_emission_copula(i, z, zk, ref, output, lambda, j_start, j_end, rng);
  }
}

void Knockoffs::sample_emission_copula(const int i, const vector<int> & z, const vector<int> & zk,
                                       const vector<int>& ref, vector<int> & output, double lambda,
                                       const int j_start, const int j_end, boost::random::taus88 & rng) const {

  std::vector<double> weights_binary(2,1.0);

  for(int j=j_start; j<j_end; j++) {
    const int motif_id = ref[z[j]];
    const int motif_id_ko = ref[zk[j]];
    const int h_predicted = H[motif_id][j];
    const int h_predicted_ko = H[motif_id_ko][j];
    const double p_orig1 = mutate(h_predicted,mut_rates[j]); //marginal probability of a 1 given latent state
    const double p_ko1 = mutate(h_predicted_ko,mut_rates[j]); //marginal probability of a 1 given latent state

    //generate copula sampling weights
    double temp = p_orig1 - (1.0 - p_ko1);
    double temp2 = 0.0;
    if(temp < 0) {
      temp2 = -temp;
      temp = 0.0; // max(0, p_orig1 - (1 - p_ko1))
    }
    if(H[i][j] == 0) {
      weights_binary[1] = 1.0 - temp2 / (1.0 - p_orig1);
      weights_binary[0] = temp2 / (1.0 - p_orig1);
    } else {
      weights_binary[1] = temp / (p_orig1);
      weights_binary[0] = 1.0 - temp / (p_orig1);
    }
    weights_binary[0] = lambda * weights_binary[0] + (1 - lambda) * (1 - p_ko1);
    weights_binary[1] = lambda * weights_binary[1] + (1 - lambda) * p_ko1;

    output[j] = weighted_choice(weights_binary, rng);
  }
}

/////////
// HMM //
/////////

void Knockoffs::init_hmm(const int K, double _rho, double _lambda) {
  // Initialize recombination rates
  rho = _rho;
  rec_rates.clear();
  rec_rates.resize(num_snps,1);
  b.resize(num_snps,0);
  for(unsigned int j=1; j<num_snps; j++) {
    rec_rates[j] = rho * genetic_dist[j];
    b[j] = std::exp(-rec_rates[j]);
  }

  // Initialize mutation rates
  mut_rates.clear();
  mut_rates.resize(num_snps);
  for(unsigned int j=0; j<num_snps; j++) {
    mut_rates[j] = _lambda;
  }
}

void Knockoffs::init_hmm(const int K, double _rho) {
  // Initialize recombination rates
  rho = _rho;
  rec_rates.clear();
  rec_rates.resize(num_snps,1);
  b.resize(num_snps,0);
  for(unsigned int j=1; j<num_snps; j++) {
    rec_rates[j] = rho * genetic_dist[j];
    b[j] = std::exp(-rec_rates[j]);
  }

  // Compute MAF
  vector<double> maf(num_snps, 0);
  compute_maf(maf);

  // Initialize mutation rates
  mut_rates.clear();
  mut_rates.resize(num_snps);
  for(unsigned int j=0; j<num_snps; j++) {
    mut_rates[j] = std::min(std::max(maf[j], 1e-5), 1e-3);
  }
}

void Knockoffs::update_hmm(const double _rho, const vector<double>& _mut_rates) {
  // Initialize recombination rates
  rho = _rho;
  for(unsigned int j=1; j<num_snps; j++) {
    rec_rates[j] = rho * genetic_dist[j];
    b[j] = std::exp(-rec_rates[j]);
  }
  // Initialize mutation rates
  for(unsigned int j=0; j<num_snps; j++) {
    mut_rates[j] = _mut_rates[j];
  }
}

void Knockoffs::print_hmm() const {
  std::cout<<"Printing parameters of HMM for "<<num_snps<<" variants"<<std::endl;
  cout << "SNP\tDist\tRecomb\tMut" << endl;
  for(unsigned int j=0; j<100; j++) {
    cout << j << "\t";
    if(genetic_dist.size()>0) {
      std::cout<<genetic_dist[j]<<"\t"<<rec_rates[j];
    }
    std::cout<<"\t"<<mut_rates[j]<<endl;
  }
  std::cout<<std::endl;
}

void Knockoffs::writeHMM(const string _out_file_name) const {
  // Write site-specific recombination and mutation rates
  string out_file_name = _out_file_name;
  out_file_name.append("_hmm.txt");
  ofstream outfile(out_file_name.c_str());
  if (!outfile.is_open()){
    cout << "Problem creating the output file: " << out_file_name;
    cout <<"Either the directory does not exist or you do not have write permissions." << endl;
  }
  outfile << "Recombination Mutation" << endl;
  for(int j=0; j<num_snps; j++) {
    outfile << rec_rates[j] << " " << mut_rates[j] << endl;
  }
  outfile.close();
  cout << "HMM parameters written to:" << endl << "  " << out_file_name << endl;
}

void Knockoffs::load_hmm(const string hmm_file) {
  int _num_snps = count_lines(hmm_file) - 1;
  assert(_num_snps==num_snps);
  string buffer;
  vector < string > tokens;
  int line = -1;
  ifile fd_hmm(hmm_file);
  while (getline(fd_hmm, buffer, '\n')) {
    if(line>=0) {
      sutils::tokenize(buffer, tokens);
      rec_rates[line] = std::stod(tokens[0]);
      b[line] = std::exp(-rec_rates[line]);
      mut_rates[line] = std::max(1e-6,std::min(std::stod(tokens[1]),1e-3));
    }
    line ++;
  }
  fd_hmm.close();
}

void Knockoffs::print_debug_log(const int wid, const vector<int> & related_cluster, const bool finished) const {
  // DEBUG: write Z to text file
  string filename = logfile + "_wid"+std::to_string(wid)+".txt";
  ofstream myfile;
  myfile.open(filename, std::fstream::app);
  if(finished) {
    myfile << "Worker " << wid << " finished." << endl;
  } else {
    for(int i : related_cluster) myfile << " " << metadata.sample_filter.ID[i/2];
    myfile << endl;
  }
  myfile.close();
}

inline double Knockoffs::mutate(const int h, const double mut_rate) const {
  if(h==1) return(1.0-mut_rate);
  else return(mut_rate);
}

void Knockoffs::CV(const int num_threads) {
  CrossVal cv(H, references_global, alpha[0], b, phys_dist, seed);
  pair<double,double> cv_best = cv.fit(num_threads);
  init_hmm(K, cv_best.first, cv_best.second);
}

#endif
