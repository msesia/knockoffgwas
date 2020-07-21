#ifndef _FAMILY_BP_CPP
#define _FAMILY_BP_CPP

#define DEBUG 0

#include "family_bp.h"

FamilyBP::FamilyBP(const vector<int>& _indices, const IbdCluster& ibd_cluster,
                   const vector <chaplotype>& _H,
                   const vector< vector<int > >& _ref,
                   const vector< vector<double> >& _alpha,
                   const vector <double>& _mut_rates, const vector<double>& _b,
                   const Metadata& _metadata, boost::random::taus88& _rng) :
  K(_alpha[0].size()), num_haps(_indices.size()), num_snps(_b.size()),
  indices(_indices), b(_b), metadata(_metadata),
  H(_H), ref(_ref), alpha(_alpha), mut_rates(_mut_rates), rng(_rng) {

  if(DEBUG) cout << endl << "[DEBUG] in FamilyBP() constructor" << endl;

  if(DEBUG) {
    cout << "Number of variants: " << num_snps << endl;
    cout << "Family size: " << num_haps << endl;
    for(int i : indices) {
      cout << " " << i;
    }
    cout << endl;
    // Show IBD cluster
    ibd_cluster.print();
  }

  // Initialize list of neighbors
  for(int i=0; i<num_haps; i++) {
    vector< vector<int> > tmp_v;
    vector<int> tmp_empty;
    for(int j=0; j<num_snps; j++) {
      tmp_v.emplace_back(tmp_empty);
    }
    neighbors.push_back(tmp_v);
  }

  // Fill list of neighbors
  for(int s=0; s<ibd_cluster.size(); s++) {
    IbdSeg segment;
    ibd_cluster.get(s,segment);
    vector<int> segment_indices;
    for(int id : segment.indices) {
      // Find id in list of ids for this family
      auto it = std::find(indices.begin(), indices.end(), id);
      int i = std::distance(indices.begin(), it);
      segment_indices.push_back(i);
    }
    for(int i1 : segment_indices) {
      for(int j=segment.j_min; j<=segment.j_max; j++) {
        for(int i2 : segment_indices) {
          if(i1 != i2) {
            neighbors[i1][j].push_back(i2);
          }
        }
      }
    }
  }

  // // DEBUG: print list of neighbors
  // for(int i1=0; i1<num_haps; i1++) {
  //   cout << "Neighbors for haplotype " << indices[i1] << ":" << endl;
  //   for(int j=0; j<num_snps; j++) {
  //     if(neighbors[i1][j].size()>0) {
  //       cout << j << " ";
  //       for(auto i2 : neighbors[i1][j]) {
  //         cout << indices[i2] << ", ";
  //       }
  //       cout << endl;
  //     }
  //   }
  // }

}

int FamilyBP::run(vector< vector<int> >& Z) {
  // Initialize output container
  Z.resize(num_haps, vector<int>(num_snps,-1));

  // Initialize messages
  init_messages();

  // Run BP until convergence
  if(DEBUG) cout << "Running BP... " << flush;
  int it = run_bp(false);
  // Check convergence
  if(it>=0) {
    if(DEBUG) cout << "converged after " << it << " iterations." << endl;
  } else {
    if(DEBUG) cout << "Warning: BP did not converge!" << endl;
  }

  vector<double> marginals(K,1);
  vector<double> f_j(K,1);

  // Sample from posterior
  sample_posterior(Z);

  // // DEBUG: Make sure Z in consistent on IBD segments
  // for(int i=0; i<num_haps; i++) {
  //   for(int j=0; j<num_snps; j++) {
  //     for(int i2 : neighbors[i][j]) {
  //       if(Z[i][j]!=Z[i2][j]) {
  //         cout << "Error: inconsistency found on Z (" << i << ", " << j << ", " << i2 << ")" << endl;
  //         throw_error_bug("Inconsistency found on Z");
  //       }
  //     }
  //   }
  // }

  // // DEBUG: Make sure H in consistent on IBD segments
  // for(int i=0; i<num_haps; i++) {
  //   int i_abs = indices[i];
  //   double agreement_num = 0;
  //   double agreement_den = 0;
  //   for(int j=0; j<num_snps; j++) {
  //     for(int i2 : neighbors[i][j]) {
  //       int i2_abs = indices[i2];
  //       if(H[i_abs][j]==H[i2_abs][j]) {
  //         agreement_num++;
  //       }
  //       agreement_den++;
  //     }
  //   }
  //   if(agreement_den==0) agreement_den = -1;
  //   const string i_name = metadata.sample_filter.ID[i_abs/2];
  //   cout << "[DEBUG] IBD agreement for haplotype " << i << "(" << i_name << "): ";
  //   cout  << agreement_num/agreement_den << endl;
  // }

  //exit(0);

  // // DEBUG: Make sure references are consistent
  // for(int i=0; i<num_haps; i++) {
  //   int i_abs = indices[i];
  //   double agreement_num = 0;
  //   double agreement_den = 0;
  //   for(int j=0; j<num_snps; j++) {
  //     for(int i2 : neighbors[i][j]) {
  //       int i2_abs = indices[i2];
  //       for(int k=0; k<K; k++) {
  //         if(H[ref[i_abs][k]][j]==H[ref[i2_abs][k]][j]) {
  //           agreement_num++;
  //         }
  //         agreement_den++;
  //       }
  //     }
  //   }
  //   if(agreement_den==0) agreement_den = -1;
  //   cout << "[DEBUG] reference agreement for haplotype " << i << "(" << i_abs << "): ";
  //   cout  << agreement_num/agreement_den << endl;
  // }

  return(0);
}

void FamilyBP::sample_posterior(vector< vector<int> >& Z) {
  if(DEBUG) cout << "Sampling junction nodes from posterior... " << endl;
  // First sample the junction nodes
  for(int i=0; i<num_haps; i++) {
    //cout << "[DEBUG] posterior sampling for i=" << i <<":" << endl;
    for(int j=0; j<num_snps; j++) {
      if(is_inactive(i,j)) continue;
      // Check whether this may be a junction
      bool is_junction = false;
      if(j>0) {
        if(neighbors[i][j] != neighbors[i][j-1]) is_junction = true;
      }
      if(j<(num_snps-1)) {
        if(neighbors[i][j] != neighbors[i][j+1]) is_junction = true;
      }
      // Verify whether this is really a junction
      vector<int> neigh_diff_left, neigh_diff_right;
      if(is_junction) {
        if(j>0) {
          if(neighbors[i][j] != neighbors[i][j-1]) {
            std::set_difference(neighbors[i][j].begin(), neighbors[i][j].end(),
                                neighbors[i][j-1].begin(), neighbors[i][j-1].end(),
                                std::inserter(neigh_diff_left, neigh_diff_left.end()));
          }
        }
        if(j<(num_snps-1)) {
          if(neighbors[i][j] != neighbors[i][j+1]) {
            std::set_difference(neighbors[i][j].begin(), neighbors[i][j].end(),
                                neighbors[i][j+1].begin(), neighbors[i][j+1].end(),
                                std::inserter(neigh_diff_right, neigh_diff_right.end()));
          }
        }
        if((neigh_diff_left.size()==0) && (neigh_diff_right.size()==0)) is_junction = false;
      }
      if(is_junction) {
        // Sample from posterior
        const int k = sample_posterior_joint(i,j,neigh_diff_left,neigh_diff_right);
        Z[i][j] = k;
        condition_on(i, j, k);
        for(int i2 : neighbors[i][j]) {
          Z[i2][j] = k;
          condition_on(i2, j, k);
        }
        if(DEBUG) cout << "Running BP again (i=" << i << ", j=" << j << ")... " << flush;
        int it = run_bp(false);
        // Check convergence
        if(it>=0) {
          if(DEBUG) cout << "converged after " << it << " iterations." << endl;
        } else {
          if(DEBUG) cout << "Warning: BP did not converge!" << endl;
        }
      }

    }
  }
  if(DEBUG) cout << "Sampling interior nodes from posterior... " << endl;
  // Then sample all other variables
  for(int i=0; i<num_haps; i++) {
    for(int j=0; j<num_snps; j++) {
      if(is_inactive(i,j)) continue;
      vector<int> neigh_diff_left, neigh_diff_right;
      const int k = sample_posterior_joint(i, j, neigh_diff_left, neigh_diff_right);
      Z[i][j] = k;
      condition_on(i, j, k, true);
      for(int i2 : neighbors[i][j]) {
        Z[i2][j] = k;
        condition_on(i2, j, k, true);
      }
    }
  }
}

void FamilyBP::condition_on(const int i, const int j, const int k) {
  condition_on(i, j, k, false);
}

void FamilyBP::condition_on(const int i, const int j, const int k, const bool forward) {
  int i_abs = indices[i];

  for(int l=0; l<K; l++) {
    // Set rightward message
    if(j<(num_snps-1)) {
      double a_jl = (1.0-b[j+1]) * alpha[i_abs][l];
      messages_right[i][j][l] = a_jl + b[j+1] * (double)(l==k);
    }
    // Set leftward message (not needed with forward sampling)
    if(!forward) {
      if(j>0) {
        double a_jk = (1.0-b[j]) * alpha[i_abs][k];
        messages_left[i][j][l] = a_jk + b[j] * (double)(l==k);
      }
    }
  }
  // Make node inactive
  inactive_nodes.insert(std::make_pair(i,j));
}

// void FamilyBP::get_marginals(const int i, const int j, vector<double>& f_j, vector<double>& marginals) const {
//   // Find absolute haplotype id
//   const int i_abs = indices[i];

//   // Find out which neighbors have just been added
//   vector<int> neighbors_diff_left;
//   if(j>0) {
//     if(neighbors[i][j] != neighbors[i][j-1]) {
//       std::set_difference(neighbors[i][j].begin(), neighbors[i][j].end(),
//                           neighbors[i][j-1].begin(), neighbors[i][j-1].end(),
//                           std::inserter(neighbors_diff_left, neighbors_diff_left.end()));
//     }
//   }
//   // Find out which neighbors are about to be removed
//   vector<int> neighbors_diff_right;
//   if(j<(num_snps-1)) {
//     if(neighbors[i][j] != neighbors[i][j+1]) {
//       std::set_difference(neighbors[i][j].begin(), neighbors[i][j].end(),
//                           neighbors[i][j+1].begin(), neighbors[i][j+1].end(),
//                           std::inserter(neighbors_diff_right, neighbors_diff_right.end()));
//     }
//   }

//   // Pre-compute emissions probabilities
//   const int h_real = H[i_abs][j];
//   for(int k=0; k<K; k++) {
//     f_j[k] = emission_prob(h_real, H[ref[i_abs][k]][j], mut_rates[j]);
//   }
//   for(int i2 : neighbors[i][j]) {
//     const int i2_abs = indices[i2];
//     const int h2_real = H[i2_abs][j];
//     for(int k=0; k<K; k++) {
//       f_j[k] *= emission_prob(h2_real, H[ref[i_abs][k]][j], mut_rates[j]);
//     }
//   }
//   // Standardize emission probabilities
//   if(neighbors[i][j].size() > 0) {
//     const double num_neighbors = neighbors[i][j].size() + 1;
//     for(int k=0; k<K; k++) {
//       f_j[k] = std::pow(f_j[k], 1.0/num_neighbors);
//     }
//   }

//   for(int k=0; k<K; k++) {
//     marginals[k] = f_j[k];
//     if(j>0) marginals[k] *= messages_right[i][j-1][k];
//     if(j<(num_snps-1)) marginals[k] *= messages_left[i][j+1][k];
//   }
//   for(int i2 : neighbors_diff_left) {
//     for(int k=0; k<K; k++) {
//       marginals[k] *= messages_right[i2][j-1][k];
//     }
//   }
//   for(int i2 : neighbors_diff_right) {
//     for(int k=0; k<K; k++) {
//       marginals[k] *= messages_left[i2][j+1][k];
//     }
//   }
// }


int FamilyBP::run_bp(bool verbose) {
  const int it_max = 20;

  // Initialize hap list
  vector<int> hap_idx(num_haps);
  std::iota(hap_idx.begin(), hap_idx.end(), 0);

  // Initialize snp list (both directions)
  vector<int> snp_idx_forw(num_snps);
  std::iota(snp_idx_forw.begin(), snp_idx_forw.end(), 0);
  vector<int> snp_idx_back = snp_idx_forw;
  std::reverse(snp_idx_back.begin(), snp_idx_back.end());
  vector<int> snp_idx;

  // Initialize direction list
  const vector<string> directions = {"right", "left"};

  bool converged = false;
  int it = 0;
  while( (!converged)&& (it++ < it_max) ) {
    if(verbose) cout << " -- Iteration " << it << " -- " << endl;
    double diff = 0;

    for(string dir : directions) {
      if(dir == "right") {
        snp_idx = snp_idx_forw;
      }
      if(dir == "left") {
        snp_idx = snp_idx_back;
      }

      vector<double> message(K,1.0); // Store new messages
      vector<double> f_j(K,1.0); // Store emission probabilities

      double diff_tmp = 0;
      for(int j : snp_idx) {
        for(int i : hap_idx) {
          if (is_inactive(i,j)) continue; // Omit messages for inactive nodes

          double is_junction = false;
          if(j<(num_snps-1)) {
            if(neighbors[i][j] != neighbors[i][j+1]) is_junction = true;
          }
          if(j>0) {
            if(neighbors[i][j] != neighbors[i][j-1]) is_junction = true;
          }

          // Find out which neighbors are changing
          vector<int> neigh_diff_left;
          vector<int> neigh_diff_right;
          if(is_junction) {
            if(j>0) {
              if(neighbors[i][j] != neighbors[i][j-1]) {
                std::set_difference(neighbors[i][j].begin(), neighbors[i][j].end(),
                                    neighbors[i][j-1].begin(), neighbors[i][j-1].end(),
                                    std::inserter(neigh_diff_left, neigh_diff_left.end()));
              }
            }
            if(j<(num_snps-1)) {
              if(neighbors[i][j] != neighbors[i][j+1]) {
                std::set_difference(neighbors[i][j].begin(), neighbors[i][j].end(),
                                    neighbors[i][j+1].begin(), neighbors[i][j+1].end(),
                                    std::inserter(neigh_diff_right, neigh_diff_right.end()));
              }
            }
            // if((neigh_diff_left.size()>0) | (neigh_diff_right.size()>0)) {
            //   // DEBUG: print neighbor differences
            //   cout << " (i="  << i <<", j=" << j << "); neigh_diff_left:";
            //   for(int i2 : neigh_diff_left) cout << " " << i2;
            //   cout << endl;
            //   cout << " (i="  << i <<", j=" << j << "); neigh_diff_right:";
            //   for(int i2 : neigh_diff_right) cout << " " << i2;
            //   cout << endl;
            // }
          }

          // Send messages
          if(dir == "right") {
            diff_tmp += send_messages_right(i, j, message, f_j, neigh_diff_left, neigh_diff_right);
            for(int k=0; k<K; k++) messages_right[i][j][k] = message[k];
          }
          if(dir == "left") {
            if(j<=0) continue;
            diff_tmp += send_messages_left(i, j, message, f_j, neigh_diff_left, neigh_diff_right);
            for(int k=0; k<K; k++) messages_left[i][j][k] = message[k];
          }
        }
      }
      diff_tmp /= (double)(snp_idx.size() * hap_idx.size());
      if(verbose) cout << "Sent " << dir << " messages: " << diff_tmp << endl;
      diff += diff_tmp;
    }

    if(diff<= 1e-3) {
      converged = true;
      if(verbose) cout << "BP converged after " << it << " iterations!" << endl;
    }
  }

  if(converged) return(it);
  else return(-1);
}

void FamilyBP::init_messages() {
  double p0 = 1.0/((double)K);
  messages_right.clear();
  messages_right.resize(num_haps, vector2d (num_snps, vector<double> (K,p0)));
  messages_left.clear();
  messages_left.resize(num_haps, vector2d (num_snps, vector<double> (K,p0)));
}

double FamilyBP::send_messages_right(const int i1, const int j, vector<double>& message, vector<double>& f_j,
                                     vector<int>& neigh_diff_left, vector<int>& neigh_diff_right) const {
  if(j>=(num_snps-1)) return(0);

  // Find absolute haplotype id
  const int i1_abs = indices[i1];

  // Pre-compute emissions probabilities
  const int h_real = H[i1_abs][j];
  for(int k=0; k<K; k++) {
    f_j[k] = emission_prob(h_real, H[ref[i1_abs][k]][j], mut_rates[j]);
  }
  for(int i2 : neighbors[i1][j]) {
    const int i2_abs = indices[i2];
    const int h2_real = H[i2_abs][j];
    for(int k=0; k<K; k++) {
      f_j[k] *= emission_prob(h2_real, H[ref[i1_abs][k]][j], mut_rates[j]);
    }
  }
  // Standardize emission probabilities
  if(neighbors[i1][j].size() > 0) {
    const double num_neighbors = neighbors[i1][j].size() + 1;
    for(int k=0; k<K; k++) {
      f_j[k] = std::pow(f_j[k], 1.0/num_neighbors);
    }
  }

  // Precompute a_{j+1} and b_{j+1} weights (assuming alpha = 1/K)
  const double b_j = b[j+1];
  const double a_j = (1.0-b_j) / (double)(K);

  // Compute message
  double const_msg = 0;
  for(int l=0; l<K; l++) {
    double new_const_msg = a_j * f_j[l];
    if(j>0) new_const_msg *= messages_right[i1][j-1][l];
    for(int i2 : neigh_diff_left) new_const_msg *= messages_right[i2][j-1][l];
    for(int i2 : neigh_diff_right) new_const_msg *= messages_left[i2][j+1][l];
    const_msg += new_const_msg;
  }
  double message_sum = 0;
  for(int k=0; k<K; k++) {
    double new_msg = b_j * f_j[k];
    if(j>0) new_msg *= messages_right[i1][j-1][k];
    for(int i2 : neigh_diff_left) new_msg *= messages_right[i2][j-1][k];
    for(int i2 : neigh_diff_right) new_msg *= messages_left[i2][j+1][k];
    message[k] = const_msg + new_msg;
    message_sum += message[k];
  }

  // Normalize message
  double diff = 0;
  const double message_sum_inv = 1.0 / message_sum;
  for(int k=0; k<message.size(); k++) {
    message[k] = message[k] * message_sum_inv;
    diff += std::abs(message[k]-messages_right[i1][j][k]);
  }

  // // DEBUG: print message
  // if((i1%2==0)&&(metadata.sample_filter.ID[indices[i1]/2]=="1059895")) {
  //   if(j<10) {
  //     cout << "[DEBUG]: in send_messages_right(" << i1 << ", " << j << "):" << endl;
  //     for(int k=0; k<message.size(); k++) {
  //       cout << " " << message[k];
  //     }
  //     cout << endl;
  //   }
  // }

  return(diff);
}

double FamilyBP::send_messages_left(const int i1, const int j, vector<double>& message, vector<double>& f_j,
                                    vector<int>& neigh_diff_left, vector<int>& neigh_diff_right) const {
  if(j<=0) return(0);

  // Find absolute haplotype id
  const int i1_abs = indices[i1];

  // Pre-compute emissions probabilities
  const int h_real = H[i1_abs][j];
  for(int k=0; k<K; k++) {
    f_j[k] = emission_prob(h_real, H[ref[i1_abs][k]][j], mut_rates[j]);
  }
  for(int i2 : neighbors[i1][j]) {
    const int i2_abs = indices[i2];
    const int h2_real = H[i2_abs][j];
    for(int k=0; k<K; k++) {
      f_j[k] *= emission_prob(h2_real, H[ref[i1_abs][k]][j], mut_rates[j]);
    }
  }
  // Standardize emission probabilities
  if(neighbors[i1][j].size() > 0) {
    const double num_neighbors = neighbors[i1][j].size() + 1;
    for(int k=0; k<K; k++) {
      f_j[k] = std::pow(f_j[k], 1.0/num_neighbors);
    }
  }

  // Precompute a_{j} and b_{j} weights (assuming alpha = 1/K)
  const double b_j = b[j];
  const double a_j = (1.0-b_j) / (double)(K);

  // Compute message
  double const_msg = 0;
  for(int l=0; l<K; l++) {
    double new_const_msg = a_j * f_j[l];
    if(j<(num_snps-1)) new_const_msg *= messages_left[i1][j+1][l];
    for(int i2 : neigh_diff_right) new_const_msg *= messages_left[i2][j+1][l];
    for(int i2 : neigh_diff_left) new_const_msg *= messages_right[i2][j-1][l];
    const_msg += new_const_msg;
  }
  double message_sum = 0;
  for(int k=0; k<K; k++) {
    double new_msg = b_j * f_j[k];
    if(j<(num_snps-1)) new_msg *= messages_left[i1][j+1][k];
    for(int i2 : neigh_diff_right) new_msg *= messages_left[i2][j+1][k];
    for(int i2 : neigh_diff_left) new_msg *= messages_right[i2][j-1][k];
    message[k] = const_msg + new_msg;
    message_sum += message[k];
  }

  // Normalize message
  double diff = 0;
  const double message_sum_inv = 1.0 / message_sum;
  for(int k=0; k<message.size(); k++) {
    message[k] = message[k] * message_sum_inv;
    diff += std::abs(message[k]-messages_left[i1][j][k]);
  }
  return(diff);
}

bool FamilyBP::is_inactive(const int i, const int j) const {
  if(inactive_nodes.count(std::make_pair(i,j))==0) {
    return(false);
  } else {
    return(true);
  }
}

int FamilyBP::sample_posterior_joint(const int i, const int j,
                                     const vector<int>& neigh_diff_left,
                                     const vector<int>& neigh_diff_right) const {
  // Find absolute haplotype id
  const int i_abs = indices[i];
  const int h_real = H[i_abs][j];

  vector<double> weights(K,1);
  for(int k=0; k<K; k++) {
    weights[k] = emission_prob(h_real, H[ref[i_abs][k]][j], mut_rates[j]);
    for(int i2 : neighbors[i][j]) {
      const int h2_real = H[indices[i2]][j];
      weights[k] *= emission_prob(h2_real, H[ref[i_abs][k]][j], mut_rates[j]);
    }
    // Standardize emission probabilities
    if(neighbors[i][j].size() > 0) {
      const double num_neighbors = neighbors[i][j].size() + 1;
      weights[k] = std::pow(weights[k], 1.0/num_neighbors);
    }
    if(j>0) {
      weights[k] *= messages_right[i][j-1][k];
      for(int i2 : neigh_diff_left) weights[k] *= messages_right[i2][j-1][k];
    }
    if(j<(num_snps-1)) {
      weights[k] *= messages_left[i][j+1][k];
      for(int i2 : neigh_diff_right) weights[k] *= messages_left[i2][j+1][k];
    }
  }
  // Sample from posterior
  int k_star = weighted_choice(weights, rng);

  return(k_star);
}

void FamilyBP::print_Z(const vector< vector<int> >& Z) const {
  // DEBUG: write Z to text file
  ofstream myfile;
  myfile.open ("tmp/Z.txt");
  for(int i=0; i<num_haps; i++) {
    for(int j=0; j<num_snps; j++) {
      if(j>0) myfile << " ";
      myfile << Z[i][j];
    }
    myfile << endl;
  }
  myfile.close();

  // DEBUG: write H to text file
  myfile.open ("tmp/H.txt");
  for(int i=0; i<num_haps; i++) {
    for(int j=0; j<num_snps; j++) {
      if(j>0) myfile << " ";
      myfile << H[indices[i]][j];
    }
    myfile << endl;
  }
  myfile.close();
}

inline double emission_prob(const int h1, const int h2, const double mut_rate) {
  if(h1==h2) return(1.0-mut_rate);
  else return(mut_rate);
}


#endif
