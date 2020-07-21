#ifndef _CROSS_VALIDATION_CPP
#define _CROSS_VALIDATION_CPP

#include "cross_validation.h"

CrossVal::CrossVal(const vector <chaplotype> & _H, const vector< vector<int> > & _ref,
                   const vector< vector<double> > & _alpha, const vector <double> & _b,
                   const vector <double> & _phys_dist, int _seed) :
  H(_H), ref(_ref), alpha(_alpha), b(_b), phys_dist(_phys_dist) {

  num_haps = H.size();
  num_snps = H[0].size();
  K = alpha[0].size();
  seed = _seed;
}

pair<double,double> CrossVal::fit(int num_threads) {
  // Parameters
  int num_rep = std::min(1000,num_haps);
  int num_it = 3;
  int num_masks = 10;
  int m_block_number = 5;
  int m_block_length = 100;
  int grid_size = 10;
  bool show_progress = true;

  vector<double> rho_grid, lambda_grid;
  double lambda_max = 1e-2;
  double rho_left=1e1, rho_right=1e-4, lambda_left=1e-6, lambda_right=lambda_max;
  log_range(rho_grid, rho_left, rho_right, grid_size);
  log_range(lambda_grid, lambda_left, lambda_right, grid_size);
  boost::random::taus88 rng;

  // Initialize best parameters
  double rho_best = rho_grid[0];
  double lambda_best = lambda_max;

  // Initialize masks
  vector< vector<bool> > mask;
  vector< vector< vector<int> > > masked_blocks;
  for(int idx_mask=0; idx_mask<num_masks; idx_mask++) {
    vector<bool> one_mask;
    vector< vector<int> > one_masked_blocks;
    double q_start = (double)(m_block_number)/(double)(num_snps);
    double q_stop = 1.0/(double)m_block_length;
    init_mask(one_mask, one_masked_blocks, q_start, q_stop, rng);
    mask.push_back(one_mask);
    masked_blocks.push_back(one_masked_blocks);
  }

  // Initialize progress bar
  int progress_nsteps = 100;
  int progress_period = std::max(1, (int)(num_it*num_masks*grid_size/progress_nsteps));
  int progress_step = 0;
  if(show_progress) {
    cout << "|" << flush;
    for(int s=0; s<progress_nsteps; s++) cout << ".";
    cout << "|" << endl;
    cout << "|" << flush;
  }

  // Print progress header
  if(!show_progress) {
    cout << setw(6) << "Iter" << setw(12) << "Error" << setw(12) << "Rho" << setw(12) << "Lambda" << endl;
  }

  for(int it=0; it<num_it; it++) {
    string optimizing;
    if(it%2==0) optimizing="rho";
    else optimizing="lambda";

    vector<double> cv_mean(grid_size,0);
    vector<double> cv_sd(grid_size,0);
    vector<double> cv_var(grid_size,0);

    for(int idx_mask=0; idx_mask<num_masks; idx_mask++) {
      for(int idx=0; idx<grid_size; idx++) {

        vector <double> b_masked;
        vector <double> mut_rates;

        // Apply mask and initialize HMM
        if(optimizing=="rho") {
          apply_mask(mask[idx_mask], rho_grid[idx], lambda_best, b_masked, mut_rates);
        } else {
          apply_mask(mask[idx_mask], rho_best, lambda_grid[idx], b_masked, mut_rates);
        }

        // Compute CV errors with this mask
        vector<double> errors(num_rep,0);
        CV_errors(errors, mask[idx_mask], masked_blocks[idx_mask], b_masked, mut_rates, num_threads, rng);

        // Store mean and SD of errors
        cv_mean[idx] += compute_mean(errors) / (double)(num_masks);
        cv_var[idx] += (std::pow(compute_sd(errors),2) / (double)(errors.size())) / (double)(num_masks);

        // Update progress
        if(show_progress && (progress_step % progress_period == 0)) {
          cout << "=" << flush;
        }
        progress_step++;
      }
    }
    for(int idx=0; idx<grid_size; idx++) {
      cv_sd[idx] = std::sqrt(cv_var[idx]);
    }

    // Find 1-SD best parameters
    int new_idx_best = min_1sd(cv_mean, cv_sd, true);

    // Find range of promising parameters
    std::pair<int, int> new_range = range_1sd(cv_mean, cv_sd);

    if(!show_progress) {
      for(int idx=0; idx<grid_size; idx++) {
        if(optimizing=="rho") {
          cout << "Parameters (" << setw(10) << setprecision(4) << rho_grid[idx];
          cout << ", " << setw(10) << setprecision(4) << lambda_best << ")";
        } else {
          cout << "Parameters (" << setw(10) << setprecision(4) << rho_best;
          cout << ", "<< setw(10) << setprecision(4)  << lambda_grid[idx] << ")";
        }
        cout << ", CV error: " << setw(6) << setprecision(4) << 100*cv_mean[idx];
        cout << " +- " << setw(6) << setprecision(4) << 100 * 2.0 * cv_sd[idx] << endl;
      }
    }

    if(optimizing=="rho") {
      rho_best = rho_grid[new_idx_best];
      rho_left = 0.5*rho_grid[new_range.first];
      rho_right = boost::algorithm::clamp(2*rho_grid[new_range.second], 1e-12, 1e-2);
    } else {
      lambda_best = lambda_grid[new_idx_best];
      lambda_left = 0.5*lambda_grid[new_range.first];
      lambda_right = boost::algorithm::clamp(2*lambda_grid[new_range.second], 1e-12, lambda_max);
    }

    // Print progress
    if(!show_progress) {
      cout << setw(6) << setprecision(3) << it << setw(12) << setprecision(3) << 100 * cv_mean[new_idx_best];
      cout << setw(12) << setprecision(3) << rho_best << setw(12) << setprecision(3) << lambda_best << endl;
    }

    // Update grids
    log_range(rho_grid, rho_left, rho_right, grid_size);
    log_range(lambda_grid, lambda_left, lambda_right, grid_size);
  }

  // Finalize progress bar
  if(show_progress) {
    for(int s=progress_step; s<progress_nsteps; s++) {
      cout << "=" << flush;
    }
    cout << "|" << endl;
    cout << "Optimal parameters found: rho = " << rho_best << ", lambda = " << lambda_best << ".";
    cout << endl << endl;
  }

  // Return optimal parameters
  return(make_pair(rho_best, lambda_best));
}

void CrossVal::CV_errors(vector<double> & errors,
                         const vector<bool> & mask, const vector< vector<int> > & masked_blocks,
                         const vector <double> & b_masked, const vector <double> & mut_rates,
                         int num_threads, boost::random::taus88 & rng) const {

  int num_rep = errors.size();

  // Randomly assign individuals to workers
  vector<int> i_list_full(num_haps);
  for(int i=0; i<num_haps; i++) i_list_full[i] = i;
  std::random_shuffle(i_list_full.begin(), i_list_full.end());

  vector< vector<double> > errors_w(num_threads);
  vector< vector<int> > i_w;
  for(int w=0; w<num_threads; w++ ) {
    int i_start = w*num_rep/num_threads;
    int i_end = (w+1)*num_rep/num_threads;
    if(w==(num_threads-1)) i_end = num_rep;
    int job_size = i_end-i_start;
    vector<int> vec(job_size, 0);
    for(int i=i_start; i<i_end; i++) vec[i-i_start] = i_list_full[i];
    i_w.push_back(vec);
    errors_w[w].resize(job_size,0);
  }

  // Random number generators for each worker
  vector<boost::random::taus88> rng_w(num_threads);

  // Create workers
  vector<boost::thread> workers;
  for(int w=0; w<num_threads; w++ ) {
    //cout << "Worker " << w << " : size " << i_w[w].size() << endl;

    // Seed the random number generator
    rng_w[w].seed(seed + w);

    // Assign task to worker
    workers.push_back(boost::thread(&CrossVal::CV_errors_worker, this, boost::ref(i_w[w]),
                                    boost::ref(errors_w[w]),
                                    boost::ref(mask), boost::ref(masked_blocks),
                                    boost::ref(b_masked), boost::ref(mut_rates),
                                    boost::ref(rng_w[w])));
  }

  // Launch workers
  for(int w=0; w<num_threads; w++ ) {
    workers[w].join();
  }

  // Combine results
  int idx=0;
  for(int w=0; w<num_threads; w++) {
    for(int i=0; i<errors_w[w].size(); i++) {
      errors[idx] = errors_w[w][i];
      idx++;
    }
  }
}

void CrossVal::CV_errors_worker(const vector<int> & i_list, vector<double> & errors,
                                const vector<bool> & mask, const vector< vector<int> > & masked_blocks,
                                const vector <double> & b_masked, const vector <double> & mut_rates,
                                boost::random::taus88 & rng) const {
  bool verbose = false;

  for(int i_idx=0; i_idx<i_list.size(); i_idx++) {
    int i = i_list[i_idx];

    // Sample posterior
    vector<int> z(num_snps, -1);
    //sample_posterior_MC(i, z, mask, rng);
    sample_posterior_MC(i, z, mask, b_masked, mut_rates, rng);

    // Impute MC
    impute_MC(i, z, mask, masked_blocks, rng);

    // Impute HMM
    vector<int> x(num_snps, 0);
    impute_HMM(i, x, z, mask, mut_rates, rng);

    // Count imputation errors
    double n_errors = 0;
    double num_missing = 0;
    for(int j=0; j<num_snps; j++) {
      if(mask[j]) {
        n_errors += (H[i][j] != x[j]);
        num_missing++;
      }
    }
    if(num_missing>0) {
      errors[i_idx] = n_errors/(double)(num_missing);
    }

    // Debug
    if(verbose) {
      std::cout<<"True and imputed values:" << endl;
      for(int j=0; j<num_snps; j++) {
        cout << j << " " << mask[j] << " ";
        if(mask[j]) cout << "NA ";
        else cout << z[j] << " ";
        cout<<z[j]<<" ";
        int k = z[j];
        cout<<H[i][j]<<" ";
        cout<<H[ref[i][k]][j]<<" ";
        cout<<x[j]<<" ";
        if(!mask[j] && (H[i][j] != H[ref[i][k]][j])) cout << "*";
        if(mask[j] && (H[i][j] != x[j])) cout << "error!";
        if(mask[j] && (H[i][j] == x[j])) cout << "ok";
        cout << endl;
      }
    }
  }
}

void CrossVal::impute_HMM(int i, vector<int> & x, const vector<int> &z, const vector<bool> & mask,
                          const vector <double> & mut_rates, boost::random::taus88 & rng) const {
  vector<double> weights_hmm(2);
  for(int j=0; j<num_snps; j++) {
    int k = z[j];
    if(mask[j]) {
      double theta_j = mutate(H[ref[i][k]][j], mut_rates[j]);
      weights_hmm[0] = 1.0 - theta_j;
      weights_hmm[1] = theta_j;
      x[j] = weighted_choice(weights_hmm, rng);
    } else {
      x[j] = H[i][j];
    }
  }
}

void CrossVal::impute_MC(int i, vector<int> & z,
                         const vector<bool> & mask, const vector< vector<int> > & masked_blocks,
                         boost::random::taus88 & rng) const {

  // FIXME: assuming extremities not masked

  for(int s=0; s<masked_blocks.size(); s++) {
    int m = masked_blocks[s].size();
    // cout << "Imputing block " << s << ": ";
    // for(int l=0; l<m; l++) {
    //   cout << masked_blocks[s][l] << " ";
    // }
    // cout << endl;
    vector< vector<double> > a_bar(m, vector<double>(K) );
    vector<double> b_bar(m);
    int j1 = masked_blocks[s].back()+1;
    int z1 = z[j1];
    // Initialize last element of a_bar, b_bar
    b_bar[m-1] = b[j1];
    for(int k=0; k<K; k++) {
      a_bar[m-1][k] = (1.0-b[j1])*alpha[i][k];
    }
    // Backward-compute the other elements of a_bar, b_bar
    for(int l=m-2; l>=0; l--) {
      int j = masked_blocks[s][l];
      b_bar[l] = b[j]*b_bar[l+1];
      double a_const = b_bar[l+1] * (1.0-b[j])*alpha[i][z1];
      for(int k=0; k<K; k++) {
        a_const += (1.0-b[j])*alpha[i][k] * a_bar[l+1][z1];
      }
      for(int k=0; k<K; k++) {
        a_bar[l][k] = a_const + b[j]*a_bar[l+1][k];
      }
    }

    // Forward-sample missing elements of Z
    vector<double> weights_mc(K);
    for(int l=0; l<m; l++) {
      int j = masked_blocks[s][l];
      for(int k=0; k<K; k++) {
        double w_left = (1.0-b[j])*alpha[i][k] + b[j] * (double)(k==z[j-1]);
        double w_right = a_bar[l][k] + (double)(k==z1);
        weights_mc[k] = w_left * w_right;
      }
      z[j] = weighted_choice(weights_mc, rng);
    }
  }
}

void CrossVal::sample_posterior_MC(int i, vector<int> & z, const vector<bool> & mask,
                                   const vector <double> & b_masked, const vector <double> & mut_rates,
                                   boost::random::taus88 & rng) const {

  vector< vector<double> > backward(num_snps, vector<double> (K));
  vector<double> weights_mc(K);

  // Compute backward weights
  std::fill(backward[num_snps-1].begin(), backward[num_snps-1].end(), 1.0);
  for(int j=num_snps-2; j>=0; j--) {
    if(mask[j]) {
      for(int k=0; k<K; k++) {
        backward[j][k] = backward[j+1][k];
      }
    } else {
      double backward_const = 0.0;
      for(int k=0; k<K; k++) {
        double a_j1_k = (1.0-b_masked[j+1]) * alpha[i][k];                     // a[j+1][k]
        double theta_j1_k = mutate(H[ref[i][k]][j+1], mut_rates[j+1]);         // theta[j+1][k]
        if ( H[i][j+1] == 1) {
          backward_const += a_j1_k * theta_j1_k * backward[j+1][k];
        }
        else {
          backward_const += a_j1_k * (1.0-theta_j1_k) * backward[j+1][k];
        }
      }
      double backwardSum = 0.0;
      for(int k=0; k<K; k++) {
        double theta_j1_k = mutate(H[ref[i][k]][j+1], mut_rates[j+1]);        // theta[j+1][k]
        if ( H[i][j+1] == 1) {
          backward[j][k] = backward_const + b_masked[j+1] * theta_j1_k * backward[j+1][k];
        }
        else {
          backward[j][k] = backward_const + b_masked[j+1] * (1.0-theta_j1_k) * backward[j+1][k];
        }
        backwardSum += backward[j][k];
      }
      for(int k=0; k<K; k++) {
        backward[j][k] /= backwardSum;
      }
    }
  }

  // Forward sampling
  if(mask[0]) {
    z[0] = -1;
  } else {
    for(int k=0; k<K; k++) {
      double a_0_k = (1.0-b_masked[0]) * alpha[i][k];                         // a[0][k]
      double theta_0_k = mutate(H[ref[i][k]][0], mut_rates[0]);               // theta[0][k]
      if ( H[i][0] == 1) {
        weights_mc[k] = a_0_k * theta_0_k * backward[0][k];
      }
      else {
        weights_mc[k] = a_0_k * (1.0-theta_0_k) * backward[0][k];
      }
    }
    z[0] = weighted_choice(weights_mc, rng);
  }

  for(int j=1; j<num_snps; j++) {
    if(mask[j]) {
      z[j] = z[j-1];
    } else {
      for(int k=0; k<K; k++) {
        double a_j_k = (1.0-b_masked[j]) * alpha[i][k];                       // a[j][k]
        double theta_j_k = mutate(H[ref[i][k]][j], mut_rates[j]);             // theta[j][k]
        if ( H[i][j] == 1) {
          weights_mc[k] = (a_j_k + b_masked[j]*(double)(k==z[j-1])) * theta_j_k * backward[j][k];
        }
        else {
          weights_mc[k] = (a_j_k + b_masked[j]*(double)(k==z[j-1])) * (1.0-theta_j_k) * backward[j][k];
        }
      }
      z[j] = weighted_choice(weights_mc, rng);
    }
  }
}

void CrossVal::init_mask(vector<bool> & mask, vector< vector<int> > & masked_blocks,
                         double q_start, double q_stop, boost::random::taus88 & rng) const {

  //cout << "[DEBUG] in init_mask(" << q_start << ", " << q_stop << ")" << endl;

  int tot_masked = 0;
  while(tot_masked==0) {
    bool masking = false;
    mask.clear();
    mask.resize(num_snps);
    masked_blocks.clear();
    vector<int> block;
    for(int j=2; j<(num_snps-2); j++) {
      mask[j] = masking;
      tot_masked += (int)(masking);
      if(masking) {
        block.push_back(j);
        if(runif(rng)<q_stop) {
          masking = false;
          masked_blocks.push_back(block);
          block.clear();
        }
      } else {
        if(runif(rng)<q_start) {
          masking = true;
        }
      }
    }
    if(masking && (block.size()>0)) {
      masked_blocks.push_back(block);
    }
  }
  // for(int l=0; l<masked_blocks.size(); l++)
  //   cout << "Block " << l <<" has size " << masked_blocks[l].size() << endl;
}

void CrossVal::apply_mask(const vector<bool> & mask, double rho, double lambda,
                          vector <double> & b_masked, vector <double> & mut_rates) {
  // cout << "Masked SNPs: ";
  // for(int j=0; j<num_snps; j++) {
  //   if(mask[j]) cout << j << " ";
  // }
  // cout << endl;

  // Initialize recombination rates
  b_masked.clear();
  b_masked.resize(num_snps);
  vector<double> rec_rates(num_snps, 1);
  vector<double> rec_rates_masked(num_snps);
  double remainder = 0;
  for(unsigned int j=0; j<num_snps; j++) {
    if(j>0) {
      rec_rates[j] = rho * (double) phys_dist[j];
    }
    if(mask[j]) {
      remainder += rec_rates[j];
      rec_rates_masked[j] = 0;
    } else {
      rec_rates_masked[j] = remainder + rec_rates[j];
      remainder = 0;
    }
    if(j>0) b_masked[j] = std::exp(-rec_rates_masked[j]);
    else b_masked[j] = 0;
  }

  // Initialize mutation rates
  mut_rates.clear();
  mut_rates.resize(num_snps);
  for(unsigned int j=0; j<num_snps; j++) {
    mut_rates[j] = lambda;
  }

  // cout << "ID" << setw(15) << "Masked" << setw(15) << "Rec" << setw(15) << "Rec_m";
  // cout << setw(15) << "b" << setw(15) << "b_m" << endl;
  // for(int j=0; j<num_snps; j++) {
  //   cout << j << "\t" << mask[j] << setw(15) << rec_rates[j] << setw(15) << rec_rates_masked[j] << setw(15) << b[j] << setw(15) << b_masked[j] << endl;
  // }

}

void log_range(vector<double> & values, double min_value, double max_value, int length) {
  values.resize(length);
  const double min_log = std::log(min_value);
  const double max_log = std::log(max_value);
  const double log_increment = ( max_log - min_log ) / length;
  double log_value = min_log;
  for(int i=0 ; i < length; i++) {
    log_value += log_increment;
    values[i] = std::exp(log_value);
  }
}

double compute_sd(const vector<double> & data) {
  int n = data.size();
  double mean = compute_mean(data);
  double sd = 0;
  for(int i=0; i<n; i++) sd += std::pow(data[i] - mean, 2);
  return std::sqrt(sd / (double)n);
}

double compute_mean(const vector<double> & data) {
  int n = data.size();
  double sum = 0;
  for(int i=0; i<n; i++) sum += data[i];
  return(sum/(double)n);
}

int min_1sd(const vector<double> & mean, const vector<double> & sd, bool prefer_larger) {
  int idx_min = std::distance(mean.begin(), std::min_element(mean.begin(), mean.end()));
  int idx_best = idx_min;

  vector<double> mean_plus = mean;
  vector<double> mean_minus = mean;
  for(int i=0; i<mean.size(); i++) {
    mean_plus[i] += sd[i];
    mean_minus[i] -= sd[i];
  }
  int idx_min_max = std::distance(mean_plus.begin(), std::min_element(mean_plus.begin(), mean_plus.end()));
  double val_min_max = mean_plus[idx_min_max];
  if(prefer_larger) {
    for(int i=idx_min; i<mean.size(); i++) {
      if((mean_minus[i]<=val_min_max)) idx_best = i;
    }
  } else {
    for(int i=idx_min; i>=0; i--) {
      if((mean_minus[i]<=val_min_max)) idx_best = i;
   }
  }
  return(idx_best);
}

pair<int,int> range_1sd(const vector<double> & mean, const vector<double> & sd) {
  vector<double> mean_plus = mean;
  vector<double> mean_minus = mean;
  for(int i=0; i<mean.size(); i++) {
    mean_plus[i] += sd[i];
    mean_minus[i] -= sd[i];
  }

  int idx_min_max = std::distance(mean_plus.begin(), std::min_element(mean_plus.begin(), mean_plus.end()));
  double val_min_max = mean_plus[idx_min_max];

  int idx_left = idx_min_max;
  int idx_right = idx_min_max;

  for(int i=idx_min_max; i<mean.size(); i++) {
    if(mean_minus[i]<val_min_max) idx_right = i;
  }
  for(int i=idx_min_max; i>0; i--) {
    if(mean_minus[i]<val_min_max) idx_left = i;
  }
  return(std::make_pair(idx_left, idx_right));
}

double CrossVal::mutate(int h, double mut_rate) const {
  if(h==1) return(1.0-mut_rate);
  else return(mut_rate);
}

double init_mut_rate(int K) {
  double theta_inv = 0.0;
  for(unsigned int k=1; k<K; k++) {
    theta_inv += 1.0/(double)k;
  }
  double theta = 1.0/theta_inv;
  return(0.5*theta/(theta+(double)K));
}

#endif
