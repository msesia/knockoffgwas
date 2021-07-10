#ifndef _FAMILY_KNOCKOFFS_CPP
#define _FAMILY_KNOCKOFFS_CPP

#define DEBUG 0
#define MAXINT (unsigned int)(-1)

#include "family_knockoffs.h"

FamilyKnockoffs::FamilyKnockoffs(const vector<unsigned int> & _indices, const IbdCluster& ibd_cluster,
                                 const vector< vector<unsigned int> > & _Z, const vector< vector<double> > & _alpha,
                                 const vector<double> & _b, boost::random::taus88 & _rng) :
  Z(_Z), alpha(_alpha), rng(_rng) {

  if(DEBUG) cout << "[DEBUG] in FamilyKnockoffs() constructor" << endl;

  // Store data
  indices = _indices;
  b = _b;
  K = alpha[0].size();
  num_haps = indices.size();
  num_snps = b.size();

  if(DEBUG) cout << "Family size: " << num_haps << endl;
  if(DEBUG) cout << "Number of variants: " << num_snps << endl;

  // Debug: Show IBD cluster
  if(DEBUG) ibd_cluster.print();

  // Initialize list of neighbors
  for(unsigned int i=0; i<num_haps; i++) {
    vector< vector<unsigned int> > tmp_v;
    vector<unsigned int> tmp_empty;
    for(unsigned int j=0; j<num_snps; j++) {
      tmp_v.push_back(tmp_empty);
    }
    neighbors.push_back(tmp_v);
  }

  // Fill list of neighbors
  for(unsigned int s=0; s<ibd_cluster.size(); s++) {
    IbdSeg segment;
    ibd_cluster.get(s,segment);
    vector<unsigned int> segment_indices;
    for (unsigned int id : segment.indices) {
      // Find id in list of ids for this family
      auto it = std::find(indices.begin(), indices.end(), id);
      unsigned int i = std::distance(indices.begin(), it);
      segment_indices.push_back(i);
    }
    for(unsigned int i1 : segment_indices) {
      for(unsigned int j=segment.j_min; j<=segment.j_max; j++) {
        for(unsigned int i2 : segment_indices) {
          if(i1 != i2) {
            neighbors[i1][j].push_back(i2);
          }
        }
      }
    }
  }

  // Fill list of IBD segments
  for(unsigned int s=0; s<ibd_cluster.size(); s++) {
    IbdSeg segment;
    ibd_cluster.get(s,segment);
    vector<unsigned int> tmp(segment.indices.size());
    std::copy(segment.indices.begin(), segment.indices.end(), tmp.begin());
    // Convert haplotype indices into local format
    for(unsigned int i=0; i<tmp.size(); i++) {
      tmp[i] = std::distance(indices.begin(), std::find(indices.begin(), indices.end(), tmp[i]));
    }
    segments_i.push_back(tmp);
    segments_j_min.push_back(segment.j_min);
    segments_j_max.push_back(segment.j_max);
  }

}

void FamilyKnockoffs::run(const vector< vector<unsigned int> > & elements, const vector<unsigned int> & groups,
                          vector< vector<unsigned int> > & Zk) {
  unsigned int num_groups = elements.size();

  // // Make sure input agree on IBD segments
  // for(unsigned int i=0; i<num_haps; i++) {
  //   for(unsigned int j=0; j<num_snps; j++) {
  //     for(unsigned int i2 : neighbors[i][j]) {
  //       if(Z[i][j]!=Z[i2][j]) {
  //         cout << "Warning: inconsistency found on Z (" << i << ", " << j << ", " << i2 << ")" << endl;
  //         throw_error_bug("Inconsistency found on Z");
  //       }
  //     }
  //   }
  // }

  // Initialize knockoff container
  Zk.resize(num_haps, vector<unsigned int>(num_snps,0));

  // Debug: Find IBD segments
  if(DEBUG) {
    cout << "[DEBUG]: IBD segments:" << endl;
    for(unsigned int s=0; s<segments_i.size(); s++) {
      cout << "(" << segments_j_min[s] << "--" << segments_j_max[s] << ") ";
      cout << "(" << groups[segments_j_min[s]] << "--" << groups[segments_j_max[s]] << "): " ;

      for(unsigned int i : segments_i[s]) {
        cout << " " << i;
      }
      cout << endl;
    }
  }

  // Find boundary groups for each individual
  vector< vector<unsigned int> > boundary_groups;
  vector< vector< vector<unsigned int> > > boundary_neighbors;
  for(unsigned int i=0; i<num_haps; i++) {
    vector<unsigned int> tmp_groups;
    vector< vector<unsigned int> > tmp_neighbors;
    for(unsigned int j=1; j<(num_snps-1); j++) {
      if(neighbors[i][j].size()==0) continue;
      if( (neighbors[i][j] != neighbors[i][j-1]) | (neighbors[i][j] != neighbors[i][j+1])) {
        // Check whether this group has already been detected
        auto it = std::find(tmp_groups.begin(), tmp_groups.end(), groups[j]);
        if(it == tmp_groups.end()) {
          tmp_groups.push_back(groups[j]);
          tmp_neighbors.push_back(neighbors[i][j]);
        } else {
          unsigned int pos = std::distance(tmp_groups.begin(),it);
          for(unsigned int k : neighbors[i][j]) {
            auto it2 = std::find(tmp_neighbors[pos].begin(), tmp_neighbors[pos].end(), k);
            if(it2 == tmp_neighbors[pos].end()) tmp_neighbors[pos].push_back(k);
          }
        }
      }
    }
    boundary_groups.push_back(tmp_groups);
    boundary_neighbors.push_back(tmp_neighbors);
  }

  // DEBUG
  if(DEBUG) {
    cout << "[DEBUG]: boundary groups:" << endl;
    for(unsigned int i=0; i<num_haps; i++) {
      for(unsigned int g=0; g<boundary_groups[i].size(); g++) {
        cout << " " << boundary_groups[i][g];
        cout << " ( ";
        for(unsigned int i2 : boundary_neighbors[i][g]) {
          cout << i2 << " ";
        }
        cout << ") ";
      }
      cout << endl;
    }
  }

  // Make trivial knockoffs in the boundary groups
  if(DEBUG) cout << "[DEBUG]: making trivial knockoffs in IBD boundary groups:" << endl;
  for(unsigned int s=0; s<segments_i.size(); s++) {
    for(unsigned int i : segments_i[s]) {
      unsigned int j_min = segments_j_min[s];
      unsigned int g_min = groups[j_min];
      for(unsigned int j : elements[g_min]) {
        Zk[i][j] = Z[i][j];
      }
      unsigned int j_max = segments_j_max[s];
      unsigned int g_max = groups[j_max];
      for(unsigned int j : elements[g_max]) {
        Zk[i][j] = Z[i][j];
      }
    }
  }

  // Generate knockoffs in the IBD groups
  if(DEBUG) cout << "[DEBUG]: generating knockoffs in IBD groups:" << endl;
  for(unsigned int s=0; s<segments_i.size(); s++) {
    unsigned int i0 = segments_i[s][0];
    unsigned int j_min = segments_j_min[s];
    unsigned int g_begin = groups[j_min]+1;
    if(groups[j_min]==0) g_begin = 0;
    unsigned int j_max = segments_j_max[s];
    unsigned int g_end = groups[j_max]-1;
    if(groups[j_max]==(num_groups-1)) g_end = groups[j_max];
    if(g_end >= g_begin) {
      // Generate knockoff for the representative individual
      sample_knockoff_MC(i0, Z[i0], Zk[i0], g_begin, g_end, elements, groups);
      // Copy the knockoffs to the other IBD-sharing individuals
      for(unsigned int l=1; l<segments_i[s].size(); l++) {
        unsigned int i = segments_i[s][l];
        for(unsigned int g=g_begin; g<=g_end; g++) {
          for(unsigned int j : elements[g]) {
            Zk[i][j] = Zk[i0][j];
          }
        }
     }
    }
  }

  // Store group endpoints for the IBD segments of each individual
  vector< vector<unsigned int> > inner_g_left(num_haps);
  vector< vector<unsigned int> > inner_g_right(num_haps);
  for(unsigned int s=0; s<segments_i.size(); s++) {
    unsigned int j_min = segments_j_min[s];
    unsigned int g_begin = groups[j_min];
    unsigned int j_max = segments_j_max[s];
    unsigned int g_end = groups[j_max];
    for(unsigned int i : segments_i[s]) {
      inner_g_left[i].push_back(g_begin);
      inner_g_right[i].push_back(g_end);
    }
  }
  for(unsigned int i=0; i<num_haps; i++) {
    std::sort(inner_g_left[i].begin(), inner_g_left[i].end());
    std::sort(inner_g_right[i].begin(), inner_g_right[i].end());
  }

  if(DEBUG) {
    cout << "[DEBUG]: IBD groups:" << endl;
    for(unsigned int i=0; i<num_haps; i++) {
      cout << i << " : ";
      for(unsigned int g=0; g<inner_g_left[i].size(); g++) {
        cout << " " << inner_g_left[i][g] << "--" << inner_g_right[i][g];
      }
      cout << endl;
    }
  }

  // Get group endpoints for the complement of the IBD segments of each individual
  vector< vector<unsigned int> > outer_g_left(num_haps);
  vector< vector<unsigned int> > outer_g_right(num_haps);
  for(unsigned int i=0; i<num_haps; i++) {
    unsigned int g_begin = 0;
    unsigned int g_end = num_groups-1;
    for(unsigned int g=0; g<inner_g_left[i].size(); g++) {
      g_end = inner_g_left[i][g]-1;
      if(g_end < MAXINT) {
        if(g_begin<=g_end) {
          outer_g_left[i].push_back(g_begin);
          outer_g_right[i].push_back(g_end);
        }
      }
      g_begin = inner_g_right[i][g]+1;
    }
    if(g_end<num_groups) {
      g_end = num_groups-1;
      if(g_begin<num_groups) {
        if( (g_begin<=g_end) && (g_end < MAXINT)) {
          outer_g_left[i].push_back(g_begin);
          outer_g_right[i].push_back(g_end);
        }
      }
    }
  }

  if(DEBUG) {
    cout << "[DEBUG]: complement of IBD groups:" << endl;
    for(unsigned int i=0; i<num_haps; i++) {
      cout << i << " : ";
      for(unsigned int s=0; s<outer_g_left[i].size(); s++) {
        cout << " " << outer_g_left[i][s] << "--" << outer_g_right[i][s];
        if(outer_g_right[i][s] > 1000) {
          cout << "DEBUG: error!" << endl;
          exit(1);
        }
      }
      cout << endl;
    }
  }

  // Generate knockoffs for groups outside the IBD segments
  if(DEBUG) cout << "[DEBUG]: generating knockoffs outside IBD groups:" << endl;
  for(unsigned int i=0; i<num_haps; i++) {
    for(unsigned int s=0; s<outer_g_left[i].size(); s++) {
      unsigned int g_begin = outer_g_left[i][s];
      unsigned int g_end = outer_g_right[i][s];
      sample_knockoff_MC(i, Z[i], Zk[i], g_begin, g_end, elements, groups);
    }
  }

  // Make sure we didn't miss anything
  for(unsigned int i=0; i<num_haps; i++) {
    for(unsigned int j=0; j<num_snps; j++) {
      if(Zk[i][j]==-1) {
        throw_error_bug("Missing knockoff variants");
      }
    }
  }

  // Make sure knockoffs agree on IBD segments
  for(unsigned int i=0; i<num_haps; i++) {
    for(unsigned int j=0; j<num_snps; j++) {
      for(unsigned int i2 : neighbors[i][j]) {
        if(Zk[i][j]!=Zk[i2][j]) {
          cout << "Warning: inconsistency found on Zk (" << i << ", " << j << ", " << i2 << ")" << endl;
          throw_error_bug("Inconsistency found on Zk");
        }
      }
    }
  }

}

void FamilyKnockoffs::sample_knockoff_MC(unsigned int i, const vector<unsigned int> & z, vector<unsigned int> & zk, unsigned int g_begin, unsigned int g_end,
                                         const vector< vector<unsigned int> > & elements, const vector<unsigned int> & groups)
  const {

  unsigned int num_groups = elements.size();

  // Intermediate variables for knockoff generation
  vector<double> weights_mc(K);
  vector<double> N(K);
  vector<double> N_old(K);
  std::vector<double> vstar;
  vector< vector<double> > vbar;

  // Initalize normalization functions
  std::fill(N.begin(), N.end(), 1.0);
  std::fill(N_old.begin(), N_old.end(), 1.0);
  double N_min = 1.0e-10;

  for(unsigned int g=g_begin; g<=g_end; g++) {
    // Compute vstar
    unsigned int groupSize = elements[g].size();
    vstar.resize(groupSize,0.0);
    for(unsigned int j=groupSize; j-- > 0; ) {
      if(j < groupSize-1) {
        vstar[j] = vstar[j+1] * b[elements[g][j+1]];
      }
      else {
        if(g < num_groups-1) {
          vstar[j] = b[elements[g][0]];
        }
        else {
          vstar[j] = 0.0;
        }
      }
    }

    // Compute vbar matrix
    vbar = vector< vector<double> >(groupSize, std::vector<double> (K,0));
    for(unsigned int j=groupSize; j-- > 0; ) {

      double sum_a = 0;
      if(j < groupSize-1) {
        for(unsigned int l=0; l<K; l++) {
          sum_a += (1.0-b[elements[g][j+1]]) * alpha[i][l];
        }
      }
      for(unsigned int z=0; z<K; z++) {
        if(j < groupSize-1) {
          vbar[j][z] = vbar[j+1][z] * (sum_a + b[elements[g][j+1]]);
          vbar[j][z] += vstar[j+1] * (1.0-b[elements[g][j+1]]) * alpha[i][z];
        }
        else {
          if(g < num_groups-1) {
            vbar[j][z] = (1.0-b[elements[g+1][0]]) * alpha[i][z];
          }
          else {
            vbar[j][z] = 1.0;
          }
        }
      }
    }

    // Precompute sum for partition function
    double N_sum = 0;
    for(unsigned int k=0; k<K; k++) {
      if(g==0) {
        N_sum += (1.0-b[elements[g][0]]) * alpha[i][k];
      }
      else if (g==g_begin) {
        unsigned int z0 = z[elements[g-1].back()];
        double tmp = ((1.0-b[elements[g][0]])* alpha[i][k] + b[elements[g][0]]*(double)(k==z0));
        N_sum += tmp / N_old[k];
      }
      else {
        unsigned int z0 = z[elements[g-1].back()];
        unsigned int z0k = zk[elements[g-1].back()];
        double tmp = ((1.0-b[elements[g][0]])* alpha[i][k] + b[elements[g][0]]*(double)(k==z0));
        tmp  = tmp * ((1.0-b[elements[g][0]])* alpha[i][k] + b[elements[g][0]]*(double)(k==z0k));
        N_sum += tmp / N_old[k];
      }
    }

    // Compute partition function
    for(unsigned int k=0; k<K; k++) {
      if(g < num_groups-1) {
        N[k] += vbar[0][k] * N_sum;
        if(g==0) {
          N[k] += vstar[0] * (1.0-b[elements[g][0]]) * alpha[i][k];
        }
        else if (g==g_begin) {
          unsigned int z0 = z[elements[g-1].back()];
          double tmp = ((1.0-b[elements[g][0]]) * alpha[i][k] + b[elements[g][0]]*(double)(k==z0));
          N[k] += vstar[0] * tmp / N_old[k];
        }
        else {
          unsigned int z0 = z[elements[g-1].back()];
          unsigned int z0k = zk[elements[g-1].back()];
          double tmp = ((1.0-b[elements[g][0]]) * alpha[i][k] + b[elements[g][0]]*(double)(k==z0));
          tmp  = tmp * ((1.0-b[elements[g][0]]) * alpha[i][k] + b[elements[g][0]]*(double)(k==z0k));
          N[k] += vstar[0] * tmp / N_old[k];
        }
      }
    }

    // Normalize partition function and make sure we avoid division by zero
    double N_norm = 0;
    for(unsigned int k=0; k<K; k++) {
      if(N[k] < N_min) {
        N[k] = N_min;
      }
      N_norm += N[k];
    }
    for(unsigned int k=0; k<K; k++) {
      N[k] /= N_norm;
      if(N[k] < N_min) {
        N[k] = N_min;
      }
    }

    // Compute sampling weights
    for(unsigned int j=0; j<groupSize; j++) {
      std::fill(weights_mc.begin(), weights_mc.end(), 1.0);
      for(unsigned int k=0; k<K; k++) {
        if(j>0) {
          unsigned int z0k = zk[elements[g][j-1]];
          weights_mc[k] *= ((1.0-b[elements[g][j]]) * alpha[i][k] + b[elements[g][j]]*(double)(k==z0k));
        }
        else {
          if( g==0 ) {
            weights_mc[k] *= (1.0-b[elements[g][0]]) * alpha[i][k];
          }
          else if (g==g_begin) {
            unsigned int z0 = z[elements[g-1].back()];
            weights_mc[k] *= ((1.0-b[elements[g][0]]) * alpha[i][k] + b[elements[g][0]]*(double)(k==z0));
          }
          else {
            unsigned int z0 = z[elements[g-1].back()];
            unsigned int z0k = zk[elements[g-1].back()];
            weights_mc[k] *= ((1.0-b[elements[g][0]]) * alpha[i][k] + b[elements[g][0]]*(double)(k==z0));
            weights_mc[k] *= ((1.0-b[elements[g][0]]) * alpha[i][k] + b[elements[g][0]]*(double)(k==z0k));
          }
          weights_mc[k] /= N_old[k];
        }
        if(g == num_groups-1) {
          weights_mc[k] *= (vbar[j][0]);
        }
        else {
          unsigned int z1 = z[elements[g+1][0]];
          weights_mc[k] *= (vbar[j][z1] + vstar[j] * (double)(k==z1));
        }
      }
      zk[elements[g][j]] = weighted_choice(weights_mc, rng);
    }

    for(unsigned int k=0; k<K; k++) {
      N_old[k] = N[k];
    }
  }

}

void FamilyKnockoffs::print_Z(vector< vector<unsigned int> > & Zk) const {
  // DEBUG: write Z to text file
  ofstream myfile;
  myfile.open("tmp/Zk.txt");
  for(unsigned int i=0; i<num_haps; i++) {
    for(unsigned int j=0; j<num_snps; j++) {
      if(j>0) myfile << " ";
      myfile << Zk[i][j];
    }
    myfile << endl;
  }
  myfile.close();
}

#endif
