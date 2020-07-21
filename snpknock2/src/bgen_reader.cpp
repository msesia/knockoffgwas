#ifndef _BGEN_READER_CPP
#define _BGEN_READER_CPP

#include "bgen_reader.h"

BgenReader::BgenReader(const string& filename_, const string& sample_filename) {
  filename = filename_;

  // Read sample IDs from ".sample" file
  BgenParser bgenParser;
  try {
    bgenParser = BgenParser(filename);
  } catch (int e) {
    std::cout << "Exception occurred while trying to open BGEN file " << filename << std::endl;
    exit(1);
  }

  string buffer;
  vector<string> tokens;
  int line = -2;
  ifile file(sample_filename);
  while (getline(file, buffer, '\n')) {
    if(line>=0) {
      sutils::tokenize(buffer, tokens);
      sample.push_back(tokens);
    }
    line++;
  }
  file.close();

  // Make sure that the number of samples matches the dimensions of the BGEN file
  if(sample.size() != bgenParser.number_of_samples()) {
    cerr << "Error loading BGEN file: the number of samples does not match." << endl;
    exit(1);
  }
}

void BgenReader::summarise() {
  BgenParser bgenParser;
  try {
    bgenParser = BgenParser(filename);
  } catch (int e) {
    std::cout << "Exception occurred while trying to open BGEN file " << filename << std::endl;
    exit(1);
  }
  bgenParser.summarise(std::cerr);
}

std::pair<bool, bool> BgenReader::prob_to_hap(const vector<double>& probs) const {
  int K = probs.size();
  int h1, h2;
  assert(K==4);

  if((probs[0]==1) && (probs[1]==0)) h1 = 0;
  else if((probs[0]==0) && (probs[1]==1)) h1 = 1;
  else throw("Incorrect input");

  if((probs[2]==1) && (probs[3]==0)) h2 = 0;
  else if((probs[2]==0) && (probs[3]==1)) h2 = 1;
  else throw("Incorrect input");

  return std::pair<bool, bool>(h1,h2);
}

void BgenReader::read(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H,
                      bool verbose) {
  read(sample_ids, snp_ids, H, verbose, 1);
}

void BgenReader::read(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H,
                      bool verbose, int nthreads) {

  // Initialize data container
  int num_samples = sample_ids.size();
  int num_snps = snp_ids.size();
  int num_haps = 2 * num_samples;
  H.clear();
  H = vector <chaplotype>(num_haps, chaplotype(num_snps));

  read_mt(sample_ids, snp_ids, H, verbose, nthreads);
}

void BgenReader::read_mt(const vector<string>& sample_ids, const vector<string>& snp_ids, vector<chaplotype>& H,
                         bool verbose, int nthreads) {
  int nsnps = snp_ids.size();

  // Sanity check
  if(nthreads > nsnps) nthreads = snp_ids.size();

  // Divide list of SNPs
  int m = nsnps / nthreads;
  vector< vector<string> > sub_snp_ids;
  for(int i=0; i<nthreads; i++) {
    int left, right;
    left = i * m;
    if(i<(nthreads-1)) {
      right = (i+1) * m;
    } else {
      right = nsnps;
    }
    std::vector<string> tmp_ids (snp_ids.begin()+left, snp_ids.begin()+right);
    sub_snp_ids.push_back(tmp_ids);
  }

  // cout << "Assignments: " << endl;
  // for(int i=0; i<nthreads; i++) {
  //   cout << "Thread " << i << " : ";
  //   for(int j=0; j<sub_snp_ids[i].size(); j++) {
  //     cout << sub_snp_ids[i][j] << " ";
  //   }
  //   cout << endl;
  // }

  vector<boost::thread> workers;
  vector<int> progress(nthreads+1,0);

  // Create workers  
  cout << "Reading BGEN file using " << nthreads;
  if(nthreads>1) cout << " threads:" << endl;
  else cout << " thread:" << endl;

  if(nthreads>1) {
    for(int w=0; w<nthreads; w++) {
      //read_worker(sample_ids, snp_ids, sub_snp_ids[i], H, verbose, w, progress);
      workers.push_back(boost::thread(&BgenReader::read_worker, this, 
                                      boost::ref(sample_ids), boost::ref(snp_ids), boost::ref(sub_snp_ids[w]), 
                                      boost::ref(H), verbose, w, boost::ref(progress)));
    }
    // Launch workers
    for(int w=0; w<nthreads; w++ ) {
      workers[w].join();
    }
  } else {
    read_worker(sample_ids, snp_ids, snp_ids, H, verbose, 0, progress);
  }

}

void BgenReader::read_worker(const vector<string>& sample_ids, const vector<string>& abs_snp_ids,
                             const vector<string>& snp_ids,
                             vector<chaplotype>& H, bool verbose, int wid, vector<int> &progress) {
  BgenParser bgenParser;
  try {
    bgenParser = BgenParser(filename);
  } catch (int e) {
    std::cout << "Exception occurred while trying to open BGEN file " << filename << std::endl;
    exit(1);
  }

  // Dimensions of BGEN data
  int num_samples = sample_ids.size();
  int num_snps = snp_ids.size();
  int num_haps = 2 * num_samples;
  int num_samples_bgen = (int)bgenParser.number_of_samples();
  int num_snps_bgen = (int)bgenParser.number_of_variants();

  // Sanity check
  assert(num_samples <= num_samples_bgen);

  // Keep track of loaded variants
  vector<bool> loaded_snp(num_snps, false);

  // Output variables for BgenParser
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector< std::string > alleles;
  std::vector< std::vector< double > > probs;

  // Find BGEN row indices of requested samples
  vector<int> sample_idx(num_samples, -1);
  match_indices(sample.ID, sample_ids, sample_idx);

  // Initialize progress bar (worker 0)
  int progress_nsteps = 100;
  int progress_max = num_snps;
  int progress_period = std::max(1, (int)(progress_max/progress_nsteps));
  if(wid==0 && verbose) {
    cout << "|" << flush;
    for(int s=0; s<progress_nsteps; s++) cout << ".";
    cout << "|" << endl;
    cout << "|" << flush;
  }

  int num_loaded = 0;
  for(int j_bgen=0; j_bgen<num_snps_bgen; j_bgen++) {
    // Read variant description from file
    bgenParser.read_variant(&chromosome, &position, &rsid, &alleles);

    // Load this column only if the variant was requested
    auto it_rel = std::find(snp_ids.begin(), snp_ids.end(), rsid);
    if(it_rel == snp_ids.end()) {
      bgenParser.ignore_probs();
    } else {
      // Read variant data from file
      bgenParser.read_probs(&probs);

      // Sanity check
      assert((int)probs.size() == num_samples_bgen);

      // Find index of variant in storage matrix
      auto it_abs = std::find(abs_snp_ids.begin(), abs_snp_ids.end(), rsid);
      int j_abs = std::distance(abs_snp_ids.begin(), it_abs);

      // Insert data into container
      for(int i=0; i<num_samples; i++) {
        int i_bgen = sample_idx[i];
        assert(i_bgen>=0);
        std::pair<bool, bool> h_pair = prob_to_hap(probs[i_bgen]);
        H[2*i].set(j_abs, h_pair.first);
        H[2*i+1].set(j_abs, h_pair.second);
      }

      // Find relative index of variant in list of requested SNPs
      int j_rel = std::distance(snp_ids.begin(), it_rel);

      // Mark the variant as loaded
      loaded_snp[j_rel] = true;
      progress[wid]++;
      num_loaded++;

      // Update progress bar (worker 0)
      if(wid==0 && verbose) {
        if((num_loaded % progress_period == 0)) {          
          int cum_progress = std::accumulate(progress.begin(),progress.end()-1,0);
          int expected_bars = (float)(cum_progress) / (float)(abs_snp_ids.size()) * progress_nsteps;          
          //cout << "Expected bars: " << expected_bars << ", drawn bars " << progress.back() << endl;
          int missing_bars = expected_bars - progress.back();
          for(int i=0; i<missing_bars; i++) {
             progress.back()++;
             cout << "=" << flush;
          }
        }
      }
    }

    // Check whether we are done
    if(num_loaded==num_snps) break;

  }

  // Make sure all requested variants have been loaded
  auto it = std::find(loaded_snp.begin(), loaded_snp.end(), false);
  if(it != loaded_snp.end()) {
    cerr << "Error loading BGEN file. Some requested variants were not found.";
    exit(1);
  }

  // Finalize progress bar (last worker to finish)
  if(verbose) {
    int cum_progress = std::accumulate(progress.begin(),progress.end()-1,0);
    if(cum_progress==abs_snp_ids.size()) {
      int missing_bars = progress_nsteps - progress.back();    
      //cout << "Missing bars: " << missing_bars << endl;
      for(int i=0; i<missing_bars; i++) {
        cout << "=" << flush;
      }
      cout << "|" << endl;
    }
  }
}

void BgenReader::print(const vector<chaplotype>& H) const {
  // Print transposed haplotypes (HAPS format)
  int num_haps = H.size();
  assert(num_haps>0);
  int num_snps = H[0].size();
  for(std::size_t j = 0; j < num_snps; j++) {
    for(std::size_t i = 0; i < num_haps; i++) {
      if(i>0) std::cout << " ";
      std::cout << H[i][j];
    }
    std::cout << "\n";
  }
}

#endif
