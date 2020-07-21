#ifndef _KINSHIP_CPP
#define _KINSHIP_CPP

#include "kinship.h"

Kinship::Kinship() {
}

Kinship::~Kinship() {
}

Kinship::Kinship(const vector<string>& ref_files_local, const vector<string>& ref_files_global) {
  // FIXME: use local references here
  assert(ref_files_local.size() == ref_files_global.size());
  num_chrs = ref_files_local.size();
  assert(num_chrs>0);

  cout << "Loading local kinship information from:" << endl;
  for(int chr=0; chr<num_chrs; chr++) {
    cout << "  " << ref_files_local[chr] << endl;
  }

  references_local.resize(num_chrs);
  load_references_local(ref_files_local);

  cout << "Loading global kinship information from:" << endl;
  for(int chr=0; chr<num_chrs; chr++) {
    cout << "  " << ref_files_global[chr] << endl;
  }

  references_global.resize(num_chrs);
  load_references_global(ref_files_global);
}

void Kinship::load_references_local(const vector<string>& ref_files) {
  references_local.resize(num_chrs);

  for(int chr=0; chr<num_chrs; chr++) {
    // Select chromosome ref file
    string buffer;
    vector <string> tokens;
    ifile fd_ref(ref_files[chr]);

    // Read header (number of haplotypes)
    getline(fd_ref, buffer, '\n');
    sutils::tokenize(buffer, tokens);
    if(tokens.size() == 2) {
      const int _num_haps = std::stoi(tokens[1]);
      if(tokens[0] != "#") throw_error_input("reading reference file, incorrect header");
      if(_num_haps > 0) {
        num_haps = _num_haps;
      } else {
        throw_error_input("reading reference file, incorrect number of haplotypes in header");
      }
    } else {
      throw_error_input("reading reference file, incorrect number of fields on line 1");
    }

    // Read header (number of windows)
    getline(fd_ref, buffer, '\n');
    sutils::tokenize(buffer, tokens);
    if(tokens.size() == 2) {
      if(tokens[0] != "#") throw_error_input("reading reference file, incorrect header");
      const int _num_windows = std::stoi(tokens[1]);
      if(_num_windows > 0) {
        num_windows = _num_windows;
      } else {
        throw_error_input("reading reference file, incorrect number of haplotypes in header");
      }
    } else {
      throw_error_input("reading reference file, incorrect number of fields on line 2");
    }

    // Read value of K
    getline(fd_ref, buffer, '\n');
    getline(fd_ref, buffer, '\n');
    sutils::tokenize(buffer, tokens);
    K = tokens.size();
    if(K<0) {
      throw_error_input("reading reference file, incorrect header");
    }
    fd_ref.close();

    cout << "Header information: " << endl;
    cout << "  number of haplotypes      : " << num_haps << endl;
    cout << "  number of genomic windows : " << num_windows << endl;
    cout << "  number of HMM references  : " << K << endl;

    // Allocate space for local references
    references_local[chr].resize(num_haps, ivector2d(num_windows, vector<int>(K)));

    // Open file again
    fd_ref.open(ref_files[chr]);

    // Read body, skipping first two lines
    int i = 0;      // Haplotype index
    int w = 0;      // Window index
    int line = -2;  // File row, starting from file body
    while (getline(fd_ref, buffer, '\n')) {
      if(line>=0) {
        sutils::tokenize(buffer, tokens);
        if(line % 2 == 0) {
          if(tokens[0] != "#") throw_error_input("reading reference file on line " + std::to_string(line));
        } else {
          if(w>=num_windows) {
            w=0;
            i++;
          }
          if(i>=num_haps) throw_error_input("reading reference file, too long.");

          K = tokens.size();
          if(K != tokens.size()) {
            string err_msg = "reading reference file on line " + std::to_string(line);
            err_msg += ", incorrect line length";
            throw_error_input(err_msg);
          }
          assert(K==tokens.size());
          for(int k=0; k<K; k++) {
            references_local[chr][i][w][k] = std::stoi(tokens[k]);
          }
          w++;
        }
      }
      line++;
    }

    if((i+1)<num_haps) {
      cout << "i = " << i << endl;
      throw_error_input("reading reference file, not enough rows.");
    }

    cout << "Kinship information loaded." << endl << endl;
    fd_ref.close();
  }
}

void Kinship::load_references_global(const vector<string>& ref_files) {
  references_global.resize(num_chrs);

  for(int chr=0; chr<num_chrs; chr++) {
    // Select chromosome ref file
    string buffer;
    vector <string> tokens;
    ifile fd_ref(ref_files[chr]);

    // Read header (number of haplotypes)
    getline(fd_ref, buffer, '\n');
    sutils::tokenize(buffer, tokens);
    if(tokens.size() == 2) {
      const int _num_haps = std::stoi(tokens[1]);
      if(tokens[0] != "#") throw_error_input("reading reference file, incorrect header");
      if(_num_haps > 0) {
        num_haps = _num_haps;
      } else {
        throw_error_input("reading reference file, incorrect number of haplotypes in header");
      }
    } else {
      throw_error_input("reading reference file, incorrect number of fields on line 1");
    }

    // Read header (number of windows)
    getline(fd_ref, buffer, '\n');
    sutils::tokenize(buffer, tokens);
    if(tokens.size() == 2) {
      if(tokens[0] != "#") throw_error_input("reading reference file, incorrect header");
      const int _num_windows = std::stoi(tokens[1]);
      if(_num_windows > 0) {
        num_windows = _num_windows;
      } else {
        throw_error_input("reading reference file, incorrect number of haplotypes in header");
      }
    } else {
      throw_error_input("reading reference file, incorrect number of fields on line 2");
    }

    // Read value of K
    getline(fd_ref, buffer, '\n');
    getline(fd_ref, buffer, '\n');
    sutils::tokenize(buffer, tokens);
    K = tokens.size();
    if(K<0) {
      throw_error_input("reading reference file, incorrect header");
    }
    fd_ref.close();

    cout << "Header information: " << endl;
    cout << "  number of haplotypes      : " << num_haps << endl;
    cout << "  number of HMM references  : " << K << endl;

    // Allocate space for global references
    references_global[chr].resize(num_haps, vector<int>(K));

    // Open file again
    fd_ref.open(ref_files[chr]);

    // Read body, skipping first two lines
    int i = 0;      // Haplotype index
    int line = -2;  // File row, starting from file body
    while (getline(fd_ref, buffer, '\n')) {
      if(line>=0) {
        sutils::tokenize(buffer, tokens);
        if(line % 2 == 0) {
          if(tokens[0] != "#") throw_error_input("reading reference file on line " + std::to_string(line));
        } else {
          if(i>=num_haps) throw_error_input("reading reference file, too long.");

          K = tokens.size();
          if(K != tokens.size()) {
            string err_msg = "reading reference file on line " + std::to_string(line);
            err_msg += ", incorrect line length";
            throw_error_input(err_msg);
          }
          assert(K==tokens.size());
          for(int k=0; k<K; k++) {
            references_global[chr][i][k] = std::stoi(tokens[k]);
          }
          i++;
        }
      }
      line++;
    }

    if((i+1)<num_haps) {
      cout << "i = " << i << endl;
      throw_error_input("reading reference file, not enough rows.");
    }

    cout << "Kinship information loaded." << endl << endl;
    fd_ref.close();
  }
}

Kinship::Kinship(const vector<Metadata>& metadata, const Covariates & covariates, vector<string> output_files,
                 int compression, int cluster_size_min, int cluster_size_max, int _K, int num_threads) {
  K = _K;

  //cout << "[DEBUG] in Kinship(with covariates)" << endl;

  // Initialize kinship computer and run K-means clustering
  KinshipComputer computer(metadata, compression, covariates, cluster_size_min, cluster_size_max, num_threads);

  // // Debug: print clusters
  computer.writeClusters(output_files);

  num_chrs = metadata.size();

  // Assign local and global references
  references_local.resize(num_chrs);
  references_global.resize(num_chrs);
  for(int chr=0; chr<num_chrs; chr++) {
    const int num_snps = metadata[chr].num_snps();

    const Windows& windows_local = metadata[chr].windows;
    references_local[chr] = computer.assign_references(K, chr, windows_local);

    if(windows_local.num_windows>1) {
      const Windows windows_global(num_snps);
      ivector3d ref_tmp = computer.assign_references(K, chr, windows_global);
      references_global[chr].resize(ref_tmp.size());
      for(int i=0; i<ref_tmp.size(); i++) {
        references_global[chr][i] = ref_tmp[i][0];
      }
    } else {
      references_global[chr].resize(references_local[chr].size());
      for(int i=0; i<references_local[chr].size(); i++) {
        references_global[chr][i] = references_local[chr][i][0];
      }
    }
  }

  // Count numbers of haplotypes and chromosomes
  assert(references_local.size() == num_chrs);
  num_haps = references_local[0].size();
  assert(num_haps > 0);
  num_windows = references_local[0][0].size();
  sanitycheck_references();
}

Kinship::Kinship(const vector<Metadata>& metadata, vector<string> output_files,
                 int compression, int cluster_size_min, int cluster_size_max, int _K, int num_threads) {
  K = _K;

  //cout << "[DEBUG] in Kinship(without covariates)" << endl;

  // Initialize kinship computer and run K-means clustering
  KinshipComputer computer(metadata, compression, cluster_size_min, cluster_size_max, num_threads);

  // Debug: print clusters
  computer.writeClusters(output_files);

  num_chrs = metadata.size();

  // Assign local and global references
  references_local.resize(num_chrs);
  references_global.resize(num_chrs);
  for(int chr=0; chr<num_chrs; chr++) {
    const int num_snps = metadata[chr].num_snps();

    const Windows& windows_local = metadata[chr].windows;
    references_local[chr] = computer.assign_references(K, chr, windows_local);

    if(windows_local.num_windows>1) {
      const Windows windows_global(num_snps);
      ivector3d ref_tmp = computer.assign_references(K, chr, windows_global);
      references_global[chr].resize(ref_tmp.size());
      for(int i=0; i<ref_tmp.size(); i++) {
        references_global[chr][i] = ref_tmp[i][0];
      }
    } else {
      references_global[chr].resize(references_local[chr].size());
      for(int i=0; i<references_local[chr].size(); i++) {
        references_global[chr][i] = references_local[chr][i][0];
      }
    }
  }

  // Count numbers of haplotypes and chromosomes
  assert(references_local.size() == num_chrs);
  assert(num_chrs > 0);
  num_haps = references_local[0].size();
  assert(num_haps > 0);
  num_windows = references_local[0][0].size();
  sanitycheck_references();
}

void Kinship::sanitycheck_references() const {
  assert(num_windows > 0);
  for(int chr=0; chr<num_chrs; chr++) {
    assert(references_local[chr].size()==num_haps);
    for(int i=0; i<num_haps; i++) {
      for(int w=0; w<num_windows; w++) {
        assert(references_local[chr][i][w].size()==K);
      }
    }
  }
}

const ivector3d & Kinship::get_references(int chr) const {
  return(references_local[chr]);
}

const ivector2d & Kinship::get_references_global(int chr) const {
  return(references_global[chr]);
}

int Kinship::get_K() const {
  return(K);
}

void Kinship::writeReferences(const vector<Metadata>& metadata, const vector<string> & out_file_names) const {
  writeReferencesLocal(metadata, out_file_names);
  writeReferencesGlobal(metadata, out_file_names);
}

void Kinship::writeReferencesLocal(const vector<Metadata>& metadata, const vector<string> & out_file_names) const {
  cout << "Individual global references written to:" << endl;
  for(int chr=0; chr<num_chrs; chr++) {
    string out_file_name = out_file_names[chr];
    out_file_name.append("_lref.txt");
    ofstream outfile(out_file_name.c_str());
    if (!outfile.is_open()){
      cout << "Problem creating the output file: " << out_file_name;
      cout <<"Either the directory does not exist or you do not have write permissions." << endl;
    }
    const int num_windows = metadata[chr].windows.num_windows;
    outfile << "# " << num_haps << endl;
    outfile << "# " << num_windows << endl;
    for(int i=0; i<num_haps; i++) {
      for(int w=0; w<num_windows; w++) {
        outfile << "# ID " << metadata[chr].sample_filter.ID[i/2] << " " << (i % 2) << " chr " << chr;
        outfile << " window " << w << " ";
        outfile << metadata[chr].windows.start[w] << " " << metadata[chr].windows.end[w] << endl;
        for(int k=0; k<K; k++) {
          outfile << references_local[chr][i][w][k];
          if(k+1<K) outfile <<" ";
        }
        outfile << endl;
      }
    }
    outfile.close();
    cout << "  " << out_file_name << endl;
  }
  cout << endl;
}

void Kinship::writeReferencesGlobal(const vector<Metadata>& metadata,
                                    const vector<string> & out_file_names) const {
  cout << "Individual local references written to:" << endl;
  for(int chr=0; chr<num_chrs; chr++) {
    string out_file_name = out_file_names[chr];
    out_file_name.append("_ref.txt");
    ofstream outfile(out_file_name.c_str());
    if (!outfile.is_open()){
      cout << "Problem creating the output file: " << out_file_name;
      cout <<"Either the directory does not exist or you do not have write permissions." << endl;
    }
    const int num_windows = 1;
    outfile << "# " << num_haps << endl;
    outfile << "# " << num_windows << endl;
    for(int i=0; i<num_haps; i++) {
      outfile << "# ID " << metadata[chr].sample_filter.ID[i/2] << " " << (i % 2) << " chr " << chr;
      outfile << " window " << 0 << " ";
      outfile << 0 << " " << metadata[chr].legend_filter.size() << endl;
      for(int k=0; k<K; k++) {
        outfile << references_global[chr][i][k];
        if(k+1<K) outfile <<" ";
      }
      outfile << endl;
    }
    outfile.close();
    cout << "  " << out_file_name << endl;
  }
  cout << endl;
}

#endif
