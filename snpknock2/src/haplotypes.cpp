#ifndef _HAPLOTYPES_CPP
#define _HAPLOTYPES_CPP

#include "haplotypes.h"

Haplotypes::Haplotypes(const Metadata& _metadata) : metadata(_metadata) {
  num_haps = metadata.num_haps();
  num_snps = metadata.num_snps();
}

Haplotypes::~Haplotypes() {
}

void Haplotypes::load_data(bool verbose, int nthreads) {
  if(metadata.df.format=="haps") {
    HapsReader hapsReader(metadata.df.data + ".haps", metadata.df.sample, metadata.df.legend);
    hapsReader.read(metadata.sample_filter.ID, metadata.legend_filter.ID, H, verbose);
  } else if (metadata.df.format=="bgen") {
    BgenReader bgenReader(metadata.df.data + ".bgen", metadata.df.sample);
    //bgenReader.read(metadata.sample_filter.ID, metadata.legend_filter.ID, H, verbose);
    bgenReader.read(metadata.sample_filter.ID, metadata.legend_filter.ID, H, verbose, nthreads);
  } else {
    std::cerr << "Error: unknown input data format " << metadata.df.format << std::endl;
    exit(1);
  }
  // Verify that the loaded data has the correct dimensions
  assert(H.size()==num_haps);
  assert(H.size()>0);
  assert(H[0].size()==num_snps);
}

void Haplotypes::push_back(const vector<int> & input) {
  H.push_back(chaplotype(input));
}

void Haplotypes::push_back(const vector<bool> & input) {
  H.push_back(chaplotype(input));
}

void Haplotypes::writeHaplotypes(string _filename) const {
  string filename = _filename;
  filename.append(".haps");
  ofstream outfile(filename.c_str());
  if (!outfile.is_open()){
    cout << "Problem creating the output file. ";
    cout <<"Either the directory does not exist or you do not have write permissions." << endl;
  }
  for(int j=0; j<num_snps; j++) {
    for(int i=0; i<num_haps; i++) {
      outfile << H[i].get(j);
      if(i+1<num_haps) outfile <<" ";
    }
    outfile << std::endl;
  }
  outfile.close();
}

int Haplotypes::rpartition(double w, vector<int> & output){

  double a = 0;
  double b = num_snps;
  double chrlen = b - a;
  w = chrlen / floor( chrlen/w );
  int nwin = chrlen/w;
  if(nwin==0) nwin=1;
  //  cout << nwin << "WINDOWS" << endl;
  int nsnp = num_snps/nwin;
  output.resize(nwin+1);
  output[0] = 0;
  output[nwin] = num_snps;
  for(int i=1;i<nwin;i++)
    output[i] = i*nsnp + 0.1 * (float)putils::getRandom(nsnp) -.05*nsnp; // some randomness in choice of windows. Do we want this?

  int prev = output[0];
  for(int i=1;i<(1+nwin);i++) {
    //  cout << i <<".\t"<< vec_pos[output[i]-1]->bp/1000000. <<"mb\t" << prev<< " - " << output[i] << "\t"  << output[i] - prev << " snps" << endl;
    assert( (output[i] - prev) > 0 );
    prev = output[i];
  }
  //cout << "rpartition finished."  << endl;
  return(nwin);

}

int Haplotypes::get_num_snps() const{
  return(num_snps);
}

int Haplotypes::get_num_haps() const{
  return(num_haps);
}

int Haplotypes::distance(const int i1, const int i2) const {
  return(H[i1].hamming(H[i2]));
}

int Haplotypes::distance(const int i1, const int i2, const int j_start, const int j_end) const {
  return(H[i1].hamming(H[i2], j_start, j_end));
}

void Haplotypes::compute_maf(vector<double>& maf) const {
  maf.resize(num_snps, 0);
  const double n_inv = 1.0 / (double)(num_haps);
  for(int j=0; j<num_snps; j++) {
    for(int i=0; i<num_haps; i++) {
      maf[j] += 1.0;
    }
    maf[j] *= n_inv;
    maf[j] = std::min(maf[j], 1.0-maf[j]);
  }  
}

#endif
