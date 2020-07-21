#ifndef PLINK_INTERFACE_CPP
#define PLINK_INTERFACE_CPP

#include "plink_interface.h"

typedef std::uint8_t byte;

using namespace std;

void writeBED_variant(ofstream& file, int j, const vector<chaplotype>& H) {
  // Dimensions
  int num_haps = H.size();
  int num_samples = num_haps/2;
  int num_blocks = ceil((double)num_samples/4.0);

  int sample = 0;
  for(int k=0; k<num_blocks; k++) {
    byte block = 0b00000000;
    for(int i=0; i<4; i++) {
      if(sample<num_samples) {
        bool b0 = H[2*sample][j];
        bool b1 = H[2*sample+1][j];
        if(b0 && b1) {
          block |= 1UL << 2*i;
          block |= 1UL << 2*i+1;
        }
        else if(b0 || b1) {
          block |= 1UL << 2*i+1;
        }
      }
      sample++;
    }
    file.write((char*)(&block), sizeof(block));
  }
}

void writeBED(const string& basename, const vector<chaplotype>& H) {
  // Dimensions
  int num_haps = H.size();
  int num_snps = H[0].size();

  // Create BED file
  string out_filename = basename;
  out_filename.append(".bed");
  ofstream file;
  file.open(out_filename, std::ios_base::binary);
  if (!file.is_open()){
    cerr << "Problem creating the output file: " << out_filename;
    cerr << "Either the directory does not exist or you do not have write permissions." << endl;
    exit(1);
  }

  // Write BED header (magic numbers)
  byte header[] = {0x6c,0x1b,0x01};
  file.write((char*)(header), sizeof(header));

  // Write BED data blocks
  int num_samples = num_haps/2;
  int num_blocks = ceil((double)num_samples/4.0);
  for(int j=0; j<num_snps; j++) {
    writeBED_variant(file, j, H);
  }
  // Close BED file
  file.close();
}

void writeBED(const string& basename, const vector<chaplotype>& H, const vector<chaplotype>& Hk,
              const vector<bool>& swap) {
  // Dimensions
  int num_haps = H.size();
  int num_snps = H[0].size();

  // Create BED file
  string out_filename = basename;
  out_filename.append(".bed");
  ofstream file;
  file.open(out_filename, std::ios_base::binary);
  if (!file.is_open()){
    cerr << "Problem creating the output file: " << out_filename;
    cerr << "Either the directory does not exist or you do not have write permissions." << endl;
    exit(1);
  }

  // Write BED header (magic numbers)
  byte header[] = {0x6c,0x1b,0x01};
  file.write((char*)(header), sizeof(header));

  // Write BED data blocks
  int num_samples = num_haps/2;
  int num_blocks = ceil((double)num_samples/4.0);
  for(int j=0; j<num_snps; j++) {
    // Write original genotypes and knockoffs, possibly swapped
    if(swap[j]) {
      writeBED_variant(file, j, Hk);
      writeBED_variant(file, j, H);
    } else {
      writeBED_variant(file, j, H);
      writeBED_variant(file, j, Hk);
    }
  }
  // Close BED file
  file.close();
}

void writeBIM(const string& basename, const Metadata& metadata, 
              bool include_genotypes, const vector<bool>& swap) {
  // Create output file
  string out_filename = basename;
  out_filename.append(".bim");
  ofstream file(out_filename);
  if (!file.is_open()){
    cerr << "Problem creating the output file: " << out_filename;
    cerr << "Either the directory does not exist or you do not have write permissions." << endl;
    exit(1);
  }

  for(int j=0; j<metadata.num_snps(); j++) {
    if(include_genotypes) {
      if(swap[j]) {
        file << metadata.legend_filter.chr[j] << "\t";
        file << metadata.legend_filter.ID[j] << ".k" << "\t";
        file << 0 << "\t";
        file << metadata.legend_filter.bp[j] << "\t";
        file << metadata.legend_filter.A0[j] << "\t";
        file << metadata.legend_filter.A1[j] << endl;

        file << metadata.legend_filter.chr[j] << "\t";
        file << metadata.legend_filter.ID[j] << "\t";
        file << 0 << "\t";
        file << metadata.legend_filter.bp[j] << "\t";
        file << metadata.legend_filter.A0[j] << "\t";
        file << metadata.legend_filter.A1[j] << endl;
      } else {
        file << metadata.legend_filter.chr[j] << "\t";
        file << metadata.legend_filter.ID[j] << "\t";
        file << 0 << "\t";
        file << metadata.legend_filter.bp[j] << "\t";
        file << metadata.legend_filter.A0[j] << "\t";
        file << metadata.legend_filter.A1[j] << endl;

        file << metadata.legend_filter.chr[j] << "\t";
        file << metadata.legend_filter.ID[j] << ".k" << "\t";
        file << 0 << "\t";
        file << metadata.legend_filter.bp[j] << "\t";
        file << metadata.legend_filter.A0[j] << "\t";
        file << metadata.legend_filter.A1[j] << endl;
      }
    } else {
      file << metadata.legend_filter.chr[j] << "\t";
      file << metadata.legend_filter.ID[j] << ".k" << "\t";
      file << 0 << "\t";
      file << metadata.legend_filter.bp[j] << "\t";
      file << metadata.legend_filter.A0[j] << "\t";
      file << metadata.legend_filter.A1[j] << endl;
    }
  }

  // Close output file
  file.close();
}

void writeFAM(const string& basename, const Metadata& metadata) {
  // Create output file
  string out_filename = basename;
  out_filename.append(".fam");
  ofstream file(out_filename);
  if (!file.is_open()){
    cerr << "Problem creating the output file: " << out_filename;
    cerr << "Either the directory does not exist or you do not have write permissions." << endl;
    exit(1);
  }

  for(int i=0; i<metadata.num_samples(); i++) {
    file << metadata.sample_filter.ID[i] << "\t";
    file << metadata.sample_filter.famID[i] << "\t";
    file << 0 << "\t";
    file << 0 << "\t";
    file << metadata.sample_filter.sex[i] << "\t";
    file << 0 << endl;
  }

  // Close output file
  file.close();
}

void writeHAPS(const string& basename, const vector<chaplotype>& Hk, const Metadata& metadata) {
  int num_haps = Hk.size();
  int num_snps = Hk[0].size();

  string filename = basename;
  filename.append(".haps");
  ofstream file(filename);
  if (!file.is_open()){
    cerr << "Problem creating the output file: " << filename;
    cerr << "Either the directory does not exist or you do not have write permissions." << endl;
  }
  for(int j=0; j<num_snps; j++) {
    // Print variant info
    file << metadata.legend_filter.chr[j] << " ";
    file << metadata.legend_filter.ID[j] << " ";
    file << metadata.legend_filter.bp[j] << " ";
    file << metadata.legend_filter.A0[j] << " ";
    file << metadata.legend_filter.A1[j] << " ";
    // Print variant values for each haplotype
    for(int i=0; i<num_haps; i++) {
      file << Hk[i].get(j);
      if(i+1<num_haps) file <<" ";
    }
    file << endl;
  }
  file.close();
}

void plink::writeSAMPLE(const string& basename, const Metadata& metadata) {
  // Create output file
  string out_filename = basename;
  out_filename.append(".sample");
  ofstream file(out_filename);
  if (!file.is_open()){
    cerr << "Problem creating the output file: " << out_filename;
    cerr << "Either the directory does not exist or you do not have write permissions." << endl;
    exit(1);
  }

  // Write header
  file << "ID_1 ID_2 missing sex" << endl;
  file << "0 0 0 D" << endl;

  // Write body
  for(int i=0; i<metadata.num_samples(); i++) {
    file << metadata.sample_filter.ID[i] << " ";
    file << metadata.sample_filter.famID[i] << " ";
    file << metadata.sample_filter.missing[i] << " ";
    file << metadata.sample_filter.sex[i] << endl;
  }

  // Close output file
  file.close();
}

void plink::write_binary(const string& basename, const Metadata& metadata, const vector<chaplotype>& Hk) {
  assert(Hk.size()>0);  
  int num_snps = Hk[0].size();
  vector<bool> swap_dummy(num_snps,0);
  writeBED(basename, Hk);
  writeBIM(basename, metadata, false, swap_dummy);
  writeFAM(basename, metadata);

  cout << "Output (binary) written to:" << endl;
  cout << "  " << basename << ".bed" << endl;
  cout << "  " << basename << ".bim" << endl;
  cout << "  " << basename << ".fam" << endl;
}

void plink::write_binary(const string& basename, const Metadata& metadata, 
                         const vector<chaplotype>& H, const vector<chaplotype>& Hk, bool random) {

  // Check dimensions
  assert(H.size()>0);
  assert(H.size()%2==0);
  assert(H.size()==Hk.size());
  for(int i=0; i<H.size(); i++) {
    assert(H[i].size()==Hk[i].size());
  }

  // Determine ramdom swaps
  int seed = 2019;
  mt19937_64 rng(seed);
  std::uniform_int_distribution<int> coin_flip(0, 1);
  int num_snps = H[0].size();
  vector<bool> swap(num_snps,0);
  for(int j=0; j<num_snps; j++) {
    swap[j] = coin_flip(rng);
  }

  writeBED(basename, H, Hk, swap);
  writeBIM(basename, metadata, true, swap);
  writeFAM(basename, metadata);

  cout << "Output (binary) written to:" << endl;
  cout << "  " << basename << ".bed" << endl;
  cout << "  " << basename << ".bim" << endl;
  cout << "  " << basename << ".fam" << endl;
}

void plink::write_haps(const string& basename, const Metadata& metadata, const vector<chaplotype>& Hk) {
  // Check dimensions
  assert(Hk.size()>0);  
  assert(Hk.size()%2==0);  

  writeHAPS(basename, Hk, metadata);
  writeSAMPLE(basename, metadata);

  cout << "Output (text) written to:" << endl;
  cout << "  " << basename << ".haps" << endl;
  cout << "  " << basename << ".sample" << endl;
}

#endif
