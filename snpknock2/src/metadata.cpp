#ifndef METADATA_CPP
#define METADATA_CPP

#define DEBUG 0

#include "metadata.h"

Metadata::Metadata(const DataFiles& _df, string _chr_id, int _window_size) {
  df = _df;
  chr_id = _chr_id;
  // Check format
  if(df.format=="bgen") {
    df.sample = df.data + ".sample";
    df.legend = df.data + ".bim";
  } else if (df.format=="haps") {
    df.sample = df.data + ".sample";
    df.legend = df.data + ".legend";
  } else {
    cout << "Error: unknown data format " << df.format << endl;
    exit(1);
  }

  // Load list of samples and variants
  load_sample();
  load_legend();
  // Apply input filters
  apply_filters();
  // Load list of partitions (if available)
  if(df.partitions != "") {
    load_partitions();
  }
  // Load genetic maps
  load_genetic_map();
  // Load list of IBD segments
  int segment_num = 0;
  if(df.ibd != "") {
    segment_num = load_ibd_segments();
  } else {
    unrelated.clear();
    for(int i=0; i<num_haps(); i++) unrelated.insert(i);
    related_samples.resize(num_samples_, vector<int>(0));
    related_families.resize(2*num_samples_, vector<int>(0));
  }

  // Define genomic windows
  define_windows(_window_size);

  // Print information
  cout << "Summary of metadata for chromosome " << chr_id << ":" << endl;
  cout << "  number of samples (after | before filtering) : " << num_samples_;
  cout << " | " << sample.size() << endl;
  cout << "  number of SNPs (after | before filtering)    : " << num_snps_;
  cout << " | " << legend.size() << endl;
  cout << "  number of variant partitions                 : " << partitions.size() << endl;
  if(window_size>0) 
    cout << "  size of genomic windows                      : " << window_size/1e6 << " Mb" << endl;
  else 
    cout << "  size of genomic windows                      : whole-chromosome" << endl;
  cout << "  number of IBD segments                       : " << segment_num << endl;
  cout << endl;

}

Metadata::Metadata(const Metadata& obj) {
  df = obj.df;
  legend = obj.legend;
  legend_filter = obj.legend_filter;
  sample = obj.sample;
  sample_filter = obj.sample_filter;
  partitions = obj.partitions;
  chr_id = obj.get_chr_id();
  num_snps_ = obj.num_snps();
  num_samples_ = obj.num_samples();
  ibd_segments = obj.ibd_segments;
  related_samples = obj.related_samples;
  related_families = obj.related_families;
  unrelated = obj.unrelated;
  related = obj.related;
  windows = obj.windows;
  window_size = obj.window_size;
}

Metadata::~Metadata() {
}

void Metadata::define_windows(int _window_size) {
  //cout << "[DEBUG] in Metadata::define_windows(" << _window_size << ")" << endl;
  window_size = _window_size;
  vector<int> start_points = {0};
  if(_window_size>0) {
    // Assume the last partition is the coarsest one
    const vector<int>& partition = partitions.back();
    int w = 0;
    int bp_cum = 0;  
    for(int j=1; j<num_snps_; j++) {
      const int bp_current = legend_filter.bp[j];
      const int bp_previous = legend_filter.bp[j-1];
      bp_cum += bp_current - bp_previous;
      if((bp_cum > window_size) && (partition[j]!=partition[j-1])) {
        w++;
        bp_cum = 0;
        start_points.push_back(j);
      }
    }
  }
  windows.load(start_points, num_snps_);
  windows.print();
}

void Metadata::load_legend() {
  // Load variant legend file
  cout << "Loading legend from:" << endl;
  cout << "  " << df.legend << endl;
  int line_num = 0;
  string buffer;
  vector<string> tokens;
  ifstream file(df.legend);
  while (std::getline(file, buffer)) {
    boost::split(tokens, buffer, boost::is_any_of(" ,\t"));
    if(df.format == "haps") {
      if(line_num > 0) legend.push_back(tokens);
    } else {
      legend.push_back(tokens);
    }
    line_num++;
  }
  file.close();
  //cout << "Legend size: " << legend.size() << endl;
}

void Metadata::load_sample() {
  // Load sample information file
  cout << "Loading sample information from:" << endl;
  cout << "  " << df.sample << endl;
  sample.clear();
  int line_num = 0;
  string buffer;
  vector<string> tokens;
  ifstream file(df.sample);
  while (std::getline(file, buffer)) {
    boost::split(tokens, buffer, boost::is_any_of(" ,\t"));
    if(line_num > 1) sample.push_back(tokens);
    line_num++;
  }
  file.close();
  //cout << "Sample size: " << sample.size() << endl;
}

int Metadata::load_ibd_segments() {
  cout << "Loading IBD segments from:" << endl;
  cout << "  " << df.ibd << endl;

  // Initialize list of IBD segments
  int num_haps = 2 * num_samples_;
  ibd_segments.clear();
  ibd_segments.resize(num_haps);
  // Initialize list of IBD families and clusters
  vector< vector<int> > related_samples_tmp(num_samples_);
  vector< vector<int> > related_families_tmp(num_haps);
  // Initialize list of unrelated individuals
  unrelated.clear();
  for(int i=0; i<num_haps; i++) unrelated.insert(i);
  // Load sample information file
  int line_num = 1;
  int segment_num = 0;
  string buffer;
  vector<string> tokens;
  ifstream file(df.ibd);
  while (std::getline(file, buffer)) {
    // Skip header line
    if(line_num==1) {
      line_num++;
      continue;
    }

    // Read line
    boost::split(tokens, buffer, boost::is_any_of(" ,\t"));

    // Make sure the format is correct
    if(tokens.size()!=12) {
      cerr << "Error on line " << line_num << endl;
      cerr << "Exiting." << endl;
      exit(-1);
    }

    // // DEBUG (print line)
    // for(int i=0; i<tokens.size(); i++) {
    //   cout << tokens[i] << " ";
    // }
    // cout << endl;

    // Parse segment information
    string chr = tokens[0];
    string ID1 = tokens[1];
    int HID1 = stoi(tokens[2]);
    string ID2 = tokens[3];
    int HID2 = stoi(tokens[4]);
    int bp_min = stoi(tokens[5]);
    int bp_max = stoi(tokens[6]);

    // Sanity checks
    if(ID1!=ID2) {
      if( ((HID1!=0)&&(HID1!=1)) | ((HID2!=0)&&(HID2!=1)) ) {
        cerr << "Skipping IBD segment on line " << line_num << " (invalid individuals)" << endl;
        line_num++;
        continue;
      }
    } else {
      if(! (((HID1==0)&&(HID2==1)) | ((HID1==1)&&(HID2==0))) ) {
        cerr << "Skipping IBD segment on line " << line_num << " (same individual)" << endl;
        line_num++;
        continue;
      }
    }

    // Find index of ID1 and ID2 in the sample list, then convert to haplotype indices
    int s1 = find_sample(ID1);
    int s2 = find_sample(ID2);

    // Skip IBD segments for individuals that we do not want to keep
    if((s1<0)|(s2<0)) {
      //cerr << "Skipping IBD segment on line " << line_num << " (filtered individuals)" << endl;
      line_num++;
      continue;
    }

    int i1 = 2 * s1 + HID1;
    int i2 = 2 * s2 + HID2;

    // Find end points of segment
    int j_min = find_variant(bp_min);
    int j_max = find_variant(bp_max);

    // Sanity check
    if(j_min >= j_max) {
      cerr << "Skipping IBD segment on line " << line_num << " (invalid variant range)" << endl;
      line_num++;
      continue;
    }

    // Add segment to list
    ibd_segments[i1].push_back(IbdSeg({i1, i2}, j_min, j_max, bp_min, bp_max));
    ibd_segments[i2].push_back(IbdSeg({i2, i1}, j_min, j_max, bp_min, bp_max));

    // Increment segment count
    segment_num++;

    // Add individuals to list of IBD-sharing families
    auto iti1 = std::find(related_samples_tmp[s1].begin(), related_samples_tmp[s1].end(), s2);
    if(iti1 == related_samples_tmp[s1].end()) {
      related_samples_tmp[s1].push_back(s2);
    }
    auto iti2 = std::find(related_samples_tmp[s2].begin(), related_samples_tmp[s2].end(), s1);
    if(iti2 == related_samples_tmp[s2].end()) {
      related_samples_tmp[s2].push_back(s1);
    }

    // Add haplotype indices to list of IBD-sharing clusters
    iti1 = std::find(related_families_tmp[i1].begin(), related_families_tmp[i1].end(), i2);
    if(iti1 == related_families_tmp[i1].end()) {
      related_families_tmp[i1].push_back(i2);
    }
    iti2 = std::find(related_families_tmp[i2].begin(), related_families_tmp[i2].end(), i1);
    if(iti2 == related_families_tmp[i2].end()) {
      related_families_tmp[i2].push_back(i1);
    }

    // Update list of related/unrelated haplotypes
    related.insert(i1);
    related.insert(i2);
    unrelated.erase(i1);
    unrelated.erase(i2);

    // Increment line counter
    line_num++;
  }
  file.close();

  cout << "Loaded " << segment_num << " IBD segments." << endl << endl;

  if(DEBUG) print_segments();

  // Find connected components and complete list of families
  connected_components(related_samples_tmp, related_samples);
  connected_components(related_families_tmp, related_families);

  // Sort list of individuals in IBD families
  for(int f=0; f<related_samples.size(); f++) {
    std::sort(related_samples[f].begin(), related_samples[f].end());
  }
  for(int f=0; f<related_families.size(); f++) {
    std::sort(related_families[f].begin(), related_families[f].end());
  }

  if(DEBUG) {
    cout << "[DEBUG] IBD clusters:" << endl;
    for(int i=0; i<related_families.size(); i++) {
      if(related_families[i].size()>1) {
        cout << "\t" << i << " : ";
        for(auto it = related_families[i].begin(); it!=related_families[i].end(); ++it) {
          cout << *it << " ";
        }
        cout << endl;
      }
    }
    cout << "[DEBUG] IBD related samples:" << endl;
    for(int i=0; i<related_samples.size(); i++) {
      if(related_samples[i].size()>1) {
        cout << "\t" << i << " : ";
        for(auto it = related_samples[i].begin(); it!=related_samples[i].end(); ++it) {
          cout << *it << " ";
        }
        cout << endl;
      }
    }
    cout << endl;
  }

  return(segment_num);
}

void Metadata::print_segments() const {
    //cout << "Sample size = " << sample_filter.size() << endl;
    cout << "[DEBUG] IBD segments:" << endl;
    for(int i=0; i<ibd_segments.size(); i++) {
      if(ibd_segments[i].size()>0) {
        cout << "Individual i=" << i << " (ID=" << sample_filter.ID[i/2] << "): " << endl;
        for(auto seg : ibd_segments[i]) {
          cout << "  ";
          seg.print();
          cout << endl;
        }
      }
    }
}

void Metadata::read_matrix(string filename, int num_columns, int skip, vector< vector<string> >& output) const {
  output.clear();
  output.resize(num_columns);
  int line_num = 0;
  string buffer;
  vector<string> tokens;
  ifstream file(filename);
  while (std::getline(file, buffer)) {
    boost::split(tokens, buffer, boost::is_any_of(" ,\t"));
    if(tokens.size()!=num_columns) {
      throw_error_input("in Metadata::read_matrix(), tokens.size()!=num_columns on row "
                        +std::to_string(line_num));
    }
    for(int col=0; col<num_columns; col++) {
      if(line_num >= skip) output[col].push_back(tokens[col]);
    }
    line_num++;
  }
  file.close();
}

void Metadata::apply_filters() {
  // cout << "Applying input filters..." << endl;
  // Filter samples
  if(df.keep != "") {
    // cout << "Reading list of samples to keep... " << flush;
    vector< vector<string> > keep;
    read_matrix(df.keep, 2, 0, keep);
    // cout << keep[0].size() << " rows." << endl;
    // cout << "Applying filter... " << flush;
    sample.filter(keep[0]);
    // cout << "done." << endl;
    // cout << "Extracting subsample... " << flush;
    sample_filter = sample.extract_filtered();
    // cout << "done." << endl;
  } else {
    sample_filter = sample;
  }

  // Count number of samples that passed filter
  num_samples_ = sample_filter.size();
  // cout << "Number of samples that passed filter: " << num_samples_ << endl;

  // Filter variants
  if(df.extract != "") {
    // cout << "Reading list of variants to extract... " << flush;
    vector< vector<string> > extract;
    read_matrix(df.extract, 1, 0, extract);
    // cout << extract[0].size() << " rows." << endl;
    // cout << "Applying filter... " << flush;
    legend.filter(extract[0]);
    // cout << "done." << endl;
    // cout << "Extracting subsample... " << flush;
    legend_filter = legend.extract_filtered();
    // cout << "done." << endl;
  } else {
    legend_filter = legend;
  }
  // Count number of variants that passed filter
  num_snps_ = legend_filter.size();
  // cout << "Number of SNPs that passed filter: " << num_snps_ << endl;
}

void Metadata::interpolate(const vector<double>& x, const vector<double>& y,
                           const vector<double>& x_new, vector<double>& y_new) const {

   alglib::real_1d_array ax, ay;
   ax.setcontent(x.size(), &(x[0]));
   ay.setcontent(y.size(), &(y[0]));
   alglib::spline1dinterpolant spline;
   alglib::spline1dbuildmonotone(ax, ay, spline);

   y_new.clear();
   y_new.resize(x_new.size());
   for(int i=0; i<x_new.size(); i++) {
     y_new[i] = alglib::spline1dcalc(spline,x_new[i]);
   }
}

void Metadata::load_genetic_map() {
  if(df.map=="") {
    vector<double> leg_cm(legend_filter.bp.begin(), legend_filter.bp.end());
    legend_filter.cm = leg_cm;
  } else {
    // cout << "Loading genetic map from " << df.map << endl;

    // Load map from file
    vector<double> map_bp;
    vector<double> map_cm;
    string buffer;
    vector <string> tokens;
    int line = 0;
    line = -1;
    ifile fd_map(df.map);
    while (getline(fd_map, buffer, '\n')) {
      if(line>=0) {
        // Divide line into columns
        boost::split(tokens, buffer, boost::is_any_of(" ,\t"));
        assert(tokens.size() == 4);
        // Read variant physical position (bp) and map (cM)
        map_bp.push_back(std::stod(tokens[1]));
        map_cm.push_back(std::stod(tokens[3]));
      }
      line ++;
    }
    fd_map.close();

    // Interpolate genetic positions
    vector<double> leg_bp(legend_filter.bp.begin(), legend_filter.bp.end());
    vector<double> leg_cm(legend_filter.size());
    interpolate(map_bp, map_cm, leg_bp, leg_cm);
    legend_filter.cm = leg_cm;

    // cout << endl << "[DEBUG] legend_filter (" << legend_filter.size() << "):" << endl;
    // for(int j=0; j<10; j++) {
    //   cout << legend_filter.ID[j] << " " << legend_filter.bp[j] << " " << legend_filter.cm[j] << " " << endl;
    // }

    //cout << "Loaded genetic map for " << legend_filter.size() << " variants." << endl;
  }
}

int Metadata::find_sample(string _id) const {
  // Search subject by ID
  auto it = std::find(sample_filter.ID.begin(), sample_filter.ID.end(), _id);
  if(it != sample_filter.ID.end()) {
    int j = std::distance(sample_filter.ID.begin(), it);
    return(j);
  } else {
    // cerr << "Error: sample " << _id << " not found." << endl;
    return(-1);
  }
}

int Metadata::find_variant(string _id) const {
  // Search variant by ID
  auto it = std::find(legend_filter.ID.begin(), legend_filter.ID.end(), _id);
  if(it != legend_filter.ID.end()) {
    int j = std::distance(legend_filter.ID.begin(), it);
    return(j);
  } else {
    cerr << "Error: variant " << _id << " not found." << endl;
    return(-1);
  }
}

int Metadata::find_variant(int _bp) const {
  // Search variant by ID
  auto it = std::lower_bound(legend_filter.bp.begin(), legend_filter.bp.end(), _bp);
  if(it != legend_filter.bp.end()) {
    int j = std::distance(legend_filter.bp.begin(), it);
    return(j);
  } else {
    cerr << "Error: variant at position " << _bp << " not found." << endl;
    return(-1);
  }
}

void Metadata::load_partitions() {
  cout << "Loading partitions from:" << endl << "  " << df.partitions << endl;
  string buffer;
  vector < string > tokens;
  ifile fd(df.partitions);

  // Read file header and determine number of partitions
  getline(fd, buffer, '\n');
  boost::split(tokens, buffer, boost::is_any_of(" ,\t"));
  int num_partitions = tokens.size()-1;
  // cout << "num_partitions = " << num_partitions << endl;
  assert(num_partitions>=1);
  partitions.resize(num_partitions);
  for(int r=0; r<num_partitions; r++) {
    partitions[r].resize(num_snps_,-1);
  }

  // Read file body
  int line = 0;
  while (getline(fd, buffer, '\n')) {
    boost::split(tokens, buffer, boost::is_any_of(" ,\t"));
    assert(tokens.size() == num_partitions+1);
    // Make sure that the SNP id matches that in the post-filter legend
    assert(tokens[0] == legend_filter.ID.at(line));
    // Store group assignments for each partition
    for(int r=0; r<num_partitions; r++) {
      partitions[r][line] = std::stoi(tokens[1+r]);
    }
    line++;
  }
  fd.close();
  assert(partitions.size()>0);
  // cout << "Read partition file." << endl;

  // Make sure we have found partition info for all variants that passed filters
  assert(line==num_snps_);

  // Make sure that group numbering starts from 0
  for(int r=0; r<partitions.size(); r++) {
    int new_group = 0;
    int old_group = 0;
    for(int j=0; j<num_snps_; j++) {
      if(j>0) {
        if(partitions[r][j]!=old_group) new_group++;
      }
      old_group = partitions[r][j];
      partitions[r][j] = new_group;
    }
  }
 //cout << "Loaded " << partitions.size() << " partitions for " << partitions[0].size() << " variants." << endl;
}

int Metadata::num_samples() const {
  return(num_samples_);
}

int Metadata::num_haps() const {
  return(2*num_samples_);
}

int Metadata::num_snps() const {
  return(num_snps_);
}

string Metadata::get_chr_id() const {
  return(chr_id);
}

Metadata Metadata::thin(int compression) const {
  assert(compression>=1);

  // Clone this object
  Metadata new_metadata = *this;

  // Thin out the SNPs in cloned object
  vector<string> new_snp_ids;
  for(int j=0; j<num_snps_; j++) {
    if(j % compression == 0) {
      new_snp_ids.push_back(legend_filter.ID[j]);
    }
  }
  new_metadata.legend.filter(new_snp_ids);
  new_metadata.legend_filter = new_metadata.legend.extract_filtered();
  new_metadata.num_snps_ = new_metadata.legend_filter.size();

  return(new_metadata);
}

Sample::Sample() {
}

Sample::Sample(const Sample& obj) {
  ID = obj.ID;
  famID = obj.famID;
  missing = obj.missing;
  sex = obj.sex;
  keep = obj.keep;
}

Sample::~Sample() {
}

void Sample::clear() {
  ID.clear();
  famID.clear();
  missing.clear();
  sex.clear();
  keep.clear();
}

void Sample::push_back(vector<string> row_) {
  assert((row_.size()==4)|(row_.size()==5));
  ID.push_back(row_[0]);
  famID.push_back(row_[1]);
  missing.push_back((bool)std::atoi(row_[2].c_str()));
  sex.push_back(std::atoi(row_[3].c_str()));
  keep.push_back(true);
}

vector<string> Sample::row(int idx) const {
  vector<string> output;
  output.push_back(ID.at(idx));
  output.push_back(famID.at(idx));
  output.push_back(std::to_string(missing.at(idx)));
  output.push_back(std::to_string(sex.at(idx)));
  return(output);
}

int Sample::size() const {
  return(ID.size());
}

void Sample::filter(const vector<string>& id_list) {
  vector<int> v;
  match_indices(ID, id_list, v);
  std::fill(keep.begin(),keep.end(),false);
  for(int i=0; i<v.size(); i++) {
    keep[v[i]] = true;
  }
}

Sample Sample::extract_filtered() const {
  Sample new_sample;
  for(int j=0; j<size(); j++) {
    if(keep[j]) {
      new_sample.push_back(row(j));
    }
  }
  return(new_sample);
}

Legend::Legend() {
}

Legend::Legend(const Legend& obj) {
  chr = obj.chr;
  ID = obj.ID;
  bp = obj.bp;
  cm = obj.cm;
  A0 = obj.A0;
  A1 = obj.A1;
  extract = obj.extract;
}

Legend::~Legend() {
}

void Legend::clear() {
  chr.clear();
  ID.clear();
  bp.clear();
  cm.clear();
  A0.clear();
  A1.clear();
  extract.clear();
}

void Legend::push_back(vector<string> row_) {
  assert((row_.size()==4)||(row_.size()==5)||(row_.size()==6));
  if(row_.size()==4) {
    // legend format, without cM
    chr.push_back("?");
    ID.push_back(row_[0]);
    bp.push_back(std::stoi(row_[1]));
    cm.push_back(-1);
    A0.push_back(row_[2]);
    A1.push_back(row_[3]);
    extract.push_back(true);
  } else if(row_.size()==5) {
    // legend format, with cM
    chr.push_back("?");
    ID.push_back(row_[0]);
    bp.push_back(std::stoi(row_[1]));
    cm.push_back(std::stod(row_[2]));
    A0.push_back(row_[3]);
    A1.push_back(row_[4]);
    extract.push_back(true);
  } else if(row_.size()==6) {
    // BIM format
    chr.push_back(row_[0]);
    ID.push_back(row_[1]);
    bp.push_back(std::stoi(row_[3]));
    cm.push_back(std::stod(row_[2]));
    A0.push_back(row_[4]);
    A1.push_back(row_[5]);
    extract.push_back(true);
  }
}

vector<string> Legend::row(int idx) const {
  vector<string> output;
  output.push_back(chr.at(idx));
  output.push_back(ID.at(idx));
  output.push_back(std::to_string(cm.at(idx)));
  output.push_back(std::to_string(bp.at(idx)));
  output.push_back(A0.at(idx));
  output.push_back(A1.at(idx));
  return(output);
}

int Legend::size() const {
  return ID.size();
}

void Legend::filter(const vector<string>& id_list) {
  // new code
  vector<int> v;
  match_indices(ID, id_list, v);
  std::fill(extract.begin(),extract.end(),false);
  for(int i=0; i<v.size(); i++) {
    extract[v[i]] = true;
  }

  // int extracted = 0;
  // std::fill(extract.begin(),extract.end(),false);
  // for(auto snp_it = id_list.begin(); snp_it != id_list.end(); ++snp_it) {
  //    auto found_it = std::find(ID.begin(), ID.end(), *snp_it);
  //    if (found_it != std::end(ID)) {
  //      int found_idx = std::distance(ID.begin(), found_it);
  //      extract[found_idx] = true;
  //      extracted++;
  //    }
  // }

}

Legend Legend::extract_filtered() const {
  Legend new_legend;
  for(int j=0; j<size(); j++) {
    if(extract[j]) {
      new_legend.push_back(row(j));
    }
  }
  return(new_legend);
}

void Legend::print() const {
  cout << "Printing legend (size " << size() << ")" << endl;
  for(int j=0; j<std::min(10,size()); j++) {
    cout << ID[j] << " " << bp[j] << " " << cm[j] << endl;
  }
}

#endif
