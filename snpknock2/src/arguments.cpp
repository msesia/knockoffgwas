#include "arguments.h"

using namespace std;

//constructor definition
Arguments::Arguments(int arg_num, const char* argslist[]) {
  n_threads = 1;
  K = K_default;
  cluster_size_min = cluster_size_min_default;
  cluster_size_max = cluster_size_max_default;
  compute_kinship = false;
  estimate_hmm = false;
  generate_knockoffs = false;
  keep_file = "";
  pc_file = "";
  debug = 0;
  hmm_rho = -1;
  hmm_lambda = -1;
  resolution = -1;
  seed = 2020;
  window_size = 0;

  parse_args(arg_num,argslist);

  num_chrs_ = data_files.size();
}

Arguments::Arguments() {
  //cout << "Default constructor called" << endl;
}

int Arguments::num_chrs() const {
  return(num_chrs_);
}

void Arguments::print_help() const {
  typedef tuple<string,string,string> argdef;
  cout << endl << "Usage manual:" << endl;
  list< argdef > options = {
    argdef("help", "", "Print this help page and exit."),
    argdef("bgen", "fileprefix", "Prefix (no extension) of haplotype input file in BGEN format.*"),
    argdef("haps", "fileprefix", "Prefix (no extension) of haplotype input file in HAPS format.*"),
    argdef("map", "filename", "File containing map with genetic distances.*"),
    argdef("part", "filename", "File containing genome partitions for knockoff generation.*"),
    argdef("extract", "filename", "File containing list of selected variant IDs; other variants will be ignored."),
    argdef("sample", "filename", "File containing list of selected sample IDs; other samples will be ignored."),
    argdef("hmm", "filename", "File containing estimated HMM parameters.*"),
    argdef("hmm-rho", "value", "Initial value of HMM recombination rate parameter."),
    argdef("hmm-lambda", "value", "Initial value of HMM mutation rate parameter."),
    argdef("ref", "filename", "File containing pre-computed list of haplotype references.*"),
    argdef("lref", "filename", "File containing pre-computed list of local haplotype references.*"),
    argdef("ibd", "filename", "File containing pre-computed list of IBD segments."),
    argdef("pc", "filename", "File containing genetic PCs. These can be optionally used to define references."),
    argdef("n_threads", "value", "Number of threads available for computations."),
    argdef("resolution", "value", "Index of partition from 'part' file. Default: all partitions will be processed."),
    argdef("windows", "value", "Width (# bases) for windows defining local reference. Default: no local windows."),
    argdef("cluster_size_min", "value", "Smallest allowed size of haplotype clusters (2-means recursive). Default: 1000."),
    argdef("cluster_size_max", "value", "Largest allowed size of haplotype clusters (2-means recursive). Default: 5000."),
    argdef("K", "value", "Number of reference for each haplotype mosaic. Default: 100."),
    argdef("seed", "value", "Seed for pseudo-random number generation. Default: 2020."),
    argdef("compute-references", "", "Request computation of haplotype references."),
    argdef("estimate-hmm", "", "Request estimation of HMM parameters with EM algorithm."),
    argdef("generate-knockoffs", "", "Request knockoff generation."),
    argdef("out", "fileprefix", "Prefix of output data files."),
    argdef("log", "fileprefix", "Prefix of output log files."),
    argdef("debug", "", "Print debugging messages. (For development only)"),    
  };

  for(auto cmd : options) {
    string cmd_name = std::get<0>(cmd);
    string cmd_value = std::get<1>(cmd);
    string cmd_description = std::get<2>(cmd);
    cout << "  " << std::internal << "--" << std::left << std::setw(18) << cmd_name << " ";    
    if(cmd_value == "") {
      cout << std::left << std::setw(15) << "[NONE]";
    } else {
      cmd_value = "[" + cmd_value + "]";
      cout << std::left << std::setw(15) << cmd_value;
    }
    cout << cmd_description << endl;
  }

  cout << endl << "* This command accepts a shortcut to process multiple chromosome files simultaneously." << endl;
  cout << "  E.g., 'filename{1:22}' will be interpreted as a list of filenames filename1,filename2,filename3,...,filename22.";
  cout << endl << endl;
}

void Arguments::parse_args(int arg_num, const char* argslist[]) {
  vector<string> args;
  if(arg_num==1){
    exit(1);
  }
  else {
    for(int i=1;i<arg_num;i++){
      args.push_back(string(argslist[i]));
    }
    cout << "Command line arguments:" << endl;
    for(vector<string>::iterator sit=args.begin();sit!=args.end();sit++) {
      //cout << *sit << endl;
      string name = *sit;
      if(name=="--help") {
        cout << "  " << name << endl;
        print_help();
        exit(0);
      }
      if(name=="--haps") {
        string value = *(++sit);
        parse_wildcard_str(value,"hap");
        data_format = "haps";
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--bgen") {
        string value = *(++sit);
        parse_wildcard_str(value,"hap");
        data_format = "bgen";
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--extract") {
        string value = *(++sit);
        parse_wildcard_str(value,"extract");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--keep") {
        string value = *(++sit);
        parse_file_str(value,"k");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--map") {
        string value = *(++sit);
        parse_wildcard_str(value,"map");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--part") {
        string value = *(++sit);
        parse_wildcard_str(value,"part");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--pc") {
        string value = *(++sit);
        parse_file_str(value,"pc");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--hmm") {
        string value = *(++sit);
        parse_wildcard_str(value,"hmm");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--hmm-rho") {
        string value = *(++sit);
        parse_hmm_rho_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--hmm-lambda") {
        string value = *(++sit);
        parse_hmm_lambda_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--ref") {
        string value = *(++sit);
        parse_wildcard_str(value,"ref");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--lref") {
        string value = *(++sit);
        parse_wildcard_str(value,"lref");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--ibd") {
        string value = *(++sit);
        parse_wildcard_str(value,"ibd");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--out") {
        string value = *(++sit);
        parse_wildcard_str(value,"out");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--log") {
        string value = *(++sit);
        parse_wildcard_str(value,"log");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name== "--n_threads") {
        string value = *(++sit);
        parse_n_threads_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name== "--resolution") {
        string value = *(++sit);
        parse_resolution_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name== "--windows") {
        string value = *(++sit);
        parse_window_size_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name== "--debug") {
        string value = *(++sit);
        parse_debug_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--cluster_size_min") {
        string value = *(++sit);
        parse_cluster_size_option(value, "min");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--cluster_size_max") {
        string value = *(++sit);
        parse_cluster_size_option(value, "max");
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--K") {
        string value = *(++sit);
        parse_K_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--seed") {
        string value = *(++sit);
        parse_seed_option(value);
        cout << "  " << name << " " << value << endl;
        continue;
      }
      if(name=="--compute-references") {
        cout << "  " << name << endl;
        compute_kinship = true;
        operations.insert("--compute-references");
        continue;
      }
      if(name=="--estimate-hmm") {
        cout << "  " << name << endl;
        estimate_hmm = true;
        operations.insert("--estimate-hmm");
        continue;
      }
      if(name=="--generate-knockoffs") {
        cout << "  " << name << endl;
        generate_knockoffs = true;
        operations.insert("--generate_knockoffs");
        continue;
      }
    }
    cout << endl;
  }
}

void Arguments::check() const {
  // Print list of tasks
  if(operations.size()>0) {
    cout << "Requested operations:" << endl;
    for(auto it = operations.begin(); it != operations.end(); ++it) {
      cout << "  " << *(it) << endl;
    }
    cout << endl;
  } else {
    cout << "Nothing to do." << endl;
    return;
  }

  // Store warnings and errors
  vector<string> options_warnings;
  vector<string> options_errors;

  // Sanity checks on number of required input files
  if(data_files.size() == 0) {
    options_errors.push_back("Missing required argument --haps or --bgen.");
  }
  if(out_files.size() == 0) {
    options_errors.push_back("Missing required argument --out.");
  }
  if(out_files.size() != data_files.size()) {
    options_errors.push_back("Arguments (--haps or --bgen) and --out should have the same length.");
  }
  if(map_files.size() > 0) {
    if(data_files.size() != map_files.size()) {
      options_errors.push_back("Arguments (--haps or --bgen) and --map should have the same length if --map is specified.");
    }
  } else {
      options_warnings.push_back("Argument --map is recommended.");
  }
  if(n_threads <= 0) {
    options_errors.push_back("Argument --n_threads should be positive.");
  }

  // Make sure parameters are compatible
  if(compute_kinship) {
    if(K <= 0) {
      options_errors.push_back("Argument --K should be positive.");
    }
    if(cluster_size_min <= 0) {
      options_errors.push_back("Argument --cluster_size_min should be positive.");
    }
    if(cluster_size_max <= 0) {
      options_errors.push_back("Argument --cluster_size_max should be positive.");
    }
    if(ref_files.size()>0) {
      options_warnings.push_back("Argument --ref will be ignored because --compute-references is specified.");
    }
    if(lref_files.size()>0) {
      options_warnings.push_back("Argument --lref will be ignored because --compute-references is specified.");
    }
  } else {
    if(K!=K_default) {
      options_warnings.push_back("Argument --K will be ignored because --ref is specified.");
    }
    if(cluster_size_min!=cluster_size_min_default) {
      options_warnings.push_back("Argument --cluster_size_min will be ignored because --ref is specified.");
    }
    if(cluster_size_max!=cluster_size_max_default) {
      options_warnings.push_back("Argument --cluster_size_max will be ignored because --ref is specified.");
    }
    if(ref_files.size()!=data_files.size()) {
      options_errors.push_back("Argument --ref should be provided (for each chromosome) if --compute-references is not requested.");
    }
  }
  if(lref_files.size()>0) {
    if(ref_files.size()!=lref_files.size()) {
      options_errors.push_back("Argument --lref should have the same length as --ref, if specified.");
    }
  }

  if(generate_knockoffs) {
    if(data_files.size() != part_files.size()) {
      options_errors.push_back("Arguments (--haps or --bgen) and --part should have the same length if --generate_knockoffs is specified.");
    }
  }

  // Make sure files exist

  // HAPS or BGEN files
  if (data_format=="bgen") {
    for(int i=0; i<data_files.size(); i++) {
      string prefix = data_files[i];
      check_file(prefix + ".bgen");
      check_file(prefix + ".bim");
      check_file(prefix + ".sample");
    }
  } else if (data_format=="haps"){
    for(int i=0; i<data_files.size(); i++) {
      string prefix = data_files[i];
      check_file(prefix + ".haps");
      check_file(prefix + ".legend");
      check_file(prefix + ".sample");
    }
  } else {
    std::cerr << "Error: unknown input data format " << data_format << std::endl;
    exit(1);
  }

  // Partition files (optional)
  if(part_files.size()>0) {
    for(int i=0; i<part_files.size(); i++) check_file(part_files[i]);
  }
  // Map files (optional)
  if(map_files.size()>0) {
    for(int i=0; i<map_files.size(); i++) check_file(map_files[i]);
  }
  // HMM files (optional)
  if(hmm_files.size()>0) {
    for(int i=0; i<hmm_files.size(); i++) check_file(hmm_files[i]);
  }
  // Reference files (optional)
  if(!compute_kinship && ref_files.size()>0) {
    for(int i=0; i<ref_files.size(); i++) check_file(ref_files[i]);
  }
  if(!compute_kinship && lref_files.size()>0) {
    for(int i=0; i<lref_files.size(); i++) check_file(lref_files[i]);
  }
  // Keep file (optional)
  if(keep_file!="") check_file(keep_file);
  // Extract files (optional)
  if(extract_files.size()>0) {
    for(int i=0; i<extract_files.size(); i++) check_file(extract_files[i]);
  }
  // PC file (optional)
  if(pc_file!="") check_file(pc_file);
  // IBD files (optional)
  if(ibd_files.size()>0) {
    for(int i=0; i<ibd_files.size(); i++) check_file(ibd_files[i]);
  }

  // Print errors (if any) and exit
  if(options_errors.size()>0) {
    cout << "Command line argument parsing errors:" << endl;
    for(auto it = options_errors.begin(); it != options_errors.end(); ++it) {
      cout << "  " << *(it) << endl;
    }
    cout << endl << "Exiting." << endl;
    exit(1);
  }

  // Print warnings
  if(options_warnings.size()>0) {
    cout << "Command line argument parsing warnings:" << endl;
    for(auto it = options_warnings.begin(); it != options_warnings.end(); ++it) {
      cout << "  " << *(it) << endl;
    }
    cout << endl;
  }
}

void Arguments::check_file(string input_file) const {
  if(!futils::isFile(input_file)) {
    cout << endl << "Problem opening '" << input_file << "'.\n";
    cout << "Make sure that this file exists and that you have read permission." << endl;
    exit(1);
  }
}

void Arguments::parse_K_option(string K_str) {
  K = atoi(K_str.c_str());
}

void Arguments::parse_seed_option(string seed_str) {
  seed = atoi(seed_str.c_str());
}

void Arguments::parse_hmm_rho_option(string input) {
  hmm_rho = atof(input.c_str());
}

void Arguments::parse_hmm_lambda_option(string input) {
  hmm_lambda = atof(input.c_str());
}

void Arguments::parse_cluster_size_option(string str, string state) {
  if(state=="min") cluster_size_min = atoi(str.c_str());
  if(state=="max") cluster_size_max = atoi(str.c_str());
}

void Arguments::parse_n_threads_option(string str) {
  n_threads = atoi(str.c_str());
}

void Arguments::parse_resolution_option(string str) {
  resolution = atoi(str.c_str());
}

void Arguments::parse_window_size_option(string str) {
  window_size = atof(str.c_str());
}

void Arguments::parse_debug_option(string str) {
  debug = atoi(str.c_str());
}

void Arguments::parse_wildcard_str(string input, string state) {
  bool parse_expression = false;
  auto open_brace = input.find("{");
  if(open_brace != std::string::npos) {
    parse_expression = true;
  }
  auto colon = input.find(":", open_brace+1);
  if(colon != std::string::npos) {
  } else {
    parse_expression = false;
  }
  auto close_brace = input.find("}", colon+1);
  if(close_brace != std::string::npos) {
  } else {
    parse_expression = false;
  }

  vector<string> filenames;

  if(parse_expression) {
    string left_str = std::string(input.begin()+open_brace+1, input.begin()+colon);
    string right_str = std::string(input.begin()+colon+1, input.begin()+close_brace);
    int left = std::stoi(left_str);
    int right = std::stoi(right_str);

    string prefix = std::string(input.begin(), input.begin()+open_brace);
    string suffix = std::string(input.begin()+close_brace+1, input.end());

    for(int chr=left; chr<=right; chr++) {
      filenames.push_back(prefix + std::to_string(chr) + suffix);
    }

  } else {
    filenames.push_back(input);
  }

  // cout << "Parsed input: " << endl;
  for(int i=0; i<filenames.size(); i++) {
    if(state=="bgen") data_files.push_back(filenames[i]);
    if(state=="hap") data_files.push_back(filenames[i]);
    if(state=="extract") extract_files.push_back(filenames[i]);
    if(state=="map") map_files.push_back(filenames[i]);
    if(state=="part") part_files.push_back(filenames[i]);
    if(state=="hmm") hmm_files.push_back(filenames[i]);
    if(state=="ref") ref_files.push_back(filenames[i]);
    if(state=="lref") lref_files.push_back(filenames[i]);
    if(state=="ibd") ibd_files.push_back(filenames[i]);
    if(state=="out") out_files.push_back(filenames[i]);
    if(state=="log") log_files.push_back(filenames[i]);
    // cout << "  --" << state << " " << filenames[i] << endl;
  }
  // cout << endl;
}

void Arguments::parse_file_str(string file_name, string state) {
  //int found = 0;
  string type;
  string temp = file_name;
  int ext_start =0;
  int length_of_name = file_name.length();
  //cout << "parsing " << file_name << endl;
  if((ext_start=file_name.rfind(".gz"))!=string::npos) {
    if(length_of_name-ext_start==3)
      temp = file_name.substr(0,ext_start);
    //cout << "temp: " << temp << endl;
  }
  if((ext_start=temp.rfind("."))!=string::npos) {
    type = temp.substr(ext_start+1);
    //found =1;
  }

  if(state=="k") {
    keep_file = file_name;
  }
  if(state=="pc") {
    pc_file = file_name;
  }
  if(state=="h") {
    data_files.push_back(file_name);
  }
  if(state=="m") {
    map_files.push_back(file_name);
  }

}

string Arguments::get_filename(string filetype) const {
  std::cerr << "Unknown file type type " << filetype << std::endl;
  return("");
}

int Arguments::get_resolution() const {
  return(resolution);
}

DataFiles Arguments::get_filenames(int chr) const {
  DataFiles df;
  df.format = data_format;
  df.data = data_files.at(chr);
  df.keep = keep_file;
  if(map_files.size()>0) df.map = map_files.at(chr);
  if(part_files.size()>0) df.partitions = part_files.at(chr);
  if(extract_files.size()>0) df.extract = extract_files.at(chr);
  if(ibd_files.size()>0) df.ibd = ibd_files.at(chr);
  return(df);
}

vector<DataFiles> Arguments::get_filenames() const {
  vector<DataFiles> df;
  for(int chr=0; chr<num_chrs_; chr++) {
    df.push_back(get_filenames(chr));
  }
  return(df);
}

string Arguments::get_hmm_file(int chr) const {
  return(hmm_files[chr]);
}

vector<string> Arguments::get_hmm_files() const {
  return(hmm_files);
}

string Arguments::get_ibd_file(int chr) const {
  return(ibd_files[chr]);
}

vector<string> Arguments::get_ibd_files() const {
  return(ibd_files);
}

string Arguments::get_ref_file(int chr) const {
  return(ref_files[chr]);
}

vector<string> Arguments::get_ref_files() const {
  return(ref_files);
}

string Arguments::get_lref_file(int chr) const {
  if(lref_files.size() == ref_files.size()) {
    return(lref_files[chr]);
  } else {
    return(get_ref_file(chr));
  }
}

vector<string> Arguments::get_lref_files() const {
  if(lref_files.size() == ref_files.size()) {
    return(lref_files);
  } else {
    return(get_ref_files());
  }
}

string Arguments::get_output_file(int chr) const {
  return(out_files[chr]);
}

vector<string> Arguments::get_output_files() const {
  return(out_files);
}

string Arguments::get_log_file(int chr) const {
  if(log_files.size()!=out_files.size()) return("tmp/hapknock_log");
  else return(log_files[chr]);
}

vector<string> Arguments::get_log_files() const {
  if(log_files.size()!=out_files.size()) {
    vector<string> emptystrings(out_files.size(), "");
    return(emptystrings);
  }
  else return(log_files);
}

int Arguments::get_cluster_size_min() const {
  return(cluster_size_min);
}

int Arguments::get_cluster_size_max() const {
  return(cluster_size_max);
}

int Arguments::num_threads() const {
  return(n_threads);
}

int Arguments::get_debug() const {
  return(debug);
}

int Arguments::get_K() const {
  return(K);
}

int Arguments::get_seed() const {
  return(seed);
}

double Arguments::get_hmm_rho() const {
  return(hmm_rho);
}

double Arguments::get_hmm_lambda() const {
  return(hmm_lambda);
}

bool Arguments::get_compute_kinship() const {
  return(compute_kinship);
}

bool Arguments::get_estimate_hmm() const {
  return(estimate_hmm);
}

bool Arguments::get_generate_knockoffs() const {
  return(generate_knockoffs);
}

string Arguments::get_pc_file() const {
  return(pc_file);
}

double Arguments::get_window_size() const {
  return(window_size);
}
