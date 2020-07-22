#include <chrono>
#include "utils.h"
#include "arguments.h"
#include "metadata.h"
#include "generator_master.h"
#include "dataset.h"
#include "kinship.h"
#include "covariates.h"

#define DELIMITER "--------------------------------------------------------------------------------"

using namespace std;

lfile LOG;

void print_header() {
  cout<<"\t" << "+----------------------+" << endl;
  cout<<"\t" << "|                      |" << endl;
  cout<<"\t" << "|  SNPKNOCK2, v0.3     |" << endl;
  cout<<"\t" << "|  July 21, 2020       |" << endl;
  cout<<"\t" << "|  Matteo Sesia        |" << endl;
  cout<<"\t" << "|                      |" << endl;
  cout<<"\t" << "+----------------------+" << endl;
  cout << endl;
  cout << "Copyright (C) 2020 Stanford University." << endl;
  cout << "Distributed under the GNU GPLv3 open source license." << endl << endl;
  cout << "Use --help for more information." << endl;
  cout << endl;
}

int main(int argc, const char * argv[]) {
  print_header();

  // //Log file LOG
  // auto now = std::chrono::system_clock::now();
  // auto in_time_t = std::chrono::system_clock::to_time_t(now);
  // std::stringstream date_ss;
  // date_ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
  // string datestr = date_ss.str();
  // string logfile = "hapknock_" + datestr;
  // logfile += "_" + putils::getRandomID() + ".log";
  // if (!LOG.open(logfile)) {
  //   cerr << "Impossible to open log file[" << logfile << "] check writing permissions!" << endl;
  //   exit(1);
  // }

  // Process input arguments
  Arguments args(argc,argv);
  args.check();

  // Extract parameters
  const int num_chrs = args.num_chrs();
  const int num_threads = args.num_threads();
  const bool compute_kinship = args.get_compute_kinship();
  const bool generate_knockoffs = args.get_generate_knockoffs();
  const bool estimate_hmm = args.get_estimate_hmm();
  const int debug = args.get_debug();
  const int requested_resolution = args.get_resolution();
  const int seed = args.get_seed();
  const int window_size = args.get_window_size();

  // Process metadata
  cout << endl << DELIMITER << endl;
  cout << "Loading metadata" << endl;
  cout << DELIMITER << endl;

  vector<Metadata> metadata;
  for(int chr=0; chr<num_chrs; chr++) {
    string chr_id = std::to_string(chr+1);
    metadata.push_back(Metadata(args.get_filenames(chr), chr_id, window_size));
  }

  // Kinship
  Kinship kinship;
  if(compute_kinship) {
    // Check whether file with covarariates was supplied
    const bool use_covariates = (args.get_pc_file()!="");
    Covariates covariates;

    // Load covariates if necessary
    // FIXME: should use metadata from all chromosomes, because IBD segments may differ
    if(use_covariates) {
      cout << endl << DELIMITER << endl;
      cout << "Loading covariates" << endl;
      cout << DELIMITER << endl;
      covariates = Covariates(args.get_pc_file(), metadata[0]);
    }

    // Print header
    cout << endl << DELIMITER << endl;
    if(use_covariates) {
      cout << "Kinship (using covariates)" << endl;
    } else {
      cout << "Kinship (using only haplotype data)" << endl;
    }
    cout << DELIMITER << endl;

    // Compute kinship
    const int K = args.get_K();
    assert(K>0);
    const int compression = 10;
    cout << "Reached 1" << endl;
    if(use_covariates) {
      // Compute kinship using covariates (principal components)
      kinship = Kinship(metadata, covariates, args.get_output_files(), compression,
                        args.get_cluster_size_min(), args.get_cluster_size_max(), K, num_threads);
    } else {
      // Compute kinship using genetic data only
      kinship = Kinship(metadata, args.get_output_files(), compression,
                        args.get_cluster_size_min(), args.get_cluster_size_max(), K, num_threads);
    }
    cout << "Reached 2" << endl;
    // Save references
    kinship.writeReferences(metadata, args.get_output_files());
  } else {
    if (generate_knockoffs) {
      // Load kinship
      kinship = Kinship(args.get_lref_files(), args.get_ref_files());
    }
  }

  if(generate_knockoffs) {
    for(int chr=0; chr<num_chrs; chr++) {
      cout << endl << DELIMITER << endl;
      cout << "Knockoffs for chromosome " << chr+1 << endl;
      cout << DELIMITER << endl;

      const string out_file = args.get_output_file(chr);
      const string log_file = args.get_log_file(chr);

      const ivector3d& references_local = kinship.get_references(chr);
      const ivector2d& references_global = kinship.get_references_global(chr);
      KnockoffGenerator knockoffs({metadata[chr]}, references_local, references_global,
                                  num_threads, debug, seed, log_file);

      if(args.get_hmm_files().size()>0) {
        // Load HMM
        cout << "Loading HMM" << endl;
        const string hmm_file = args.get_hmm_file(chr);
        knockoffs.load_hmm({hmm_file});
      }
      else {
        const double hmm_rho = args.get_hmm_rho();
        const double hmm_lambda = args.get_hmm_lambda();
        if(estimate_hmm) {
          // Estimate HMM parameters
          knockoffs.fit_HMM(hmm_rho);

        } else {
          cout << "Initializing HMM with user-supplied hyperparameters: ";
          cout << "rho = " << hmm_rho << ", lambda = " << hmm_lambda << "." << endl << endl;
          knockoffs.init_hmm(hmm_rho, hmm_lambda);
        }
        // Save HMM
        knockoffs.writeHMM({out_file});
      }

      // Generate knockoffs at each resolution
      for(int r=0; r<knockoffs.num_partitions(); r++) {
        if(requested_resolution!=-1) {
          if(r!=requested_resolution) continue;
        }
        const string out_file_res = out_file + "_res" + std::to_string(r);

        // Initialize partition
        knockoffs.set_partition(r);

        // Write information
        knockoffs.writeGroups({out_file_res});
        knockoffs.writeWindows({out_file_res});

        // Generate group knockoffs
        knockoffs.generate();

        // Write knockoffs to file
        knockoffs.writeKnockoffs({out_file_res});

        // Write additional information
        knockoffs.writeAncestries({out_file_res});
        // knockoffs.writeZ({out_file_res});
     }
    }
  } else {
    cout << "Skipping knockoff generation." << endl;
  }

  cout << endl << "Finished." << endl << endl;

  return(0);
}
