#ifndef _BGEN_HAPS_CPP
#define _BGEN_HAPS_CPP

//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "bgen_utils.h"
#include "bgen_reader.h"

//////////
// Main //
//////////

// This example program reads data from a bgen file specified as the first argument
// and outputs it as a VCF file.
int main( int argc, char** argv ) {
  if( argc != 2 ) {
    std::cerr << "You must supply an argument, the name of the bgen file to process.\n";
    exit(1);
  }
  std::string const filename = argv[1];
  try {
    BgenReader bgenReader(filename);
    
    bgenReader.summarise();

    int num_snps = 30;
    int num_samples = 10;
    vector<chaplotype> H;
    bgenReader.read(num_snps, num_samples, H);

    bgenReader.print(H);

    return 0;
  }
  catch( genfile::bgen::BGenError const& e ) {
    std::cerr << "!! Uh-oh, error parsing bgen file.\n";
    return -1;
  }
}

#endif
