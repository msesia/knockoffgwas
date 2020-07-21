//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BGEN_UTILS_H
#define _BGEN_UTILS_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "genfile/bgen/bgen.hpp"
#include "utils.h"

using namespace std;

// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct ProbSetter {
  typedef vector< vector< double > > Data;
  ProbSetter(Data* result):
    m_result(result),
    m_sample_i(0)
  {}
		
  // Called once allowing us to set storage.
  void initialise(size_t number_of_samples, size_t number_of_alleles) {
    m_result->clear();
    m_result->resize(number_of_samples);
  }
	
  // If present with this signature, called once after initialise()
  // to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
  // This enables us to set up storage for the data ahead of time.
  void set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries) {
    for(size_t i = 0; i < m_result->size(); ++i) {
      m_result->at(i).reserve(max_entries);
    }
  }
	
  // Called once per sample to determine whether we want data for this sample
  bool set_sample(size_t i) {
    m_sample_i = i;
    // Yes, here we want info for all samples.
    return true;
  }
	
  // Called once per sample to set the number of probabilities that are present.
  void set_number_of_entries(size_t ploidy, size_t number_of_entries, 
                             genfile::OrderType order_type, genfile::ValueType value_type) {
    assert(value_type == genfile::eProbability);
    m_result->at(m_sample_i).resize(number_of_entries);
    m_entry_i = 0;
  }

  // Called once for each genotype (or haplotype) probability per sample.
  void set_value(uint32_t, double value) {
    m_result->at(m_sample_i).at(m_entry_i++) = value;
  }

  // Ditto, but called if data is missing for this sample.
  void set_value(uint32_t, genfile::MissingValue value) {
    // Here we encode missing probabilities with -1
    m_result->at(m_sample_i).at(m_entry_i++) = -1;
  }

  // If present with this signature, called once after all data has been set.
  void finalise() {
    // nothing to do in this implementation.
  }

private:
  Data* m_result;
  size_t m_sample_i;
  size_t m_entry_i;
};

enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 };

class BgenParser {
public:
  BgenParser();
  ~BgenParser();
  BgenParser(string filename);
  ostream& summarise(ostream& o) const;
  size_t number_of_samples() const;
  size_t number_of_variants() const;
  template< typename Setter >
  void get_sample_ids(Setter setter);
  vector<string> get_sample_ids() const;
  bool read_variant(string* chromosome, uint32_t* position, string* rsid, vector< string >* alleles);
  void read_probs(vector< vector<double> >* probs);
  void ignore_probs();

private:
  string m_filename;
  shared_ptr< istream > m_stream;

  // bgen::Context object holds information from the header block,
  // including bgen flags
  genfile::bgen::Context m_context;

  // offset byte from top of bgen file.
  uint32_t m_offset;

  // We keep track of our state in the file.
  // Not strictly necessary for this implentation but makes it clear that
  // calls must be read_variant() followed by read_probs() (or ignore_probs())
  // repeatedly.
  State m_state;

  // If the BGEN file contains samples ids, they will be read here.
  bool m_have_sample_ids;
  vector< string > m_sample_ids;
	
  // Buffers, these are used as working space by bgen implementation.
  vector< genfile::byte_t > m_buffer1, m_buffer2;
};

#endif
