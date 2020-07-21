#ifndef _BGEN_UTILS_CPP
#define _BGEN_UTILS_CPP

#include "bgen_utils.h"

BgenParser::BgenParser() {}

BgenParser::~BgenParser() {}

BgenParser::BgenParser(string filename) {
  
  m_filename = filename;
  m_state = e_NotOpen;
  m_have_sample_ids = false;

  // Open the stream
  m_stream.reset(new ifstream(filename, ifstream::binary));
  if(!*m_stream) {
    throw invalid_argument(filename);
  }
  m_state = e_Open;

  // Read the offset, header, and sample IDs if present.
  genfile::bgen::read_offset(*m_stream, &m_offset);
  genfile::bgen::read_header_block(*m_stream, &m_context);
  if(m_context.flags & genfile::bgen::e_SampleIdentifiers) {
    genfile::bgen::read_sample_identifier_block(*m_stream, m_context, 
                                                [this](string id) { m_sample_ids.push_back(id); });
    m_have_sample_ids = true;
  }
		
  // Jump to the first variant data block.
  m_stream->seekg(m_offset + 4);

  // We keep track of state (though it's not really needed for this implementation.)
  m_state = e_ReadyForVariant;
}

ostream& BgenParser::summarise(ostream& o) const {
  o << "BgenParser: bgen file ("
    << (m_context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout")
    << ", "
    << (m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed") << ")"
    << " with " 
    << m_context.number_of_samples << " " << (m_have_sample_ids ? "named" : "anonymous") << " samples and "
    << m_context.number_of_variants << " variants.\n";
  return o;
}

size_t BgenParser::number_of_samples() const {
  return m_context.number_of_samples;
}

size_t BgenParser::number_of_variants() const {
  return m_context.number_of_variants;
}

// Report the sample IDs in the file using the given vector of strings
vector<string> BgenParser::get_sample_ids() const {
  if(m_have_sample_ids) {
    return(m_sample_ids);
  } else {
    return(vector<string>(m_context.number_of_samples, "unknown"));
  }
}

// Report the sample IDs in the file using the given setter object
// (If there are no sample IDs in the file, we report a dummy identifier).
template< typename Setter >
void BgenParser::get_sample_ids(Setter setter) {
  if(m_have_sample_ids) {
    for(size_t i = 0; i < m_context.number_of_samples; ++i) {
      setter(m_sample_ids[i]);
    }
  } else {
    for(size_t i = 0; i < m_context.number_of_samples; ++i) {
      setter("(unknown_sample_" + std::to_string(i+1) + ")");
    }
  }
}

// Attempt to read identifying information about a variant from the bgen file, returning
// it in the given fields.
// If this method returns true, data was successfully read, and it should be safe to call read_probs()
// or ignore_probs().
// If this method returns false, data was not successfully read indicating the end of the file.
bool BgenParser::read_variant(string* chromosome, uint32_t* position, string* rsid, vector< string >* alleles) {
  assert(m_state == e_ReadyForVariant);
  string SNPID; // read but ignored in this toy implementation
		
  if(
     genfile::bgen::read_snp_identifying_data(*m_stream, m_context, &SNPID, rsid, chromosome, position,
                                              [&alleles](size_t n) { alleles->resize(n); },
                                              [&alleles](size_t i, string const& allele) {alleles->at(i) = allele; }
                                             )
    ) {
    m_state = e_ReadyForProbs;
    return true;
  } else {
    return false;
  }
}

// Read genotype probability data for the SNP just read using read_variant()
// After calling this method it should be safe to call read_variant() to fetch
// the next variant from the file.
void BgenParser::read_probs(vector< vector< double > >* probs) {
  assert(m_state == e_ReadyForProbs);
  ProbSetter setter(probs);
  genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(*m_stream, m_context, setter, &m_buffer1, &m_buffer2);
  m_state = e_ReadyForVariant;
}

// Ignore genotype probability data for the SNP just read using read_variant()
// After calling this method it should be safe to call read_variant()
// to fetch the next variant from the file.
void BgenParser::ignore_probs() {
  genfile::bgen::ignore_genotype_data_block(*m_stream, m_context);
  m_state = e_ReadyForVariant;
}


#endif
