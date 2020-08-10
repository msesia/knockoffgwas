#ifndef _UTILS_H
#define _UTILS_H

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.14159265358979323846

#define LOW_POS_DOUBLE 1e-300
#define BIG_POS_DOUBLE 1e300
#define LOW_NEG_DOUBLE -1e-300
#define BIG_NEG_DOUBLE -1e300
#define LOW_POS_FLOAT 1e-30
#define BIG_POS_FLOAT 1e30
#define LOW_NEG_FLOAT -1e-30
#define BIG_NEG_FLOAT -1e30
#define BIG_POS_INT 1000000000
#define BIG_NEG_INT -1000000000

#include <string>
#include <vector>
#include <queue>
#include <map>
#include <bitset>
#include <list>
#include <set>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>

//#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/pending/integer_log2.hpp>

using namespace std;
namespace bio = boost::iostreams;
namespace bpo = boost::program_options;
namespace bid = boost::uuids;

/******************************************************/
/*                  UTILS STATISTICS                  */
/******************************************************/
namespace putils {
  void initRandom(long s);
  double getRandom();
  string getRandomID();
  int getRandom(int);
  long getSeed();
  void normalise(vector < double > & v);
  int sample(vector< double > & v, double sum);
  double entropy(vector < double > & v);
  double KLdistance(vector < double > & P, vector < double > & Q);
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
  int max(vector < double > & v);
  int max(vector < int > & v);
  void findUniqueSet(vector < bool > & B, vector < int > & U);
  void decompose(int min, vector < vector < int > > & B, vector < vector < vector < int > > > & BB);
  int checkDuo (int pa1, int pa2, int ca1, int ca2);
  int checkTrio (int fa1, int fa2, int ma1, int ma2, int ca1, int ca2);
};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
  int tokenize(string &, vector < string > &);
  int tokenize(string &, vector < string > &, int);
  string int2str(int n);
  string int2str(vector < int > & v);
  string long2str(long int n);
  string double2str(double n, int prc = 4);
  string double2str(vector < double > &v, int prc = 4);
  string bool2str(vector<bool> & v);
  string date2str(time_t * t, string format);
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
  bool isFile(string f);
  bool createFile(string f);
  string extensionFile(string & filename);
  void bool2binary(vector < bool > & V, ostream &fd);
  bool binary2bool(vector < bool > & V, istream & fd);
};


/******************************************************/
/*                  EXCEPTIONS                        */
/******************************************************/
class myException : public exception {
public:
  explicit myException(std::string msg) : msg_(msg) {}

  virtual ~myException() throw() {}

  virtual const char* what() const throw() {
    return msg_.c_str();
  }

private:
  std::string msg_;
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
class ifile : public bio::filtering_istream {
private:
  string file;
  ifstream fd;

public:
  ifile();
  ifile(string filename , bool binary = false);
  ~ifile();
  string name();
  bool open(string filename, bool binary = false);
  bool readString(string &);
  void close();
};

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
class ofile : public bio::filtering_ostream {
private:
  string file;
  ofstream fd;

public:
  ofile();
  ofile(string filename , bool binary = false);
  ~ofile();
  string name();
  bool open(string filename, bool binary = false);
  void writeString(string &);
  void close();
};

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
class lfile {
private:
  string file;
  ofstream fd;
  bool verboseC;
  bool verboseL;

public:
  lfile();
  ~lfile();
  string name();
  bool open(string filename = "file.log");
  void close();
  string getPrefix();
  void muteL();
  void unmuteL();
  void muteC();
  void unmuteC();
  void print(string s);
  void printC(string s);
  void printL(string s);
  void println(string s);
  void printlnC(string s);
  void printlnL(string s);
  void warning(string s);
  void error(string s);
};

////////////////////
// Matteo's utils //
////////////////////

int count_lines(string filename);
void match_indices(const vector<string>& v1, const vector<string>& v2_, vector<int>& v3);
void right_join(const vector<string>& v1, const vector<string>& v2_, vector<int>& v3);
void arg_sort(const vector<string>& v, vector<int>& out);
set<int> set_intersection(const set<int> & s1, const set<int> & s2);
set<int> set_union(const set<int> & s1, const set<int> & s2);
pair<int,int> segment_intersection(int j1_min, int j1_max, int j2_min, int j2_max);

const boost::random::uniform_real_distribution<> runif;
int weighted_choice(const std::vector<double> & weights, boost::random::taus88 & rng);

void throw_error_bug(const string message);
void throw_error_input(const string message);

#endif
