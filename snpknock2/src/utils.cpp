#ifndef _UTILS_CPP
#define _UTILS_CPP

#include "utils.h"

/******************************************************/
/*                  UTILS STATISTICS                  */
/******************************************************/
long seed = -123456789;

namespace putils {
  void initRandom(long s) {
    seed = - s;
  }

  double getRandom() {
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    if (seed <= 0 || !iy) {
      if (-(seed) < 1) seed=1;
      else seed = -(seed);
      for (j=NTAB+7;j>=0;j--) {
        k=(seed)/IQ;
        seed=IA*(seed-k*IQ)-IR*k;
        if (seed < 0) seed += IM;
        if (j < NTAB) iv[j] = seed;
      }
      iy=iv[0];
    }
    k=(seed)/IQ;
    seed=IA*(seed-k*IQ)-IR*k;
    if (seed < 0) seed += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = seed;
    temp=AM*iy;
    if (temp > RNMX) return RNMX;
    else return temp;
  }

  string getRandomID(){
    return boost::lexical_cast<string>(boost::uuids::random_generator()());
  }

  long getSeed() {
    return seed;
  }

  int getRandom(int n) {
    return (int)floor(getRandom() * n);
  }

  void normalise(vector < double > & v) {
    double sum = 0.0;
    for (int i = 0 ; i < v.size() ; i++) sum += v[i];
    if (sum != 0.0) for (int i = 0 ; i < v.size() ; i++) v[i] /= sum;
  }

  int sample(vector< double > & v, double sum) {
    double csum = v[0];
    double u = getRandom() * sum;
    for (int i = 0; i < v.size() - 1; ++i) {
      if ( u < csum ) return i;
      csum += v[i+1];
    }
    return v.size() - 1;
  }

  double entropy(vector < double > & v) {
    double e = 0.0;
    for (int i = 0 ; i < v.size() ; i ++) {
      if (v[i] > 0.0) e += v[i] * log(v[i]);
    }
    return e;
  }

  double KLdistance(vector < double > & P, vector < double > & Q) {
    assert(P.size() == Q.size());
    double d = 0.0;
    for (int i = 0 ; i < Q.size() ; i ++) {
      if (Q[i] > 0.0 && P[i] > 0.0) d += P[i] * (log(P[i]) - log(Q[i]));
    }
    return d;
  }
};

/******************************************************/
/*                  UTILS ALGORITHM                   */
/******************************************************/
namespace autils {
  int max(vector < double > & v) {
    double max = -1e300;
    int index_max = 0;
    for (int i = 0 ; i < v.size() ; i ++)
      if (v[i] > max) {
        max = v[i];
        index_max = i;
      }
    return index_max;
  }

  int max(vector < int > & v) {
    int max = -1000000000;
    int index_max = 0;
    for (int i = 0 ; i < v.size() ; i ++)
      if (v[i] > max) {
        max = v[i];
        index_max = i;
      }
    return index_max;
  }

  void findUniqueSet(vector < bool > & B, vector < int > & U) {
    U.clear();
    for (int b = 0 ; b < B.size() ; b ++)
      if (B[b]) U.push_back(b);
    //sort( U.begin(), U.end() );
    //U.erase(unique( U.begin(), U.end() ), U.end());
  }

  void decompose(int min, vector < vector < int > > & B, vector < vector < vector < int > > > & BB) {
    if (B.size() < 2 * min || B.size() == 2) {
      BB.push_back(B);
      return;
    }

    int l = putils::getRandom(B.size() - 2 * min) + min;
    vector < vector < int > > LB = vector < vector < int > > (B.begin(), B.begin() + l + 1);
    vector < vector < int > > RB = vector < vector < int > > (B.begin() + l - 1, B.end());
    vector < vector < vector < int > > > LBB;
    vector < vector < vector < int > > > RBB;
    decompose(min, LB, LBB);
    decompose(min, RB, RBB);
    BB = LBB;
    BB.insert(BB.end(), RBB.begin(), RBB.end());
  }

  int checkDuo (int pa1, int pa2, int ca1, int ca2) {
    if (pa1 == 0 && pa2 == 0 && ca1 == 0 && ca2 == 0) { return 0; }
    if (pa1 == 0 && pa2 == 0 && ca1 == 0 && ca2 == 1) { return 0; }
    if (pa1 == 0 && pa2 == 0 && ca1 == 1 && ca2 == 0) { return 0; }
    if (pa1 == 0 && pa2 == 0 && ca1 == 1 && ca2 == 1) { return -1; }
    if (pa1 == 0 && pa2 == 1 && ca1 == 0 && ca2 == 0) { return 1; }
    if (pa1 == 0 && pa2 == 1 && ca1 == 0 && ca2 == 1) { return 0; }
    if (pa1 == 0 && pa2 == 1 && ca1 == 1 && ca2 == 0) { return 0; }
    if (pa1 == 0 && pa2 == 1 && ca1 == 1 && ca2 == 1) { return 1; }
    if (pa1 == 1 && pa2 == 0 && ca1 == 0 && ca2 == 0) { return 1; }
    if (pa1 == 1 && pa2 == 0 && ca1 == 0 && ca2 == 1) { return 0; }
    if (pa1 == 1 && pa2 == 0 && ca1 == 1 && ca2 == 0) { return 0; }
    if (pa1 == 1 && pa2 == 0 && ca1 == 1 && ca2 == 1) { return 1; }
    if (pa1 == 1 && pa2 == 1 && ca1 == 0 && ca2 == 0) { return -1; }
    if (pa1 == 1 && pa2 == 1 && ca1 == 0 && ca2 == 1) { return 0; }
    if (pa1 == 1 && pa2 == 1 && ca1 == 1 && ca2 == 0) { return 0; }
    if (pa1 == 1 && pa2 == 1 && ca1 == 1 && ca2 == 1) { return 0; }
    return 0;
  }

  int checkTrio (int fa1, int fa2, int ma1, int ma2, int ca1, int ca2) {
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 0; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return -1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return -1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return -1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
    if (fa1 == 0 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return -1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 2; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 2; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 2; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 2; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 0 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return 2; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 0; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 0; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 2; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return 2; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 2; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 1 && fa2 == 0 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return -1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 0; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 0; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return -1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 0 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 0 ) { return -1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 0 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 0 ) { return 1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 0 && ca1 == 1 && ca2 == 1 ) { return 1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 0 ) { return -1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 0 && ca2 == 1 ) { return -1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 0 ) { return -1; }
    if (fa1 == 1 && fa2 == 1 && ma1 == 1 && ma2 == 1 && ca1 == 1 && ca2 == 1 ) { return 0; }
    return 0;
  }
};

/******************************************************/
/*                  UTILS STRING                      */
/******************************************************/
namespace sutils {
  int tokenize (string & str, vector < string > & tokens) {
    tokens.clear();
    string::size_type p_last = str.find_first_not_of(" 	", 0);
    string::size_type p_curr = str.find_first_of(" 	", p_last);
    while (string::npos != p_curr || string::npos != p_last) {
      tokens.push_back(str.substr(p_last, p_curr - p_last));
      p_last = str.find_first_not_of("	 ", p_curr);
      p_curr = str.find_first_of("	 ", p_last);
    }
    if (tokens.back()[tokens.back().size()-1] == '\r')
      tokens.back() = tokens.back().substr(0, tokens.back().size()-1);
    return tokens.size();
  }

  int tokenize(string & str, vector < string > & tokens, int n_max_tokens) {
    tokens.clear();
    string::size_type p_last = str.find_first_not_of(" 	", 0);
    string::size_type p_curr = str.find_first_of(" 	", p_last);
    while ((string::npos != p_curr || string::npos != p_last) && tokens.size() < n_max_tokens) {
      tokens.push_back(str.substr(p_last, p_curr - p_last));
      p_last = str.find_first_not_of("	 ", p_curr);
      p_curr = str.find_first_of("	 ", p_last);
    }
    return tokens.size();
  }

  string int2str(int n) {
    ostringstream s2( stringstream::out );
    s2 << n;
    return s2.str();
  }

  string int2str(vector < int > & v) {
    ostringstream s2( stringstream::out );
    for (int l = 0 ; l < v.size() ; l++) {
      if (v[l] < 0) s2 << "-";
      else s2 << v[l] ;
    }
    return s2.str();
  }

  string long2str(long int n) {
    ostringstream s2( stringstream::out );
    s2 << n;
    return s2.str();
  }

  string double2str(double n, int prc) {
    ostringstream s2;
    s2 << setiosflags( ios::fixed );
    if ( prc > 0 ) s2.precision(prc);
    s2 << n;
    return s2.str();
  }

  string double2str(vector < double > &v, int prc) {
    ostringstream s2;
    s2 << setiosflags( ios::fixed );
    if ( prc >= 0 ) s2.precision(prc);
    for (int l = 0 ; l < v.size() ; l++) {
      s2 << v[l] << " ";
    }
    return s2.str();
  }

  string bool2str(vector<bool> & v) {
    ostringstream s2( stringstream::out );
    for (int l = 0 ; l < v.size() ; l++) {
      if (v[l]) s2 << "1";
      else s2 << "0";
    }
    return s2.str();
  }

  string date2str(time_t * t, string format) {
    struct tm * timeinfo = localtime(t);
    char buffer[128];
    strftime(buffer, 128, format.c_str(), timeinfo);
    ostringstream s2( stringstream::out );
    s2 << buffer;
    return s2.str();
  }
};

/******************************************************/
/*                  UTILS FILE                        */
/******************************************************/
namespace futils {
  bool isFile(string f) {
    ifstream inp;
    inp.open(f.c_str(), ifstream::in);
    if(inp.fail()) {
      inp.clear(ios::failbit);
      inp.close();
      return false;
    }
    inp.close();
    return true;
  }

  bool createFile(string f) {
    ofstream out;
    out.open(f.c_str(), ofstream::out);
    if(out.fail()) {
      out.clear(ios::failbit);
      out.close();
      return false;
    }
    out.close();
    return true;
  }

  string extensionFile(string & filename) {
    if (filename.find_last_of(".") != string::npos)
      return filename.substr(filename.find_last_of(".") + 1);
    return "";
  }


  void bool2binary(vector < bool > & V, ostream &fd) {
    int nb = V.size();
    fd.write((char*)&nb, sizeof(int));
    int cpt_byte = 0;
    int cpt_bin = 0;
    int nb_byte = (int)ceil( (V.size() * 1.0) / 8);
    while (cpt_byte < nb_byte) {
      bitset<8> byte_bitset;
      for (int l = 0; l < 8 && cpt_bin < V.size() ; l++) {
        byte_bitset[l] = V[cpt_bin];
        cpt_bin ++;
      }
      char byte_char = (char)byte_bitset.to_ulong();
      fd.write(&byte_char, 1);
      cpt_byte++;
    }
  }

  bool binary2bool(vector < bool > & V, istream & fd) {
    int nb;
    fd.read((char*)&nb, sizeof(int));
    if (!fd) return false;
    int cpt_byte = 0;
    int cpt_bin = 0;
    int nb_byte = (int)ceil( (nb * 1.0) / 8);
    V = vector < bool >(nb);
    while (cpt_byte < nb_byte) {
      char byte_char;
      fd.read(&byte_char, 1);
      if (!fd) return false;
      bitset<8> byte_bitset = byte_char;
      for (int l = 0; l < 8 && cpt_bin < nb ; l++) {
        V[cpt_bin] = byte_bitset[l];
        cpt_bin++;
      }
      cpt_byte++;
    }
    return true;
  }
};

/******************************************************/
/*                  INPUT FILE                        */
/******************************************************/
ifile::ifile() {
}

ifile::ifile(string filename , bool binary) {
  open(filename, binary);
}

ifile::~ifile() {
  close();
}

string ifile::name() {
  return file;
}

bool ifile::open(string filename, bool binary) {
  file = filename;
  string ext = futils::extensionFile(filename);
  if (ext == "gz") {
    fd.open(file.c_str(), ios::in | ios::binary);
    push(bio::gzip_decompressor());
  } else if (ext == "bz2") {
    fd.open(file.c_str(), ios::in | ios::binary);
    push(bio::bzip2_decompressor());
  } else if (binary) {
    fd.open(file.c_str(), ios::in | ios::binary);
  } else  {
    fd.open(file.c_str());
  }
  if (fd.fail()) return false;
  push(fd);
  return true;
}

bool ifile::readString(string & str) {
  int s;
  if (!read((char*)&s, sizeof(int))) return false;
  char  * buffer = new char [s + 1];
  if (!read(buffer, s)) return false;
  buffer[s] = '\0';
  str = string(buffer);
  delete [] buffer;
  return true;
}

void ifile::close() {
  if (!empty()) reset();
  fd.close();
}

/******************************************************/
/*                  OUTPUT FILE                       */
/******************************************************/
ofile::ofile() {
}

ofile::ofile(string filename , bool binary) {
  open(filename, binary);
}

ofile::~ofile() {
  close();
}

string ofile::name() {
  return file;
}

bool ofile::open(string filename, bool binary) {
  file = filename;
  string ext = futils::extensionFile(filename);
  if (ext == "gz") {
    fd.open(file.c_str(), ios::out | ios::binary);
    push(bio::gzip_compressor());
  } else if (ext == "bz2") {
    fd.open(file.c_str(), ios::out | ios::binary);
    push(bio::bzip2_compressor());
  } else if (binary) {
    fd.open(file.c_str(), ios::out | ios::binary);
  } else  {
    fd.open(file.c_str());
  }
  if (fd.fail()) return false;
  push(fd);
  return true;
}

void ofile::writeString(string & str) {
  int s = str.size();
  write((char*)&s, sizeof(int));
  write(str.c_str(), s * sizeof(char));
}

void ofile::close() {
  if (!empty()) reset();
  fd.close();
}

/******************************************************/
/*                  LOG FILE                          */
/******************************************************/
lfile::lfile() {
  verboseC = true;
  verboseL = true;
}

lfile::~lfile() {
  close();
}

string lfile::name() {
  return file;
}

bool lfile::open(string filename) {
  file = filename;
  if (futils::extensionFile(file) != "log") file += ".log";
  fd.open(file.c_str());
  if (fd.fail()) return false;
  return true;
}

void lfile::close() {
  fd.close();
}

string lfile::getPrefix() {
  return file.substr(0, file.find_last_of("."));
}

void lfile::muteL() {
  verboseL = false;
}

void lfile::unmuteL() {
  verboseL = true;
}

void lfile::muteC() {
  verboseC = false;
}

void lfile::unmuteC() {
  verboseC = true;
}


void lfile::print(string s) {
  if (verboseL) { fd << s; fd.flush(); }
  if (verboseC) { cout << s; cout.flush(); }
}

void lfile::printC(string s) {
  if (verboseC) { cout << s; cout.flush(); }
}

void lfile::printL(string s) {
  if (verboseL) { fd << s; fd.flush(); }
}

void lfile::println(string s) {
  if (verboseL) { fd << s << endl; }
  if (verboseC) { cout << s << endl; }
}

void lfile::printlnC(string s) {
  if (verboseC) { cout << s << endl; }
}

void lfile::printlnL(string s) {
  if (verboseL) { fd << s << endl; }
}

void lfile::warning(string s) {
  cout << "\33[33mWARNING:\33[0m " << s << endl;
  if (verboseL) fd << "WARNING: " << s << endl;
}

void lfile::error(string s) {
  cout << "\33[33mERROR:\33[0m " << s << endl;
  if (verboseL) fd << "ERROR: " << s << endl;
  close();
  exit(1);
}

////////////////////
// Matteo's utils //
////////////////////

int count_lines(string filename) {
  int num_lines = 0;
  string line;
  ifstream file(filename);
  while (std::getline(file, line)) num_lines++;
  file.close();
  return(num_lines);
}

void match_indices(const vector<string>& v1, const vector<string>& v2_, vector<int>& v3) {
  // Returns indices of v1 elements that match with at least one element in v2
  vector<string> v2(v2_);
  std::sort(v2.begin(), v2.end());
  v3.clear();
  v3.reserve(v2.size());
  for (int i=0; i<v1.size(); i++) {
    if (std::binary_search(v2.begin(), v2.end(), v1[i])) v3.push_back(i);
  }
}

void right_join(const vector<string>& v1, const vector<string>& v2_, vector<int>& v3) {
  // For each value in v1, return the index of the matching v2 element
  vector<string> v2(v2_);
  std::sort(v2.begin(), v2.end());
  vector<int> indices(v2.size());
  arg_sort(v2_, indices);
  v3.clear();
  v3.resize(v1.size());
  for (int i=0; i<v1.size(); i++) {
    auto lower = std::lower_bound(v2.begin(), v2.end(), v1[i]);
    // check that value has been found
    const bool found = lower != v2.end() && *lower == v1[i];
    if(found) {
      int j = (int) std::distance(v2.begin(), lower);
      v3[i] = indices[j];
    } else {
      v3[i] = -1;
    }
  }
  // // DEBUG
  // cout << "[DEBUG]" << endl;
  // for (int i=0; i<10; i++) {
  //   cout << "v1[" << i << "]: " << v1[i] << "\t";
  //   cout << "v1[" << v3[i] << "]: " << v2_[v3[i]] << endl;
  // }
}

void arg_sort(const vector<string>& v, vector<int>& out) {
    // Vector to store element with respective present index
    vector< pair<string,int> > vp;
    // Insert element in pair vector to keep track of previous indexes
    for (int i = 0; i < v.size(); i++) {
      vp.push_back(make_pair(v[i], i));
    }
    // Sort pair vector
    sort(vp.begin(), vp.end());
    // Return indices
    out.clear();
    out.resize(v.size());
    for (int i = 0; i < v.size(); i++) {
      out[i] = vp[i].second;
    }
}

set<int> set_intersection(const set<int> & s1, const set<int> & s2) {
  vector<int> v(std::min(s1.size(),s2.size()));
  auto it = std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());
  v.resize(it-v.begin());
  set<int> output;
  for(auto i : v) output.insert(i);
  return(output);
}

set<int> set_union(const set<int> & s1, const set<int> & s2) {
  vector<int> v(s1.size()+s2.size());
  auto it = std::set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());
  v.resize(it-v.begin());
  set<int> output;
  for(auto i : v) output.insert(i);
  return(output);
}

pair<int,int> segment_intersection(int j1_min, int j1_max, int j2_min, int j2_max) {
  // Make sure the first segment has the left-most start
  if(j1_min > j2_min) {
    int j1_min_tmp = j1_min;
    j1_min = j2_min;
    j2_min = j1_min_tmp;
    int j1_max_tmp = j1_max;
    j1_max = j2_max;
    j2_max = j1_max_tmp;
  }
  int j_min, j_max;
  // Check whether the two segments overlap
  if(j1_max>=j2_min) {
    j_min = j2_min;
    j_max = std::min(j1_max, j2_max);
  } else {
    j_min = -1;
    j_max = -1;
  }
  return(std::make_pair(j_min, j_max));
}

int weighted_choice(const std::vector<double> & weights, boost::random::taus88 & rng) {
  double weights_sum = accumulate(weights.begin(),weights.end(),0.0);
  double R = runif(rng) * weights_sum;
  for(unsigned int i=0; i<weights.size(); i++) {
    if(R < weights[i]) return(i);
    R -= weights[i];
  }
  return(weights.size()-1);
}

void throw_error_bug(const string message) {
  cerr << "Error: " << message << endl;
  cerr << "This is likely to be due to a bug. Please report it." << endl;
  exit(-1);
}

void throw_error_input(const string message) {
  cerr << "Error: " << message << endl;
  cerr << "This is likely to be due to incorrectly formatted input." << endl;
  exit(-1);
}

#endif
