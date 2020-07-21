#ifndef _IBD_CPP
#define _IBD_CPP

#define DEBUG 0

#include "ibd.h"

IbdCluster::IbdCluster(const set<IbdSeg> & segments_) : segments(segments_) {
  if(DEBUG) cout << endl << endl << "[DEBUG] in IbdCluster::IbdCluster()" << endl;

  // Show cluster
  if(DEBUG) print();

  // Tidy
  if(DEBUG) cout << "[DEBUG] tidy()" << endl;
  tidy();

  // Show cluster
  if(DEBUG) print();

}

void IbdCluster::get(int k, IbdSeg & placeholder) const {
  auto it = std::next(segments.begin(), k);
  placeholder = *it;
}

int IbdCluster::size() const {
  return(segments.size());
}

void IbdCluster::print() const {
  // Print list of IBD segments within this cluster
  cout << "[DEBUG] This cluster contains " << segments.size() << " IBD segments: " << endl;
  for(auto ibd : segments) {
    cout << "  ";
    ibd.print();
    cout << endl;
  }
}

set<IbdSeg>::iterator IbdCluster::find_short_segments(int min_length) const {
  // Look for short segments (in BP)
  for(int i1=0; i1<segments.size(); i1++) {
    auto it = std::next(segments.begin(), i1);
    if((*it).bp_max - (*it).bp_min < min_length) {
      return(it);
    }
  }
  return(segments.end());
}

pair<set<IbdSeg>::iterator,set<IbdSeg>::iterator> IbdCluster::find_overlapping_segments() const {
  // Look for a pairwise overlaps between IBDs in this cluster
  pair<set<IbdSeg>::iterator, set<IbdSeg>::iterator > output;
  for(int i1=0; i1<segments.size(); i1++) {
    for(int i2=i1+1; i2<segments.size(); i2++) {
      auto it1 = std::next(segments.begin(), i1);
      auto it2 = std::next(segments.begin(), i2);
      if((*it1).overlaps_with(*it2)>0) {
        return(std::make_pair(it1, it2));
      }
    }
  }
  return(std::make_pair(segments.end(),segments.end()));
}

pair<set<IbdSeg>::iterator,set<IbdSeg>::iterator> IbdCluster::find_touching_segments() const {
  // Look for a pairwise overlaps between IBDs in this cluster
  pair<set<IbdSeg>::iterator, set<IbdSeg>::iterator > output;
  for(int i1=0; i1<segments.size(); i1++) {
    for(int i2=i1+1; i2<segments.size(); i2++) {
      auto it1 = std::next(segments.begin(), i1);
      auto it2 = std::next(segments.begin(), i2);
      if((*it1).touches_with(*it2)>0) {
        return(std::make_pair(it1, it2));
      }
    }
  }
  return(std::make_pair(segments.end(),segments.end()));
}

pair<set<IbdSeg>::iterator,set<IbdSeg>::iterator>
IbdCluster::find_close_segments(int min_dist, bool same_indices) const {
  // Look for pairs of segments (involving the same haplotypes) that are close to each other
  pair<set<IbdSeg>::iterator, set<IbdSeg>::iterator > output;
  for(int i1=0; i1<segments.size(); i1++) {
    for(int i2=i1+1; i2<segments.size(); i2++) {
      auto it1 = std::next(segments.begin(), i1);
      auto it2 = std::next(segments.begin(), i2);
      if((*it1).close_to(*it2, min_dist, same_indices)) {
        return(std::make_pair(it1, it2));
      }
    }
  }
  return(std::make_pair(segments.end(),segments.end()));
}

void IbdCluster::tidy() {
  const bool do_tidy_overlapping = true;
  const bool do_tidy_touching = true;
  const bool do_tidy_short = true;
  const bool do_merge_close = true;
  const bool do_expand_close = false;

  const int min_length = (int) 0.1e6;
  const int min_dist = (int) 0;

  bool found_overlap = do_tidy_overlapping;
  bool found_contact = do_tidy_touching;
  bool found_short = do_tidy_short;
  bool found_merge_close = do_merge_close;
  bool found_expand_close = do_expand_close;
  bool complete = false;

  while(!complete) {

    // Look for pairwise overlaps and tidy them until none is left
    if(do_tidy_overlapping) {
      auto it = find_overlapping_segments();
      if( (it.first != segments.end()) && (it.second != segments.end()) ) {
        found_overlap = true;
        tidy_overlapping (*(it.first), *(it.second), min_length);
        // print(); // DEBUG
      } else {
        found_overlap = false;
      }
    }

    // Look for segments touching each other and separate them until none is left
    if(do_tidy_touching) {
      auto it = find_touching_segments();
      if( (it.first != segments.end()) && (it.second != segments.end()) ) {
        found_contact = true;
        tidy_touching(*(it.first), *(it.second), min_length);
        // print(); // DEBUG
      } else {
        found_contact = false;
      }
    }

    // Remove short segments
    if(do_tidy_short) {
      auto it = find_short_segments(min_length);
      if(it != segments.end()) {
        found_short = true;
        if(DEBUG) {
          cout << "Removing segment ";
          (*it).print();
          cout << endl;
        }
        segments.erase(it);
        //print(); // DEBUG
      } else {
        found_short = false;
      }
    }

    // Look for nearby segments involving the same haplotypes and merge them until none is left
    if(do_merge_close) {
      auto it = find_close_segments(min_dist, true);
      if( (it.first != segments.end()) && (it.second != segments.end()) ) {
        found_merge_close = true;
        merge_close(*(it.first), *(it.second));
        //print(); // DEBUG
      } else {
        found_merge_close = false;
      }
    }

    // Look for nearby segments involving different haplotypes and expand one of them to remove gap
    if(do_expand_close) {
      auto it = find_close_segments(min_dist, false);
      if( (it.first != segments.end()) && (it.second != segments.end()) ) {
        found_expand_close = true;
        expand_close(*(it.first), *(it.second));
        //print(); // DEBUG
      } else {
        found_expand_close = false;
      }
    }

    if( (!found_overlap) && (!found_contact) && (!found_short) &&
        (!found_merge_close) && (!found_expand_close) ) {
      complete = true;
    }

  }
}

void IbdCluster::tidy_touching(const IbdSeg & s1, const IbdSeg & s2, int min_length) {
  //assert(s1.size()>1);
  //assert(s2.size()>1);
  // Find longer segment and shorten it
  IbdSeg s1_new = s1;
  IbdSeg s2_new = s2;
  if(s1.j_max==s2.j_min) {
    if(s1.size() >= s2.size()) {
      s1_new = s1.shave_from_end();
      s2_new = s2;
    } else {
      s1_new = s1;
      s2_new = s2.shave_from_start();
    }
  } else {
    if(s1.size() >= s2.size()) {
      s1_new = s1.shave_from_start();
      s2_new = s2;
    } else {
      s1_new = s1;
      s2_new = s2.shave_from_end();
    }
  }

  // Collect all the pieces
  vector<IbdSeg> new_segments = {s1_new, s2_new};

  // Debug
  if(DEBUG) {
    cout << "Separating segments ";
    s1.print();
    s2.print();
    cout << " into:" << endl;
    for(auto s : new_segments) {
      cout << "  ";
      s.print();
      cout << endl;
    }
  }

  // Erase old segments
  segments.erase(s1);
  segments.erase(s2);

  // Insert new segments
  for(auto s : new_segments) {
    segments.insert(s);
  }
}

void IbdCluster::tidy_overlapping(const IbdSeg & s1, const IbdSeg & s2, int min_length) {
  // Find shared segment
  IbdSeg shared = s1.find_shared(s2);
  // See what's left when the shared segment is removed
  vector<IbdSeg> remaining1 = s1.subtract(shared);
  vector<IbdSeg> remaining2 = s2.subtract(shared);
  // Collect all the pieces
  vector<IbdSeg> new_segments;
  if(shared.length()>min_length) new_segments.push_back(shared);
  for(auto s : remaining1) {
    if(s.length()>min_length) new_segments.push_back(s);
  }
  for(auto s : remaining2) {
    if(s.length()>min_length) new_segments.push_back(s);
  }

  // Debug
  if(DEBUG) {
    cout << "Splitting segments ";
    s1.print();
    s2.print();
    cout << " into:" << endl;
    for(auto s : new_segments) {
      cout << "  ";
      s.print();
      cout << endl;
    }
  }

  // Erase old segments
  segments.erase(s1);
  segments.erase(s2);

  // Insert new segments
  for(auto s : new_segments) {
    segments.insert(s);
  }
}

void IbdCluster::expand_close(const IbdSeg & s1, const IbdSeg & s2) {
  int j1_min, j1_max, bp1_min, bp1_max;
  int j2_min, j2_max, bp2_min, bp2_max;

  if(s1.j_min <= s2.j_min) {
    // Update first segment
    j1_min = s1.j_min;
    bp1_min = s1.bp_min;
    j1_max = (s1.j_max+s2.j_min)/2;
    bp1_max = (s1.bp_max+s2.bp_min)/2;
    // Update second segment
    j2_min = (s1.j_max+s2.j_min)/2;
    bp2_min = (s1.bp_max+s2.bp_min)/2;
    j2_max = s2.j_max;
    bp2_max = s2.bp_max;
  } else {
    expand_close(s2, s1);
    return;
  }
  IbdSeg new_s1(s1.indices, j1_min, j1_max, bp1_min, bp1_max);
  IbdSeg new_s2(s2.indices, j2_min, j2_max, bp2_min, bp2_max);

  // Debug
  if(DEBUG) {
    cout << "Expanding segments ";
    s1.print();
    s2.print();
    cout << " into:" << endl;
    cout << "  ";
    new_s1.print();
    new_s2.print();
    cout << endl;
  }

  // Erase old segments
  segments.erase(s1);
  segments.erase(s2);

  // Insert new segments
  segments.insert(new_s1);
  segments.insert(new_s2);
}

void IbdCluster::merge_close(const IbdSeg & s1, const IbdSeg & s2) {
  // Merge indices
  set<int> new_indices = set_union(s1.indices, s2.indices);
  // Merge variant ranges, inserting gap
  int j_min = std::min(s1.j_min, s2.j_min);
  int bp_min = std::min(s1.bp_min, s2.bp_min);
  int j_max = std::max(s1.j_max, s2.j_max);
  int bp_max = std::max(s1.bp_max, s2.bp_max);
  IbdSeg new_segment(new_indices, j_min, j_max, bp_min, bp_max);

  // Debug
  if(DEBUG) {
    cout << "Merging segments ";
    s1.print();
    s2.print();
    cout << " into:" << endl;
    cout << "  ";
    new_segment.print();
    cout << endl;
  }

  // Erase old segments
  segments.erase(s1);
  segments.erase(s2);

  // Insert new segments
  segments.insert(new_segment);
}

IbdSeg::IbdSeg() {
}

IbdSeg::IbdSeg(const set<int> & indices_, int j_min_, int j_max_, int bp_min_, int bp_max_) {
  indices = indices_;
  j_min = j_min_;
  j_max = j_max_;
  bp_min = bp_min_;
  bp_max = bp_max_;
}

IbdSeg::IbdSeg(const IbdSeg& obj) {
  indices = obj.indices;
  j_min = obj.j_min;
  j_max = obj.j_max;
  bp_min = obj.bp_min;
  bp_max = obj.bp_max;
}

void IbdSeg::print() const {
  cout << "( ";
  for(auto i : indices) cout << i << " ";
  cout << "; " << j_min << "--" << j_max << " ; ";
  cout << std::setprecision(4) << (double)bp_min/1e6 << "--";
  cout << std::setprecision(4) << (double)bp_max/1e6 << " Mb ) ";
}

IbdSeg IbdSeg::find_shared(const IbdSeg & obj) const {
  set<int> indices_union = set_union(indices, obj.indices);
  pair<int,int> j_range = segment_intersection(j_min, j_max, obj.j_min, obj.j_max);
  pair<int,int> bp_range = segment_intersection(bp_min, bp_max, obj.bp_min, obj.bp_max);
  IbdSeg intersection(indices_union, j_range.first, j_range.second, bp_range.first, bp_range.second);
  return(intersection);
}

vector<IbdSeg> IbdSeg::subtract(const IbdSeg & obj) const {
  vector<IbdSeg> output;

  // TODO: make sure this is correct in general
  if(j_min == obj.j_min) {
    // If the two segments share the left-hand side
    output.emplace_back(indices, obj.j_max, j_max, obj.bp_max, bp_max);
  } else if(j_max == obj.j_max) {
    // If the two segments share the right-hand side
    output.emplace_back(indices, j_min, obj.j_min, bp_min, obj.bp_min);
  } else {
    // If the second segment is in the middle of the first segment
    output.emplace_back(indices, j_min, obj.j_min, bp_min, obj.bp_min);
    output.emplace_back(indices, obj.j_max, j_max, obj.bp_max, bp_max);
  }
  return(output);
}

bool IbdSeg::close_to(const IbdSeg & obj, int min_dist, bool same_indices) const {
  if(same_indices) {
    // Check whether two segments share all haplotype indices
    if(indices.size() != obj.indices.size()) return(false);
    const int intersect_size = set_intersection(indices, obj.indices).size();
    if((intersect_size!=indices.size())|(intersect_size!=obj.indices.size())) return(false);
    // Check whether two segments are physically close to each other
    if(bp_min <= obj.bp_min) {
      if(bp_max >= obj.bp_min - min_dist) return(true);
    } else {
      if(obj.bp_max >= bp_min - min_dist) return(true);
    }
  } else {
    // Check whether two segments are physically close to each other, but not touching
    if(bp_min <= obj.bp_min) {
      if( (bp_max >= obj.bp_min - min_dist) && (bp_max < obj.bp_min) ) return(true);
    } else {
      if( (obj.bp_max >= bp_min - min_dist) && (obj.bp_max < bp_min) ) return(true);
    }
  }
  return(false);
}

int IbdSeg::length() const {
  return(bp_max-bp_min+1);
}

int IbdSeg::size() const {
  return(j_max-j_min+1);
}

int IbdSeg::overlaps_with(const IbdSeg & obj) const {
  // Check whether two segments share at least one haplotype index and overlap physically
  if(set_intersection(indices, obj.indices).size()==0) {
    // If the haplotype indices are all different, there can be no intersection
    return(-1);
  } else {
    // If at least one of the haplotype indices matches, check the variant ranges
    if(j_min <= obj.j_min) {
      if(j_max < obj.j_min) {
        return(-1);
      } else {
        return(j_max - obj.j_min);
      }
    } else {
      if(obj.j_max < j_min) {
        return(-1);
      } else {
        return(obj.j_max - j_min);
      }
    }
  }
}

int IbdSeg::touches_with(const IbdSeg & obj) const {
  // Check whether two segments share at least one haplotype index and overlap physically
  if(set_intersection(indices, obj.indices).size()==0) {
    // If the haplotype indices are all different, there can be no intersection
    return(-1);
  } else {
    // If at least one of the haplotype indices matches, check the variant ranges
    if(j_min <= obj.j_min) {
      if(j_max == obj.j_min) {
        return(1);
      } else {
        return(0);
      }
    } else {
      if(obj.j_max == j_min) {
        return(1);
      } else {
        return(0);
      }
    }
  }
}

bool IbdSeg::operator< (const IbdSeg & obj) const {
  if((*this)==obj) return(false);
  if(indices.size()<obj.indices.size()) {
    return(true);
  }
  else if(indices.size()>obj.indices.size()) {
    return(false);
  }
  else {
    if(indices==obj.indices) {
      // Compare the variant indices
      if(j_min < obj.j_min) {
        return(true);
      } else if (j_min > obj.j_min) {
        return(false);
      } else {
        if(j_max < obj.j_max) {
          return(true);
        } else if (j_max > obj.j_max) {
          return(false);
        }
        return(false);
      }
    } else {
      return(indices<obj.indices);
    }
  }
}

bool IbdSeg::operator== (const IbdSeg & obj) const {
  if(indices.size()!=obj.indices.size()) return(false);
  if(indices!=obj.indices) return(false);
  if((j_min==obj.j_min)&&(j_max==obj.j_max)) return(true);
  else return(false);
}

IbdSeg IbdSeg::shave_from_start() const {
  IbdSeg new_obj(indices, j_min+1, j_max, bp_min, bp_max);
  return(new_obj);
}

IbdSeg IbdSeg::shave_from_end() const {
  IbdSeg new_obj(indices, j_min, j_max-1, bp_min, bp_max);
  return(new_obj);
}

#endif
