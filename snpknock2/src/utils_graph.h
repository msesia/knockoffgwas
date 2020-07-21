#ifndef UTILS_GRAPH_H
#define UTILS_GRAPH_H

#include <vector>
#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/pending/disjoint_sets.hpp>

using namespace std;

void connected_components(const vector< vector<int> > & input, vector< vector<int> > & output);

#endif
