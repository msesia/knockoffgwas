#ifndef UTILS_GRAPH_CPP
#define UTILS_GRAPH_CPP

#include "utils_graph.h"

using namespace boost;

void connected_components(const vector< vector<int> > & input, vector< vector<int> > & output) {
  typedef adjacency_list <vecS, vecS, undirectedS> Graph;
  typedef graph_traits<Graph>::vertex_descriptor Vertex;
  typedef graph_traits<Graph>::vertices_size_type VertexIndex;

  const int VERTEX_COUNT = input.size();
  Graph graph(VERTEX_COUNT);

  std::vector<VertexIndex> rank(num_vertices(graph));
  std::vector<Vertex> parent(num_vertices(graph));

  typedef VertexIndex* Rank;
  typedef Vertex* Parent;

  disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);

  initialize_incremental_components(graph, ds);
  incremental_components(graph, ds);

  graph_traits<Graph>::edge_descriptor edge;
  bool flag;

  for(int i=0; i<input.size(); i++) {
    if(input[i].size()>0) {
      for(auto it = input[i].begin(); it!=input[i].end(); ++it) {
        int j = *it;
        if(i<=j) {
          boost::tie(edge, flag) = add_edge(i, j, graph);
          ds.union_set(i,j);
        }
      }
    }
  }

  // std::cout << "An undirected graph:" << std::endl;
  // print_graph(graph, get(boost::vertex_index, graph));
  // std::cout << std::endl;

  typedef component_index<VertexIndex> Components;
  Components components(parent.begin(), parent.end());

  // Create output
  output.clear();
  output.resize(input.size());

  // Iterate through the component indices
  BOOST_FOREACH(VertexIndex current_index, components) {
    int component_size = 0;
    BOOST_FOREACH(VertexIndex child_index, components[current_index]) {
      component_size++;
    }

    // Iterate through the child vertex indices for [current_index]
    BOOST_FOREACH(VertexIndex i, components[current_index]) {
      BOOST_FOREACH(VertexIndex j, components[current_index]) {
        output[i].push_back(j);
      }
    }

  }

}

#endif
