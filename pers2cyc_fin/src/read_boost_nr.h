//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

/*
  Reads maximal flow problem in extended DIMACS format.
  This works, but could use some polishing.
*/

/* ----------------------------------------------------------------- */

//#ifndef BOOST_GRAPH_READ_DIMACS_HPP
//#define BOOST_GRAPH_READ_DIMACS_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string.h>
#include <queue>
#include <limits.h>
#include <string.h>
#include <algorithm>
#include <boost/graph/graph_traits.hpp>
using namespace std;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;

map <int,int> threesimplex_to_node; // maps index of tetrahedrons from simplex->first to dual-graph-nodes->second
map <int,int> twosimplex_to_tri; // maps index of traingles from simplex->first to dual-graph-edges->second
map <int,int> threesimplex_to_node_rev; // maps index of tetrahedrons from dual-graph-nodes->first to simplex->second
map <int,int> twosimplex_to_tri_rev; // maps index of traingles from dual-graph-edges->first to simplex->second
map<size_t, vector<int>>  sizescale;
namespace boost {

  namespace detail {

template <class Graph, class CapacityMap, class ReverseEdgeMap>
int read_dimacs_max_flow_internal(Graph& g,
                                  CapacityMap capacity,
                                  ReverseEdgeMap reverse_edge,
                                  typename graph_traits<Graph>::vertex_descriptor& src,
                                  typename graph_traits<Graph>::vertex_descriptor& sink,
                                  bool require_source_and_sink,
                                  const std::string& problem_type,
                                  Bitmap_cubical_complex &cub,
                                  size_t birthi, size_t deathi)
{
 
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    double birth = cub.filtration(birthi);
    double death = cub.filtration(deathi);
  std::vector<vertex_descriptor> verts;

    long m=sizescale[deathi][0]*2, n= sizescale[deathi][1]+1,                    /*  number of edges and nodes */
    i, head, tail, cap;
    bool in1, in2;
    edge_descriptor e1, e2;
    
    for (long vi = 0; vi < n; ++vi)
          verts.push_back(add_vertex(g));
  
    src = threesimplex_to_node[deathi];
    sink = sizescale[deathi][1]-1;
    int ccc =6170175;

    for (size_t ind=0;ind<cub.num_simplices(); ind++){
        edge_descriptor e1, e2;
        bool in1, in2;
        if(cub.dimension(ind)==2 && cub.filtration(ind)<=birth){//Triangle:Push purple or green edges
            auto cof = cub.coboundary_range(ind);
            
            if(cof.size()==2 && cub.filtration(cof[0])<=death  && cub.filtration(cof[1])<=death){//G-G, G-P, P-P
                
                tail = threesimplex_to_node[cof[0]];
                head = threesimplex_to_node[cof[1]];
                cap = 1;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                
                tail = threesimplex_to_node[cof[1]];
                head = threesimplex_to_node[cof[0]];
                cap = 1;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
            }
            else if(cof.size()==2 && cub.filtration(cof[0])<=death){ // G-R, P-R
                
                tail = threesimplex_to_node[cof[0]];
                head = sizescale[deathi][1]-1;
                cap = 1;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                
                tail = sizescale[deathi][1]-1;
                head = threesimplex_to_node[cof[0]];
                cap = 1;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
            }
            else if(cof.size()==2 && cub.filtration(cof[1])<=death){ // G-R, P-R
                
                tail = threesimplex_to_node[cof[1]];
                head = sizescale[deathi][1]-1;
                cap = 1;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                
                tail = sizescale[deathi][1]-1;
                head = threesimplex_to_node[cof[1]];
                cap = 1;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
            }
            else if(cof.size()==1 && cub.filtration(cof[0])<=death){
                tail = threesimplex_to_node[cof[0]];
                head = sizescale[deathi][1]-1;
                cap = 1;

                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                
                tail = sizescale[deathi][1]-1;
                head = threesimplex_to_node[cof[0]];
                cap = 1;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
            }// Else both red, no need
        }
        else if(cub.dimension(ind)==2 && cub.filtration(ind)<=death){//Triangle:Push purple or green edges
            auto cof = cub.coboundary_range(ind);
            if(cof.size()==2 && cub.filtration(cof[0])<=death  && cub.filtration(cof[1])<=death){//G-G, G-P, P-P

                tail = threesimplex_to_node[cof[0]];
                head = threesimplex_to_node[cof[1]];
                if(deathi==2377095 && ind==ccc){cout<<"1: "<<verts.at(tail)<<" "<<verts.at(head);getchar();}
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                tail = threesimplex_to_node[cof[1]];
                head = threesimplex_to_node[cof[0]];
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                

            }
            else if(cof.size()==2 && cub.filtration(cof[0])<=death){ // G-R, P-R
                
                tail = threesimplex_to_node[cof[0]];
                head = sizescale[deathi][1]-1;
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                
                tail = sizescale[deathi][1]-1;
                head = threesimplex_to_node[cof[0]];
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
            }
            else if(cof.size()==2 && cub.filtration(cof[1])<=death){ // G-R, P-R
                
                tail = threesimplex_to_node[cof[1]];
                head = sizescale[deathi][1]-1;
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                
                tail = sizescale[deathi][1]-1;
                head = threesimplex_to_node[cof[1]];
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
            }
            else if(cof.size()==1 && cub.filtration(cof[0])<=death){
                tail = threesimplex_to_node[cof[0]];
                head = sizescale[deathi][1]-1;
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
                
                tail = sizescale[deathi][1]-1;
                head = threesimplex_to_node[cof[0]];
                cap = INT_MAX;
                boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
                boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
                capacity[e1] = cap;
                capacity[e2] = 0;
                reverse_edge[e1] = e2;
                reverse_edge[e2] = e1;
//                cout<<"a "<<tail<<" "<<sizescale[deathi][1]-1<<" "<<INT_MAX<<"\n";
//                ft<<"a "<<sizescale[deathi][1]-1<<" "<<graphInd1<<" "<<INT_MAX<<"\n";
            }// Else both red, no need
        }
//        if(deathi==2377095 && ind==ccc)
//        {cout<<"out of else if";getchar();}
    }
    
//          std::sscanf ( in_line.c_str(),"%*c %ld %ld %ld",
//                       &tail, &head, &cap );
//        edge_descriptor e1, e2;
//        bool in1, in2;
//        boost::tie(e1, in1) = add_edge(verts.at(tail), verts.at(head), g);
//        boost::tie(e2, in2) = add_edge(verts.at(head), verts.at(tail), g);
//
//        if (!in1 || !in2) {
//          std::cerr << "unable to add edge (" << head << "," << tail << ")"
//                    << std::endl;
//          return -1;
//        }
//        capacity[e1] = cap;
//        capacity[e2] = 0;
//        reverse_edge[e1] = e2;
//        reverse_edge[e2] = e1;
    return 0;
}
/* --------------------   end of parser  -------------------*/

  } // namespace detail

template <class Graph, class CapacityMap, class ReverseEdgeMap>
int read_my_dimacs_max_flow(Graph& g,
                         CapacityMap capacity,
                         ReverseEdgeMap reverse_edge,
                         typename graph_traits<Graph>::vertex_descriptor& src,
                         typename graph_traits<Graph>::vertex_descriptor& sink,
                            Bitmap_cubical_complex &cub,
                            size_t birthi, size_t deathi
                         ) {
//    std::filebuf fb;
//    fb.open ("output.txt",std::ios::in);
//    std::istream is (&fb);
    return detail::read_dimacs_max_flow_internal(g, capacity, reverse_edge, src, sink, true, "max", cub, birthi, deathi);
}



} // namespace boost
