///////////////////////////////////////////////////////////////////////////////
//
// THIS SOFTWARE IS PROVIDED "AS-IS". THERE IS NO WARRANTY OF ANY KIND.
// NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE
// FOR ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY.
//
// Copyright (c) 2010 Jyamiti Research Group.
// CS&E Department of the Ohio State University, Columbus, OH.
// All rights reserved.
//
// Author: Sayan Mandal
//
///////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_set>
#include <string>
#include <sys/stat.h>
#include <cstddef>

#ifndef INFINITY
#define INFINITY std::numeric_limits< double >::infinity()
#endif // INFINITY

#include <boost/config.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <boost/geometry.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
//#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>
// standard stuff

#include<ParseCommand.h>
#include <Legal.h>
#include "Graph.h"

using namespace boost;
using namespace std;
namespace bg = boost::geometry;

//extern vector<int> complexSizes;
//extern vector<int> accumulativeSizes;
//extern int accumulativeSize;
//extern double fThreshold;
//extern vector<double> vecFiltrationScale;
//extern SimplicialTree<bool> domain_complex;


typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;

extern map <int,int> threesimplex_to_node; // maps index of tetrahedrons from simplex->first to dual-graph-nodes->second
extern map <int,int> twosimplex_to_tri; // maps index of traingles from simplex->first to dual-graph-edges->second
extern map <int,int> threesimplex_to_node_rev; // maps index of tetrahedrons from dual-graph-nodes->first to simplex->second
extern map <int,int> twosimplex_to_tri_rev; // maps index of traingles from dual-graph-edges->first to simplex->second
extern map<size_t, vector<int>>  sizescale;//0- edge, 1- vertex, for each deathindex, count of vertex and edges
extern vector<double> x;
extern vector<double> y;
extern vector<double> z;
extern map<unsigned long,unsigned long> vertind; // vertex index in cubical complex(first) to actual index in x,y,z, filter


struct cmp_intervals_by_dim_then_length {
    explicit cmp_intervals_by_dim_then_length(Bitmap_cubical_complex * sc)
    : sc_(sc) { }
    template<typename Persistent_interval>
    bool operator()(const Persistent_interval & p1, const Persistent_interval & p2) {
        if (sc_->dimension(get < 0 > (p1)) == sc_->dimension(get < 0 > (p2)))
            return (sc_->filtration(get < 1 > (p1)) - sc_->filtration(get < 0 > (p1))
                    > sc_->filtration(get < 1 > (p2)) - sc_->filtration(get < 0 > (p2)));
        else
            return (sc_->dimension(get < 0 > (p1)) > sc_->dimension(get < 0 > (p2)));
    }
    Bitmap_cubical_complex* sc_;
};

Bitmap_cubical_complex buildComplex(const char* file) {
    
    bool dbg = false;
    std::ifstream inFiltration;
    inFiltration.open(file);
    
    unsigned dimensionOfData;
    inFiltration >> dimensionOfData;
    
    if (dbg) {
        std::cerr << "dimensionOfData : " << dimensionOfData << std::endl;
    }
    
    std::vector<unsigned> sizes;
    sizes.reserve(dimensionOfData);
    // all dimensions multiplied
    std::size_t dimensions = 1;
    size_t top_dim = 1;
    for (std::size_t i = 0; i != dimensionOfData; ++i) {
        unsigned size_in_this_dimension;
        inFiltration >> size_in_this_dimension;
        sizes.push_back(size_in_this_dimension - 1);
        dimensions *= (size_in_this_dimension);
        top_dim *= (size_in_this_dimension - 1);
        if (dbg) {
            std::cerr << "size_in_this_dimension : " << size_in_this_dimension << std::endl;
        }
    }
    if (dbg) {
        std::cerr << "Top dimension names : " << top_dim << std::endl;
    }
    
    //Create a blank Cubical Complex
    vector<double> top_dim_cells(top_dim, 0);
    Bitmap_cubical_complex b(sizes, top_dim_cells);
    double filtrationLevel;
    
    vector<size_t> indices_to_consider;
    
    //Assign values to vertices
    //Bitmap_cubical_complex::Skeleton_simplex_iterator it(&b, 0);
    for (Bitmap_cubical_complex::Skeleton_simplex_iterator it = b.skeleton_simplex_range(0).begin();
         it != b.skeleton_simplex_range(0).end(); it++) {
        if (!(inFiltration >> filtrationLevel) || (inFiltration.eof())) {
            throw std::ios_base::failure("Bad Perseus file format.");
        }
        if (dbg) {
            std::cerr << "Cell of an index : " << (*it)
            << " and dimension: " << b.get_dimension_of_a_cell(*it)
            << " get the value : " << filtrationLevel << std::endl;
        }
        b.get_cell_data(*it) = filtrationLevel;
        size_t s = (*it);
        indices_to_consider.push_back(s);
        
    }
    vector<bool> is_this_cell_considered(b.num_simplices(), false);
    while (indices_to_consider.size()) {
        if (dbg) {
            std::cerr << "indices_to_consider in this iteration \n";
            for (std::size_t i = 0; i != indices_to_consider.size(); ++i) {
                std::cout << indices_to_consider[i] << "  ";
            }
        }
        std::vector<std::size_t> new_indices_to_consider;
        for (std::size_t i = 0; i != indices_to_consider.size(); ++i) {
            std::vector<std::size_t> bd = b.get_coboundary_of_a_cell(indices_to_consider[i]);
            for (std::size_t boundaryIt = 0; boundaryIt != bd.size(); ++boundaryIt) {
                if (dbg) {
                    std::cerr << "filtration of a cell : " << bd[boundaryIt] << " is : " << b.filtration(bd[boundaryIt])
                    << " while of a cell: " << indices_to_consider[i] << " is: " << b.filtration(indices_to_consider[i])
                    << std::endl;
                }
                if (b.filtration(bd[boundaryIt]) < b.filtration(indices_to_consider[i])) {
                    b.get_cell_data(bd[boundaryIt]) = b.filtration(indices_to_consider[i]);
                    if (dbg) {
                        std::cerr << "Setting the value of a cell : " << bd[boundaryIt]
                        << " to : " << b.filtration(indices_to_consider[i]) << std::endl;
                    }
                }
                if (is_this_cell_considered[bd[boundaryIt]] == false) {
                    new_indices_to_consider.push_back(bd[boundaryIt]);
//                    if (b.get_dimension_of_a_cell(bd[boundaryIt]) == dimensionOfData - 1) {
//                        edges->push_back(bd[boundaryIt]);
//                    }
                    is_this_cell_considered[bd[boundaryIt]] = true;
                }
            }
        }
        indices_to_consider.swap(new_indices_to_consider);
    }
    
    
    
    b.initialize_simplex_associated_to_key();
//cub = b;
    return b;
}


int main(int argc, char *argv[] )
{
    int dimensions, noPoints, currentEdgeSize = 0, barcode_count = 0, scalecount = 0;
    int filtration = 0, totinsCount = 0, numbars = 0, index = 0, gap = 0;
  
  string point_file, tr_file, to_file;
  std::vector<bg::model::point<double, 3, bg::cs::cartesian>> vPoint;
  bool verbose;
    ParseCommand(argc, argv,  point_file, numbars, verbose, index, gap);
//  Simplex_tree simplexTree;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"Point_file: "<<point_file; //getchar();
    // THIS LINE GENERATES THE CUBICAL COMPLEX
    Bitmap_cubical_complex cub(point_file.c_str());//buildComplex(point_file.c_str());//(//
//    c = ;
    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh(cub);
    
    int p = 2;
    double min_persistence = 0;
    vector<size_t> birthiv, deathiv;
    pcoh.init_coefficients(p);  // initializes the coefficient field for homology
    pcoh.compute_persistent_cohomology(min_persistence);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ofstream spers(point_file+"_pers");
    ofstream fpers(point_file+"_filtpers");
    
    
    cmp_intervals_by_dim_then_length cmp(&cub);
    auto persistent_pairs = pcoh.get_persistent_pairs();
    cout<<" persistent_pairs: "<<persistent_pairs.size()<<"\n";
    if(verbose)
        getchar();
    sort(std::rbegin(persistent_pairs), std::rend(persistent_pairs), cmp);

    for (auto pair : persistent_pairs) {
        if(cub.dimension(get<0>(pair))!=2)
            continue;
        
        birthiv.push_back(get<0>(pair));
        deathiv.push_back(get<1>(pair));
    }
    
    if(birthiv.size()<numbars)
    {
        cout<<"Not enough 2-cycles: "<<birthiv.size()<<" "<<numbars<<endl;
        getchar();
        exit(0);
    }
//        birthiv.erase(birthiv.begin(),birthiv.end()-numbars);
//        deathiv.erase(deathiv.begin(),deathiv.end()-numbars);
    if(verbose){
        cout<<"birth size: "<<birthiv.size()<<" "<<deathiv.size()<<" numbers: "<<numbars<<" gaps: "<<gap;
    getchar();
    }
    birthiv.erase(birthiv.begin(),birthiv.end()-numbars-gap);
    deathiv.erase(deathiv.begin(),deathiv.end()-numbars-gap);
    birthiv.erase(birthiv.end()-gap,birthiv.end());
    deathiv.erase(deathiv.end()-gap,deathiv.end());
    buildFullCubicSetting(cub, point_file, deathiv);        // GRAPH BUILDING
    cout<<"Sizes: "<<birthiv.size()<<" "<<deathiv.size()<<" "<<numbars<<" "<<(birthiv.size()-numbars);
    if(verbose)
     getchar();

    for(int gotonumbars=birthiv.size()-1; gotonumbars >= 0; gotonumbars--){

        size_t birthi = birthiv[gotonumbars];
        size_t deathi = deathiv[gotonumbars];
//        cout<<"birthi: "<<birthi<<" deathi:"<<deathi<<endl;
//        cout.setf(ios::fixed);
        
        spers<<"2 "<<to_string(int(birthi))<<" "<<to_string(int(deathi))<<"\n";
        fpers<<"2 "<<to_string(float(cub.filtration(birthi)))<<" "<<to_string(float(cub.filtration(deathi)))<<"\n";
        if(verbose){
            cout<<"birthi: "<<birthi<<" deathi:"<<deathi<<endl;
            getchar();
            
        }
        
        
        
        vector<vector<int>> rGraph;
        vector<vector<int>> vectri
        = boostCut( cub, birthi, deathi, verbose);
//        = minCut(adjMatrix, threesimplex_to_node[death], adjMatrix.size()-1, rGraph );
        cout<<" Found Boykov Komolov cut.\n"; //getchar();
        string ind = point_file;

        string dir_path = ind+"loops/";

        boost::filesystem::path dir(dir_path.c_str());
        boost::filesystem::create_directory(dir_path.c_str());
        string ff = ind+"loops/"+to_string(int(birthi))+"_"+to_string(int(deathi))+".off";
        tr_file = ind+"loops/"+to_string(int(birthi))+"_"+to_string(int(deathi))+".txt";
        to_file = ind+"loops/"+to_string(int(birthi))+"_"+to_string(int(deathi))+".off";
        string cu_file = ind+"loops/"+to_string(int(birthi))+"_"+to_string(int(deathi))+"cub.txt";
        cout<<"tr_file:"<<tr_file<<endl;
//        writeToFile( cub, vectri, birthi, deathi, tr_file);//
        writeToFileOFF( cub, vectri, birthi, deathi, to_file);
        cout<<"Written to TEXT: "<<gotonumbars<<" vec: "<<birthiv.size()<<" current index:"<<gotonumbars<<endl; //getchar();

    }
  
    threesimplex_to_node.clear();
    twosimplex_to_tri.clear();
    threesimplex_to_node_rev.clear();
    twosimplex_to_tri_rev.clear();
    sizescale.clear();
    x.clear();
    y.clear();
    z.clear();
    vertind.clear();



  return EXIT_SUCCESS;
}//end of main
 
