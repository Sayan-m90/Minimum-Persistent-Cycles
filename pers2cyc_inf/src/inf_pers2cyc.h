// https://people.sc.fsu.edu/~jburkardt/cpp_src/geometry/geometry.html
// https://www.geometrictools.com/Source/Intersection3D.html

// TODO: use std::vector::shrink_to_fit()

#ifndef _INF_PERS2CYC_H_
#define _INF_PERS2CYC_H_

#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>
// #include "cube_complex.h"
#include "cube_complex_bm.h"
#include "arr_label_id_map.hpp"
#include "graph.hpp"
#include "geom_utils.h"


struct ExecOptions {
    bool verbose;
    bool write_intem;
};


typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> 
    Bitmap_cubical_complex_base;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> 
    Bitmap_cubical_complex;
using Gudhi::persistent_cohomology::Field_Zp;
typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> 
    Persistent_cohomology;


/* definitions for the map from graph edges to 2-cells */

template<typename T1, typename T2>
struct StdPairHash { 
    size_t operator()(const std::pair<T1,T2> &p) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        
        return seed;
    }
};

typedef std::unordered_map<std::pair<int,int>, 
    vector<int>, StdPairHash<int,int>> GEdgeCellMap;


enum VisitColor {
    CL_WHITE, CL_GRAY
};

enum IntervSortType {
    BY_FILT_INDEX, BY_FUNC_VAL
};


void computePersistenceCube(
    const std::string& perseus_fname, 
    const IntervSortType interv_sort_type, 
    const ExecOptions& ops,
    std::string* filt_fname);

void writeFiltCube(
    Bitmap_cubical_complex* bm_cube_cmplx, 
    const int x_res, const int y_res, const int z_res, 
    vector<std::pair<int,int>> pers_pairs_dim2,
    const std::string& filename);

void infPers2CycsFromFile(
    const std::string& filt_fname, 
    const ComplexType cmplx_type, 
    const int start_interval, 
    const int num_intervals, 
    const ExecOptions& ops);

void loadComplex(
    const std::string& filt_fname, 
    const ComplexType cmplx_type, 
    CellComplex* complex, 
    const int interval_id,
    int* interv_start_key,
    int* interv_start_id,
    vector<int>* pos_cell_2d_verts);

// note that the orig_complex is to be deleted in this function
// pos_cell_2d_verts should be in canon order
void infPers2Cyc(
    const ComplexType cmplx_type, CellComplex* orig_complex, 
    const vector<int>& pos_cell_2d_verts, const ExecOptions& ops, 
    const std::string& file_prefix);

void pruneComplex(CellComplex* complex);

void cellConnectedComponent(
    CellComplex* in_complex, const int start_cell_2d_id, CellComplex* out_complex);

void reconVoidBound(
    const ComplexType cmplx_type, 
    CellComplex* complex, 
    vector<array<int,2>>* cell_2d_void_map,
    const ExecOptions& ops, 
    const std::string& file_prefix,
    int* void_cnt);

template<typename WeightType>
void computeMinCyc(
    CellComplex* complex, 
    const GEdgeCellMap* gedge_cell_map, 
    const int void_cnt,
    const int src,
    const int sink, 
    const ExecOptions& ops, 
    const std::string& file_prefix);

void get2DCellPtCubic(
    CellComplex* complex, const int cell_2d_id, 
    const vector<int>& cell_1d_verts, Vector3d* pt);

// get the canonical oriented 2d cell vertex sequence where the orientation is
// from the from_vert to the to_vert. the canonical vertex sequence then starts
// from the smallest vertex following the orientation.
void getCanonOrien2DCellVertsCubic(const vector<int>& verts, 
    const int& from_vert, const int& to_vert, vector<int>* canon_verts);

void addOrien2DCellPairs(
    const ComplexType cmplx_type, 
    CellComplex* complex, 
    const int cell_1d_id,
    const vector<int>& cell_1d_verts, 
    VertArrIdMap* orien_cell_2d_id_map, 
    Graph<int>* orien_cell_2d_graph);

void getBitmapComplexCellVerts(const int cell_handle, const int cell_dim, 
    const int x_res, const int y_res, const int z_res, vector<int>* verts);

void getCell2dVoidIds(
    const int cell_2d_id, CellComplex* complex, 
    vector<array<int,2>>* cell_2d_void_map, vector<int>* void_ids);

inline void toCanonEdge(std::pair<int,int>* e) {
    if (e->first > e->second) {
        std::swap(e->first, e->second);
    }
}

void printGEdgeCellMap(const GEdgeCellMap* gedge_cell_map, CellComplex* complex);

#endif
